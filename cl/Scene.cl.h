#ifndef _SCENE_CL_H_
#define _SCENE_CL_H_

#include"common.cl.h"
#include"kernels.cl.h"

CLDEM_NAMESPACE_BEGIN();

/* possibly we won't need different materials in one single simulation; stay general for now */
struct ElastMat{ Real density, young; };
inline struct ElastMat ElastMat_new(){ struct ElastMat ret; ret.density=NAN; ret.young=NAN; return ret; }

enum _mat_enum { Mat_None=0, Mat_ElastMat=1, };

struct Material{
	int flags;
	union{
		struct ElastMat elast;
	}  AMD_UNION_ALIGN_BUG_WORKAROUND() mat;
};

#define MAT_LEN_matT 3
#define MAT_OFF_matT 0
constant flagSpec mat_flag_matT={MAT_OFF_matT,MAT_LEN_matT};
#define MATERIAL_FLAG_GET_SET(what) \
	inline int mat_##what##_get_global(global const struct Material *m){ return flags_get(m->flags,mat_flag_##what); } \
	inline int mat_##what##_get(const struct Material *m){ return flags_get(m->flags,mat_flag_##what); } \
	inline void mat_##what##_set_global(global struct Material *m, int val){ flags_set_global(&m->flags,mat_flag_##what,val); } \
	inline void mat_##what##_set(struct Material *m, int val){ flags_set_local(&m->flags,mat_flag_##what,val); }
MATERIAL_FLAG_GET_SET(matT);

#define MATT2_COMBINE(m1,m2) ((m1) | (m2)<<(MAT_LEN_matT))


enum _energy{ENERGY_Ekt=0,ENERGY_Ekr,ENERGY_grav,ENERGY_damp,ENERGY_elast,SCENE_ENERGY_NUM };
struct  EnergyProperties {
	const char name[16];
	long int index; /* ! must be long, AMD otherwise reads garbage for incremental?! !! */ 
	cl_bool incremental;
};
static constant struct EnergyProperties energyDefinitions[]={
	{"Ekt",  ENERGY_Ekt,  false},
	{"Ekr",  ENERGY_Ekt,  false},
	{"grav", ENERGY_grav, true }, 
	{"damp", ENERGY_damp, true },
	{"elast",ENERGY_elast,false},
};

// interrupt flags
enum _int_flags {
	// null flags, for clarity of the code
	INT_NOT_IMMEDIATE=0, INT_NOT_DESTRUCTIVE=0, 
	INT_IMMEDIATE=1,
	INT_DESTRUCTIVE=2,
	INT_ARRAYS=4,
	INT_BBOXES=8,
};
// dynamic arrays indices
enum _dynarrays { ARR_CON=0, ARR_CONFREE, ARR_POT, ARR_POTFREE, ARR_CJOURNAL, ARR_CLUMPS, ARR_NUM };

#ifdef __cplusplus
static string arrName(int arrIx){
	switch(arrIx){
		case ARR_CON: return "CON";
		case ARR_CONFREE: return "CONFREE";
		case ARR_POT: return "POT";
		case ARR_POTFREE: return "POTFREE";
		case ARR_CJOURNAL: return "CJOURNAL";
		case ARR_CLUMPS: return "CLUMPS";
		default: throw std::logic_error("arrName("+lexical_cast<string>(arrIx)+") not known.");
	};
};
#endif

#define SCENE_MAT_NUM 8
constant int SCENE_MAT_NUM_=SCENE_MAT_NUM;

struct Scene{
	Real t;
	Real dt;
	long step;
	cl_int rollback; // used to detect very first step at adjust arrays, not used inside kernels
	struct Interrupt {
		cl_int step; // when -1, no interrupt; // FIXME: should be long, but there are no atomics on longs!
		cl_int substep;
		cl_int flags;
	} interrupt;
	Vec3 gravity;
	Real damping;
	Real verletDist;
	/*
	groupmask for groups which don't interact mutually;
	(p1.group & p2.group & scene.loneGroups)!=0 ⇒ no contact between p1 & p2
	*/
	cl_int loneGroups;
	cl_bool updateBboxes;
	struct Material materials[SCENE_MAT_NUM];
	float energy[SCENE_ENERGY_NUM];
	int energyMutex[SCENE_ENERGY_NUM]; // use mutexes until atomics work for floats
	cl_int arrSize[ARR_NUM];  // must be ints since we need atomics for them
	cl_int arrAlloc[ARR_NUM];
	#ifdef __cplusplus
		Scene();
	#endif
};


#ifdef __cplusplus

inline void Scene_interrupt_reset(struct Scene* s){
	s->interrupt.step=-1;
	s->interrupt.substep=s->interrupt.flags=0;
}

inline void Scene_init(struct Scene* s){
	s->t=0;
	s->dt=1e-8;
	s->step=-1;
	s->gravity=Vec3_set(0,0,0);
	s->damping=0.;
	s->verletDist=-1.; // no interrupts due to spheres getting out of bboxes
	s->loneGroups=0;
	Scene_interrupt_reset(s);
	s->updateBboxes=false;
	s->rollback=0;
	// no materials
	for(int i=0; i<SCENE_MAT_NUM; i++) mat_matT_set(&s->materials[i],0); 
	//Scene_energyReset(&s): but the pointer is not global, copy here instead
	for(int i=0; i<SCENE_ENERGY_NUM; i++){ s->energy[i]=0; s->energyMutex[i]=0; }
	// no allocated arrays
	for(int i=0; i<ARR_NUM; i++) s->arrAlloc[i]=s->arrSize[i]=0;
}
inline Scene::Scene(){ Scene_init(this); }

static py::dict Scene_arr_get(struct Scene* s){
	py::dict ret;
	ret["con"]=py::make_tuple(s->arrSize[ARR_CON],s->arrAlloc[ARR_CON]);
	ret["conFree"]=py::make_tuple(s->arrSize[ARR_CONFREE],s->arrAlloc[ARR_CONFREE]);
	ret["pot"]=py::make_tuple(s->arrSize[ARR_POT],s->arrAlloc[ARR_POT]);
	ret["potFree"]=py::make_tuple(s->arrSize[ARR_POTFREE],s->arrAlloc[ARR_POTFREE]);
	ret["cJournal"]=py::make_tuple(s->arrSize[ARR_CJOURNAL],s->arrAlloc[ARR_CJOURNAL]);
	ret["clumps"]=py::make_tuple(s->arrSize[ARR_CLUMPS],s->arrAlloc[ARR_CLUMPS]);
	return ret;
}


#endif

/*
Decide whether particles may create contact at all:
1. (not yet implemented:) clumps can't collide (only their member particles can)
2. particles which have a common group which is in Scene::loneGroups may not collide
3. particles which share no common group may not collide

This functions is called by the collider; non-interacting particles should not
be inserted in potential contacts at all.
 */
inline bool Scene_particles_may_collide(global struct Scene* s, global struct Particle* p1, global struct Particle* p2){
	// check if one of the particles is a clump, in that case there is no collision possible
	if((par_groups_get_global(p1) & par_groups_get_global(p2) & s->loneGroups)!=0) return false;
	if((par_groups_get_global(p1) & par_groups_get_global(p2))==0) return false;
	return true;
}

/* try to allocate an additional item in array with index arrIx;
if the array would overflow, returns -1, but the arrSize will
nevertheless be increased (the arrSize thus designates required array
size). The new item should be written at the index returned, unless
negative.
*/
inline long Scene_arr_append(global struct Scene* s, int arrIx){
	long oldSize=atom_inc(&(s->arrSize[arrIx]));
	if(oldSize>=s->arrAlloc[arrIx]) return -1; // array size not sufficient
	return oldSize;
};

// in addition, try to claim the new element by flipping it to 0 
inline long Scene_arr_append_claim(global struct Scene* s, int arrIx, global int* arr){
	while(true){
		long oldSize=atom_inc(&(s->arrSize[arrIx]));
		if(/* new size*/oldSize+1 > s->arrAlloc[arrIx]) return -1; // array size not sufficient
		// flip the element so that it is not unclaimed;
		// if that fails, another thread was faster and we need to try again
		//printf("Trying %d at %ld\n",arr[oldSize],oldSize);
		if(atom_cmpxchg(&arr[oldSize],-1,0)==-1){
			//printf("Returning: %d at %ld\n",arr[oldSize],oldSize);
			// for(int i=0; i<s->arrAlloc[arrIx]; i++){ printf("->%d %d: %d\n",oldSize,i,arr[arrIx]); }
			return oldSize;
		}
	}
};
/*
traverse an array of ints from the beginning;
change the first negative value atomically to 0 and return its index;
if no negative value is found, append it to the array.

This function should be used with potFree and conFree; it muse be tried
whether it is better for the performance to write predictably at the end
of those arrays (trading speed for size) or going through them first
(trading size for speed)
*/
inline long Scene_arr_findNegative_or_append(global struct Scene* s, int arrIx, global int* arr){
	// disable this block to prefer speed at the expense of size
	#if 1
		for(int i=0; i<s->arrSize[arrIx]; i++){
			// the element is negative and gets atomically changed to 0
			#ifdef __cplusplus
				//cerr<<"At element "<<i<<"="<<arr[i]<<endl;
			#endif
			if(atom_cmpxchg(&arr[i],-1,0)==-1){
				#ifdef __cplusplus
					//cerr<<"Changed to 0, returning "<<i<<endl;
				#endif
				return i;
			}
		}
	#endif
	return Scene_arr_append_claim(s,arrIx,arr);
}

/*
return free index in array with index arrIx,
trying first non-negative elements of arrFree (which has the index arrFreeIx);
if no free slot is found, try to enlarge arr; returns -1 if the reallocation fails;
in that case, the caller is responsible for setting an appropriate interrupt.
If *shrunk* is true and the last element of arrFree is used, the array is traversed
backwards and shrunk so that there are no trailing invalid values.

NB: the shrinking logic does not work well with concurrent access and is disabled now
 */
inline long Scene_arr_fromFreeArr_or_append(global struct Scene* scene, global int* arrFree, int arrFreeIx, int arrIx, bool shrink){
	int ix=-1;
	for(int i=0; i<scene->arrSize[arrFreeIx]; i++){
		//printf("Trying arrFree[%d]=%d: ",i,arrFree[i]);

		// optimize: don't xchg values which are already negative
		if(arrFree[i]<0) continue;

		ix=atom_xchg(&(arrFree[i]),-1); // read old value and reset to -1
		//printf("%d\n",ix);
		if(ix>=0){ return ix; }
	}
	//printf("Allocating new slot, arrSize %d\n.",scene->arrSize[arrIx]);
	// no free slot found, allocate a new one; handle possible allocation failure
	ix=Scene_arr_append(scene,arrIx);
	// -1 if appending failed, otherwise a valid index
	// return in both cases
	return ix;
}

#ifdef __OPENCL_VERSION__
bool Scene_interrupt_set(global struct Scene* s, int substep, int flags){
	if(atom_cmpxchg(&s->interrupt.step,-1,s->step)==-1){
		// enabling this printf makes some bboxes NaN?!!!
		// printf("%d/%d: Setting interrupt %d\n",/*make AMD happy*/(int)s->step,substep,what);
		s->interrupt.substep=substep;
		s->interrupt.flags=flags;
		return true;
	}
	// this is in itself OK, since another interrupt was set already;
	// use for debugging now
	//printf("%d/%d: Interrupt %d not set!\n",(int)s->step,substep,what);
	return false;
}

bool Scene_skipKernel(global struct Scene* s, int substep){
	return
		s->interrupt.step>=0 /* valid interrupt */
		&& (
			s->interrupt.step<s->step /* interrupted in previous step already */
			||
			(s->interrupt.step==s->step /* interrupted in this step */
			&& (
				// non-destructive: other workitems within the same kernel must complete
				(s->interrupt.substep<substep)
				||
				// immediate interrupt: other workitems with the current kernel can be skipped
				(s->interrupt.substep==substep && (s->interrupt.flags&INT_IMMEDIATE))
			)
		)
	);
}

#endif

// exposed to python, hence must be visible to the host
inline void Scene_energyReset(global struct Scene* s){
	for(int i=0; i<SCENE_ENERGY_NUM; i++) s->energy[i]=0;
}

CLDEM_NAMESPACE_END();

#ifdef __cplusplus
namespace clDem{
	static
	py::object Material_mat_get(const Material* m){
		int matT=mat_matT_get(m);
		switch(matT){
			case 0: return py::object();
			case Mat_ElastMat: return py::object(m->mat.elast);
			default: throw std::runtime_error("Material has mat with unknown index "+lexical_cast<string>(matT));
		}
	}
	static
	void Material_mat_set(Material *m, py::object mat){
		if(mat==py::object()){ mat_matT_set(m,0); return; }
		py::extract<ElastMat> elast(mat);
		if(elast.check()){ m->mat.elast=elast(); mat_matT_set(m,Mat_ElastMat); }
		else throw std::runtime_error("Unknown mat object.");
	}
	static
	py::list Scene_mats_get(Scene* self){
		py::list ret;
		for(const Material& m: self->materials) ret.append(Material_mat_get(&m));
		return ret;
	}
	static
	void Scene_mats_set(Scene* self, const py::list mm){
		if(py::len(mm)>=SCENE_MAT_NUM) throw std::runtime_error("Too many materials, to be defined (maximum "+lexical_cast<string>(SCENE_MAT_NUM)+", given "+lexical_cast<string>(py::len(mm)));
		for(int i=0; i<py::len(mm); i++){ Material_mat_set(&(self->materials[i]),mm[i]); }
	};
	static
	py::dict Scene_energy_get(Scene* self/*, bool omitZero=true*/){
		const bool omitZero=true;
		py::dict ret;
		for(int i=0; i<SCENE_ENERGY_NUM; i++){
			if(self->energy[i]!=0. || !omitZero) ret[energyDefinitions[i].name]=self->energy[i];
		}
		return ret;
	}
	static
	Real Scene_energyTotal(Scene* s){
		Real ret=0;
		for(int i=0; i<SCENE_ENERGY_NUM; i++) ret+=s->energy[i];
		return ret;
	}
	static
	Real Scene_energyError(Scene* s){
		Real sum=0,absSum=0;
		for(int i=0; i<SCENE_ENERGY_NUM; i++){ sum+=s->energy[i]; absSum+=std::abs(s->energy[i]); }
		return sum/absSum;  // return NaN for absSum==0
	}
	static
	py::dict Scene_interrupt_get(Scene *s){
		py::dict ret;
		ret["step"]=s->interrupt.step;
		ret["substep"]=s->interrupt.substep;
		ret["destructive"]=(bool)(s->interrupt.flags&INT_DESTRUCTIVE);
		ret["immediate"]=(bool)(s->interrupt.flags&INT_IMMEDIATE);
		ret["arrays"]=(bool)(s->interrupt.flags&INT_ARRAYS);
		ret["bboxes"]=(bool)(s->interrupt.flags&INT_BBOXES);
		return ret;
	}
	#ifndef YADE_CLDEM
		static
		void Scene_cl_h_expose(){
			VECTOR_SEQ_CONV(Material);
			py::class_<ElastMat>("ElastMat").def("__init__",ElastMat_new).PY_RW(ElastMat,density).PY_RW(ElastMat,young);

			py::class_<Scene>("Scene")
				.PY_RW(Scene,t).PY_RW(Scene,dt).PY_RW(Scene,step).PY_RWV(Scene,gravity).PY_RW(Scene,damping).PY_RW(Scene,verletDist).PY_RW(Scene,loneGroups)
				.add_property("arr",Scene_arr_get)
				.add_property("materials",Scene_mats_get,Scene_mats_set)
				.add_property("energy",Scene_energy_get) //,py::arg("omitZero")=true)
				.add_property("interrupt",Scene_interrupt_get)
				.def("energyReset",Scene_energyReset)
				.def("energyTotal",Scene_energyTotal)
				.def("energyError",Scene_energyError)
			;
		}
	#endif
};
#endif

#endif
