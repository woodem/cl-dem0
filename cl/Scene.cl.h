#ifndef _SCENE_CL_H_
#define _SCENE_CL_H_

#include"common.cl.h"
#include"kernels.cl.h"

CLDEM_NAMESPACE_BEGIN();

/* possibly we won't need different materials in one single simulation; stay general for now */
struct ElastMat{ Real density, young; };
inline struct ElastMat ElastMat_new(){ struct ElastMat ret; ret.density=NAN; ret.young=NAN; return ret; }

enum _mat_enum { Mat_ElastMat=1, };

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
	inline int mat_##what##_get(global const struct Material *m){ return flags_get(m->flags,mat_flag_##what); } \
	inline void mat_##what##_set(global struct Material *m, int val){ flags_set(&m->flags,mat_flag_##what,val); } \
	inline void mat_##what##_set_local(struct Material *m, int val){ flags_set_local(&m->flags,mat_flag_##what,val); }
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

// interrupt codes
enum _interrupts{ INT_BBOXES_UPDATED=0, INT_ARR_CON, INT_ARR_CONFREE, INT_ARR_POT, INT_ARR_POTFREE, INT_NUM };
// interrupt flags
enum _int_flags {
	INT_NOT_IMMEDIATE=0, INT_NOT_DESTRUCTIVE=0, // null flags, for clarity of code
	INT_IMMEDIATE=1, INT_DESTRUCTIVE=2 };
// dynamic arrays indices
enum _dynarrays { ARR_CON=0, ARR_CONFREE, ARR_POT, ARR_POTFREE, ARR_NUM };

#define SCENE_MAT_NUM 8
struct Scene{
	Real t;
	Real dt;
	long step;
	struct Interrupt {
		cl_int step; // when -1, no interrupt; // FIXME: should be long, but there are no atomics on longs!
		cl_int substep;
		cl_int what;
		cl_int flags;
	} interrupt;
	Vec3 gravity;
	Real damping;
	Real verletDist;
	cl_bool updateBboxes;
	struct Material materials[SCENE_MAT_NUM];
	float energy[SCENE_ENERGY_NUM];
	int energyMutex[SCENE_ENERGY_NUM]; // use mutexes until atomics work for floats
	cl_int arrSize[ARR_NUM];  // must be ints since we need atomics for them
	cl_int arrAlloc[ARR_NUM];
};


#ifdef __cplusplus
inline void Scene_interrupt_reset(struct Scene* s){
	s->interrupt.step=-1;
	s->interrupt.substep=s->interrupt.what=s->interrupt.flags=0;
}
inline struct Scene Scene_new(){
	struct Scene s;
	s.t=0;
	s.dt=1e-8;
	s.step=-1;
	s.gravity=Vec3_set(0,0,0);
	s.damping=0.;
	s.verletDist=-1.; // no interrupts due to spheres getting out of bboxes
	Scene_interrupt_reset(&s);
	s.updateBboxes=false;
	// no materials
	for(int i=0; i<SCENE_MAT_NUM; i++) mat_matT_set_local(&s.materials[i],0); 
	//Scene_energyReset(&s): but the pointer is not global, copy here instead
	for(int i=0; i<SCENE_ENERGY_NUM; i++){ s.energy[i]=0; s.energyMutex[i]=0; }
	// no allocated arrays
	for(int i=0; i<ARR_NUM; i++) s.arrAlloc[i]=s.arrSize[i]=0;
	return s;
}
#endif


#ifdef __OPENCL_VERSION__
/* try to allocate an additional item in array with index arrIx;
if the array would overflow, returns -1, but the arrSize will
nevertheless be increased (the arrSize thus designates required array
size). The new item should be written at the index returned, unless
negative. */
long Scene_arr_append(global struct Scene* s, int arrIx){
	long oldSize=atom_inc(&(s->arrSize[arrIx]));
	if(oldSize>=s->arrAlloc[arrIx]) return -1; // array size not sufficient
	return oldSize;
};

/*
return free index in array with index arrIx,
trying first non-negative elements of arrFree (which has the index arrFreeIx);
if no free slot is found, try to enlarge arr; returns -1 if the reallocation fails;
in that case, the caller is responsible for setting an appropriate interrupt.
If *shrunk* is true and the last element of arrFree is used, the array is traversed
backwards and shrunk so that there are no trailing invalid values.
 */
long Scene_arr_free_or_append(global struct Scene* scene, global int* arrFree, int arrFreeIx, int arrIx, bool shrink){
	int ix=-1;
	for(int i=0; i<scene->arrSize[arrFreeIx]; i++){
		//printf("Trying arrFree[%d]=%d: ",i,arrFree[i]);
		ix=atom_xchg(&(arrFree[i]),-1); // read old value and reset to -1
		//printf("%d\n",ix);
		if(ix>=0){
			// in case we just grabbed the last element in arrFree, shrink its size (go to the first valid item)
			if(shrink && ix==scene->arrSize[arrFreeIx]-1){
				int last; for(last=ix-1; arrFree[last]<0 && last>=0; last--);
				//printf("Array shrunk to %d\n",last+1);
				scene->arrSize[arrFreeIx]=last+1;
			}
			return ix;
		}
	}
	//printf("Allocating new slot, arrSize %d\n.",scene->arrSize[arrIx]);
	// no free slot found, allocate a new one; handle possible allocation failure
	ix=Scene_arr_append(scene,arrIx);
	// -1 if appending failed, otherwise a valid index
	// return in both cases
	return ix;
}

bool Scene_interrupt_set(global struct Scene* s, int substep, int what, int flags){
	if(atom_cmpxchg(&s->interrupt.step,-1,s->step)==-1){
		// enabling this printf makes some bboxes NaN?!!!
		// printf("%d/%d: Setting interrupt %d\n",/*make AMD happy*/(int)s->step,substep,what);
		s->interrupt.substep=substep;
		s->interrupt.what=what;
		s->interrupt.flags=flags;
		return true;
	}
	// this is in itself OK, since another interrupt was set already;
	// use for debugging now
	printf("%d/%d: Interrupt %d not set!\n",(int)s->step,substep,what);
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
		ret["what"]=s->interrupt.what;
		ret["destructive"]=(bool)(s->interrupt.flags&INT_DESTRUCTIVE);
		ret["immediate"]=(bool)(s->interrupt.flags&INT_IMMEDIATE);
		return ret;
	}
	static
	void Scene_cl_h_expose(){
		VECTOR_SEQ_CONV(Material);
		py::class_<ElastMat>("ElastMat").def("__init__",ElastMat_new).PY_RW(ElastMat,density).PY_RW(ElastMat,young);

		py::class_<Scene>("Scene").def("__init__",Scene_new)
			.PY_RW(Scene,t).PY_RW(Scene,dt).PY_RW(Scene,step).PY_RWV(Scene,gravity).PY_RW(Scene,damping).PY_RW(Scene,verletDist)
			.add_property("materials",Scene_mats_get,Scene_mats_set)
			.add_property("energy",Scene_energy_get) //,py::arg("omitZero")=true)
			.add_property("interrupt",Scene_interrupt_get)
			.def("energyReset",Scene_energyReset)
			.def("energyTotal",Scene_energyTotal)
			.def("energyError",Scene_energyError)
		;
	}
#endif

#endif