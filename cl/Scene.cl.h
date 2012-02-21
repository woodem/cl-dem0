#ifndef _SCENE_CL_H_
#define _SCENE_CL_H_

#include"common.cl.h"

CLDEM_NAMESPACE_BEGIN();

/* possibly we won't need different materials in one single simulation; stay general for now */
struct ElastMat{ Real density, young; };
struct ElastMat ElastMat_new(){ struct ElastMat ret; ret.density=NAN; ret.young=NAN; return ret; }

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
	int mat_##what##_get(global const struct Material *m){ return flags_get(m->flags,mat_flag_##what); } \
	void mat_##what##_set(global struct Material *m, int val){ flags_set(&m->flags,mat_flag_##what,val); } \
	void mat_##what##_set_local(struct Material *m, int val){ flags_set_local(&m->flags,mat_flag_##what,val); }
MATERIAL_FLAG_GET_SET(matT);

#define MATT2_COMBINE(m1,m2) ((m1) | (m2)<<(MAT_LEN_matT))


enum _energy{ENERGY_Ekt=0,ENERGY_Ekr,ENERGY_grav,ENERGY_damp,ENERGY_elast,SCENE_ENERGY_NUM };
struct  EnergyProperties {
	const char name[16];
	long int index; /* ! must be long, AMD otherwise reads garbage for incremental?! !! */ 
	cl_bool incremental;
};
constant struct EnergyProperties energyDefinitions[]={
	{"Ekt",  ENERGY_Ekt,  false},
	{"Ekr",  ENERGY_Ekt,  false},
	{"grav", ENERGY_grav, true }, 
	{"damp", ENERGY_damp, true },
	{"elast",ENERGY_elast,false},
};

// substep numbers for kernels (they must be always run in this order)
enum _substeps{ SUB_nextTimestep=0, SUB_integrator, SUB_updateBboxes, SUB_contCompute, SUB_forcesToParticles, };
// interrupt codes
enum _interrupts{ INT_BBOXES_UPDATED=0, };

#define SCENE_MAT_NUM 8
struct Scene{
	Real t;
	Real dt;
	long step;
	struct Interrupt {
		int step; // when -1, no interrupt; // FIXME: should be long, but there are no atomics on longs!
		int substep;
		int what;
		cl_bool destructive;
	} interrupt;
	Vec3 gravity;
	Real damping;
	Real verletDist;
	cl_bool updateBboxes;
	struct Material materials[SCENE_MAT_NUM];
	float energy[SCENE_ENERGY_NUM];
	int energyMutex[SCENE_ENERGY_NUM]; // use mutexes until atomics work for floats
};


struct Scene Scene_new(){
	struct Scene s;
	s.t=0;
	s.dt=1e-8;
	s.step=-1;
	s.gravity=Vec3_set(0,0,0);
	s.damping=0.;
	s.verletDist=0.;
	s.interrupt.step=-1;
	s.interrupt.substep=s.interrupt.what=s.interrupt.destructive=0;
	s.updateBboxes=false;
	// no materials
	for(int i=0; i<SCENE_MAT_NUM; i++) mat_matT_set_local(&s.materials[i],0); 
	//Scene_energyReset(&s): but the pointer is not global, copy here instead
	for(int i=0; i<SCENE_ENERGY_NUM; i++){ s.energy[i]=0; s.energyMutex[i]=0; }
	return s;
}

#ifdef __OPENCL_VERSION__
bool Scene_interrupt_set(global struct Scene* s, int substep, int what, bool destructive){
	if(atom_cmpxchg(&s->interrupt.step,-1,s->step)==-1){
		printf("%d/%d: Setting interrupt %d\n",/*make AMD happy*/(int)s->step,substep,what);
		s->interrupt.substep=substep;
		s->interrupt.what=what;
		s->interrupt.destructive=destructive;
		return true;
	}
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
				// destructive interrupt: other workitems with the current kernel can be skipped
				(s->interrupt.destructive && s->interrupt.substep<=substep)
				||
				// non-destructive: other workitems within the same kernel must complete
				(!s->interrupt.destructive && s->interrupt.substep<substep)
			)
		)
	);
}
#endif

// exposed to python, hence must be visible to the host
void Scene_energyReset(global struct Scene* s){
	for(int i=0; i<SCENE_ENERGY_NUM; i++) s->energy[i]=0;
}

CLDEM_NAMESPACE_END();

#ifdef __cplusplus
	py::object Material_mat_get(const Material* m){
		int matT=mat_matT_get(m);
		switch(matT){
			case 0: return py::object();
			case Mat_ElastMat: return py::object(m->mat.elast);
			default: throw std::runtime_error("Material has mat with unknown index "+lexical_cast<string>(matT));
		}
	}

	void Material_mat_set(Material *m, py::object mat){
		if(mat==py::object()){ mat_matT_set(m,0); return; }
		py::extract<ElastMat> elast(mat);
		if(elast.check()){ m->mat.elast=elast(); mat_matT_set(m,Mat_ElastMat); }
		else throw std::runtime_error("Unknown mat object.");
	}

	py::list Scene_mats_get(Scene* self){
		py::list ret;
		for(const Material& m: self->materials) ret.append(Material_mat_get(&m));
		return ret;
	}

	void Scene_mats_set(Scene* self, const py::list mm){
		if(py::len(mm)>=SCENE_MAT_NUM) throw std::runtime_error("Too many materials, to be defined (maximum "+lexical_cast<string>(SCENE_MAT_NUM)+", given "+lexical_cast<string>(py::len(mm)));
		for(int i=0; i<py::len(mm); i++){ Material_mat_set(&(self->materials[i]),mm[i]); }
	};

	py::dict Scene_energy_get(Scene* self/*, bool omitZero=true*/){
		const bool omitZero=true;
		py::dict ret;
		for(int i=0; i<SCENE_ENERGY_NUM; i++){
			if(self->energy[i]!=0. || !omitZero) ret[energyDefinitions[i].name]=self->energy[i];
		}
		return ret;
	}

	Real Scene_energyTotal(Scene* s){
		Real ret=0;
		for(int i=0; i<SCENE_ENERGY_NUM; i++) ret+=s->energy[i];
		return ret;
	}

	Real Scene_energyError(Scene* s){
		Real sum=0,absSum=0;
		for(int i=0; i<SCENE_ENERGY_NUM; i++){ sum+=s->energy[i]; absSum+=std::abs(s->energy[i]); }
		return sum/absSum;  // return NaN for absSum==0
	}

	void Scene_cl_h_expose(){
		VECTOR_SEQ_CONV(Material);
		py::class_<ElastMat>("ElastMat").def("__init__",ElastMat_new).PY_RW(ElastMat,density).PY_RW(ElastMat,young);

		py::class_<Scene>("Scene").def("__init__",Scene_new)
			.PY_RW(Scene,t).PY_RW(Scene,dt).PY_RW(Scene,step).PY_RWV(Scene,gravity).PY_RW(Scene,damping)
			.add_property("materials",Scene_mats_get,Scene_mats_set)
			.add_property("energy",Scene_energy_get) //,py::arg("omitZero")=true)
			.def("energyReset",Scene_energyReset)
			.def("energyTotal",Scene_energyTotal)
			.def("energyError",Scene_energyError)
		;
	}
#endif

#endif
