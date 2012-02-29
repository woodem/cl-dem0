#ifndef _PARTICLE_CL_H_
#define _PARTICLE_CL_H_

#include"common.cl.h"

CLDEM_NAMESPACE_BEGIN();

struct Sphere{	Real radius; };
inline struct Sphere Sphere_new(){ struct Sphere ret; ret.radius=NAN; return ret; }

enum _shape_enum { Shape_None=0, Shape_Sphere, };


struct Particle;
static void Particle_init(global struct Particle*);

// http://gpu.doxos.eu/trac/wiki/OpenCLDataStructures
struct Particle{
	int flags;
	Vec3 pos, vel, angVel, bboxPos;
	Quat ori;
	Vec3 inertia;
	Real mass;
	Vec3 force, torque;
	union {
		struct Sphere sphere;
	} AMD_UNION_ALIGN_BUG_WORKAROUND() shape;
	int mutex;
	// needed for boost::python
	#ifdef __cplusplus
		Particle(){ Particle_init(this); }
	#endif
};

// inline struct Particle Particle_new(){ struct Particle p; Particle_init(&p); return p; }

#define PAR_LEN_shapeT  2
#define PAR_LEN_clumped 1
#define PAR_LEN_stateT  2
#define PAR_LEN_dofs    6
#define PAR_LEN_groups  4
#define PAR_LEN_matId   6

#define PAR_OFF_shapeT  0
#define PAR_OFF_clumped PAR_OFF_shapeT + PAR_LEN_shapeT
#define PAR_OFF_stateT  PAR_OFF_clumped + PAR_LEN_clumped
#define PAR_OFF_dofs    PAR_OFF_stateT + PAR_LEN_stateT
#define PAR_OFF_groups  PAR_OFF_dofs + PAR_LEN_dofs
#define PAR_OFF_matId   PAR_OFF_groups + PAR_LEN_groups

constant flagSpec par_flag_shapeT ={PAR_OFF_shapeT,PAR_LEN_shapeT};
constant flagSpec par_flag_clumped={PAR_OFF_clumped,PAR_LEN_clumped};
constant flagSpec par_flag_stateT ={PAR_OFF_stateT,PAR_LEN_stateT};
constant flagSpec par_flag_dofs   ={PAR_OFF_dofs,PAR_LEN_dofs};
constant flagSpec par_flag_groups ={PAR_OFF_groups,PAR_LEN_groups};
constant flagSpec par_flag_matId  ={PAR_OFF_matId,PAR_LEN_matId};
/* END */

// static_assert(par_flag_groups.x+par_flag_groups.y<=32); // don't overflow int

#define PARTICLE_FLAG_GET_SET(what) \
	inline int par_##what##_get(global const struct Particle *p){ return flags_get(p->flags,par_flag_##what); } \
	inline void par_##what##_set(global struct Particle *p, int val){ flags_set(&(p->flags),par_flag_##what,val); } \
	inline void par_##what##_set_local(struct Particle *p, int val){ flags_set_local(&(p->flags),par_flag_##what,val); }
PARTICLE_FLAG_GET_SET(shapeT);
PARTICLE_FLAG_GET_SET(clumped);
PARTICLE_FLAG_GET_SET(stateT);
PARTICLE_FLAG_GET_SET(dofs);
PARTICLE_FLAG_GET_SET(groups);
PARTICLE_FLAG_GET_SET(matId);

inline int dof_axis(int axis, int rot){ return 1<<(axis+(rot?3:0)); }
constant int par_dofs_trans=7; // 0b0000111
constant int par_dofs_rot=56;  // 0b0111000
constant int par_dofs_all=63;  // 0b0111111

static void Particle_init(global struct Particle* p){
	p->flags=0;
	p->pos=Vec3_set(NAN,NAN,NAN);
	p->bboxPos=Vec3_set(NAN,NAN,NAN);
	p->ori=Quat_identity();
	p->inertia=Vec3_set(1,1,1);
	p->mass=1.;
	p->vel=p->angVel=Vec3_set(0,0,0);
	p->mutex=0;
	p->force=p->torque=Vec3_set(0,0,0);

	par_shapeT_set(p,0);
	par_clumped_set(p,0);
	par_stateT_set(p,0);
	par_dofs_set(p,par_dofs_all);
	par_groups_set(p,0);
	par_matId_set(p,0);
}


CLDEM_NAMESPACE_END();

#ifdef __cplusplus
namespace clDem{
	// needed for py::indexing_suite
		inline bool operator==(const Particle& a, const Particle& b){ return memcmp(&a,&b,sizeof(Particle))==0; }	

	static
	py::object Particle_shape_get(Particle* p){
		int shapeT=par_shapeT_get(p);
		switch(shapeT){
			case 0: return py::object(); // None
			case Shape_Sphere: return py::object(p->shape.sphere);
			default:	throw std::runtime_error("Particle has shape with unknown index "+lexical_cast<string>(shapeT));
		}
	}
	static
	void Particle_shape_set(Particle *p, py::object sh){
		if(sh==py::object()){ par_shapeT_set(p,0); return; }
		py::extract<Sphere> sphere(sh);
		if(sphere.check()){
			p->shape.sphere=sphere();
			par_shapeT_set(p,Shape_Sphere);
		}
		else throw std::runtime_error("Unknown shape object.");
	}
	// not needed inside yade
	#ifndef YADE_CLDEM
		static
		void Particle_cl_h_expose(){
			py::class_<Sphere>("Sphere").def("__init__",Sphere_new).PY_RW(Sphere,radius);

			py::class_<Particle>("Particle")//.def("__init__",py::make_constructor(Particle_new))
				.PY_RWV(Particle,pos).PY_RWV(Particle,ori).PY_RWV(Particle,inertia).PY_RW(Particle,mass).PY_RWV(Particle,vel).PY_RWV(Particle,force).PY_RWV(Particle,torque).PY_RWV(Particle,bboxPos)
				.add_property("shape",Particle_shape_get,Particle_shape_set)
				// flags
				.def_readonly("flags",&Particle::flags)
				.add_property("shapeT",par_shapeT_get)
				.add_property("clumped",par_clumped_get,par_clumped_set)
				.add_property("stateT",par_stateT_get)
				.add_property("dofs",par_dofs_get,par_dofs_set)
				.add_property("groups",par_groups_get,par_groups_set)
				.add_property("matId",par_matId_get,par_matId_set)
			;

			// from-python as list
			custom_vector_from_seq<Particle>();
			// to-python as ParticleList proxy
			//py::class_<std::vector<Particle>>("ParticleList").def(py::vector_indexing_suite<std::vector<Particle>>());
			py::class_<std::vector<Particle>>("ParticleList").def(py::vector_indexing_suite<std::vector<Particle>>());
		}
	#endif

};

#endif

#endif
