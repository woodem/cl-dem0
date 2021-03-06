#ifndef _PARTICLE_CL_H_
#define _PARTICLE_CL_H_

#include"common.cl.h"
#include"serialization.cl.h"

CLDEM_NAMESPACE_BEGIN();

struct Sphere{
	Real radius;
	#ifdef GCC46
		Sphere(): radius(NAN) {}
	#endif
	CLDEM_SERIALIZE_ATTRS((radius),/* */)
};

struct Wall{
	cl_short axis;
	cl_short sense;
	#ifdef GCC46
		Wall(): axis(-1), sense(0){}
	#endif
	CLDEM_SERIALIZE_ATTRS((axis)(sense),/* */)
};

struct Clump{
	cl_long ix; // index in the clumps array
	#ifdef GCC46
		Clump(): ix(-1){};
	#endif
	CLDEM_SERIALIZE_ATTRS((ix),/* */)
};

// stored in the clumps array
struct ClumpMember{
	par_id_t id;
	Vec3 relPos;
	Quat relOri;
	#ifdef GCC46
		ClumpMember(): id(-1){}
	#endif
	CLDEM_SERIALIZE_ATTRS((id)(relPos)(relOri),/* */)
};

// Sphere should come at the end
enum _shape_enum { Shape_None=0, Shape_Clump, Shape_Wall, Shape_Sphere };

struct Particle;
static void Particle_init(struct Particle*);
int par_shapeT_get(const struct Particle* p);

// http://gpu.doxos.eu/trac/wiki/OpenCLDataStructures
struct Particle{
	int flags;
	// negative for particles which are not within clump at all
	par_id_t clumpId; 
	Vec3 pos, vel, angVel, bboxPos;
	Quat ori;
	Vec3 inertia, angMom; // angular momentum for aspherical particles
	Real mass;
	Vec3 force, torque;
	union _shape {
		UNION_BITWISE_CTORS(_shape)
		struct Sphere sphere;
		struct Wall wall;
		struct Clump clump;
	} AMD_UNION_ALIGN_BUG_WORKAROUND() shape;
	int mutex;
	CLDEM_SERIALIZE_ATTRS((flags)(pos)(vel)(angVel)(bboxPos)(ori)(inertia)(mass)(force)(torque),
		switch(par_shapeT_get(this)){
			case Shape_None: break;
			case Shape_Sphere: ar & boost::serialization::make_nvp("sphere",shape.sphere); break;
			case Shape_Wall: ar & boost::serialization::make_nvp("wall",shape.wall); break;
			case Shape_Clump: ar & boost::serialization::make_nvp("clump",shape.clump); break;
			default: throw std::runtime_error("Invalid shapeT value at (de)serialization.");
		};
	)
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
	inline int par_##what##_get_local(const struct Particle *p) { return flags_get(p->flags, par_flag_##what); } \
	inline int par_##what##_get_global(global const struct Particle *p) { return flags_get(p->flags, par_flag_##what); } \
	inline int par_##what##_get(const struct Particle *p){ return flags_get(p->flags,par_flag_##what); } \
	inline void par_##what##_set_global(global struct Particle *p, int val){ flags_set_global(&(p->flags),par_flag_##what,val); } \
	inline void par_##what##_set(struct Particle *p, int val){ flags_set_local(&(p->flags),par_flag_##what,val); }
PARTICLE_FLAG_GET_SET(shapeT);
PARTICLE_FLAG_GET_SET(clumped);
PARTICLE_FLAG_GET_SET(stateT);
PARTICLE_FLAG_GET_SET(dofs);
PARTICLE_FLAG_GET_SET(groups);
PARTICLE_FLAG_GET_SET(matId);

inline int dof_axis(int axis, bool rot){ return 1<<(axis+(rot?3:0)); }
constant int par_dofs_trans=7; // 0b0000111
constant int par_dofs_rot=56;  // 0b0111000
constant int par_dofs_all=63;  // 0b0111111

// all groups
constant int par_groups_all=(1<<PAR_LEN_groups)-1;

static void Particle_init(struct Particle* p){
	p->flags=0;
	p->clumpId=-1;
	p->pos=Vec3_set(NAN,NAN,NAN);
	p->bboxPos=Vec3_set(NAN,NAN,NAN);
	p->ori=Quat_identity();
	p->inertia=Vec3_set(1,1,1);
	p->angMom=Vec3_set(0,0,0);
	p->mass=1.;
	p->vel=p->angVel=Vec3_set(0,0,0);
	p->mutex=0;
	p->force=p->torque=Vec3_set(0,0,0);

	par_shapeT_set(p,0);
	par_clumped_set(p,0);
	par_stateT_set(p,0);
	par_dofs_set(p,par_dofs_all);
	par_groups_set(p,par_groups_all);
	par_matId_set(p,0);
}

#ifdef __OPENCL_VERSION__

void Clump_collectFromMembers(struct Particle* C, global struct ClumpMember* clumps, global struct Particle* par){
	// go through the clumps array from given index until an invalid item is found
	for(size_t i=C->shape.clump.ix; clumps[i].id>=0; i++){
		struct Particle cp=par[clumps[i].id];
		C->force+=cp.force;
		C->torque+=cp.torque+cross(cp.pos-C->pos,cp.force);
	}
}

void Clump_applyToMembers(struct Particle* C, global struct ClumpMember* clumps, global struct Particle* par, global bool* updateBboxes, Real verletDist){
	for(size_t i=C->shape.clump.ix; clumps[i].id>=0; i++){
		// work on the local copy instead?
		const struct ClumpMember cm=clumps[i];
		global struct Particle* cp=&(par[cm.id]);
		cp->pos=C->pos+Quat_rotate(C->ori,cm.relPos);
		cp->ori=Quat_multQ(C->ori,cm.relOri);
		cp->vel=C->vel+cross(C->angVel,cp->pos-C->pos);
		cp->angVel=C->angVel;
		cp->force=cp->torque=(Vec3)0.;
		if(/*non-NULL*/ updateBboxes && !isnan(verletDist) && (isnan(cp->bboxPos.s0) ||  Vec3_sqNorm(cp->pos-cp->bboxPos)>pown(verletDist,2))) *updateBboxes=true;
	}
}

#endif


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
			case Shape_Wall: return py::object(p->shape.wall);
			default:	throw std::runtime_error("Particle has shape with unknown index "+lexical_cast<string>(shapeT));
		}
	}
	static
	void Particle_shape_set(Particle *p, py::object sh){
		if(sh==py::object()){ par_shapeT_set(p,0); return; }
		py::extract<Sphere> sphere(sh);
		py::extract<Wall> wall(sh);
		if(sphere.check()){
			p->shape.sphere=sphere();
			par_shapeT_set(p,Shape_Sphere);
		}
		else if(wall.check()){
			p->shape.wall=wall();
			par_shapeT_set(p,Shape_Wall);
		}
		else throw std::runtime_error("Unknown shape object.");
	}
	// not needed inside woo
	#ifndef WOO_CLDEM
		static
		void Particle_cl_h_expose(){
			py::class_<Sphere>("Sphere").PY_RW(Sphere,radius);
			py::class_<Wall>("Wall").PY_RW(Wall,axis).PY_RW(Wall,sense);
			py::class_<Clump>("Clump").add_property("ix",&Clump::ix);

			py::class_<ClumpMember>("ClumpMember").add_property("id",&ClumpMember::id).PY_RWV(ClumpMember,relPos).PY_RWV(ClumpMember,relOri);
			VECTOR_SEQ_CONV(ClumpMember);

			py::class_<Particle>("Particle")//.def("__init__",py::make_constructor(Particle_new))
				.PY_RWV(Particle,pos).PY_RWV(Particle,ori).PY_RWV(Particle,inertia).PY_RWV(Particle,angMom).PY_RW(Particle,mass).PY_RWV(Particle,vel).PY_RWV(Particle,angVel).PY_RWV(Particle,force).PY_RWV(Particle,torque).PY_RWV(Particle,bboxPos)
				.add_property("shape",Particle_shape_get,Particle_shape_set)
				// flags
				.def_readonly("mutex", &Particle::mutex)
				.def_readonly("flags",&Particle::flags)
				.add_property("shapeT",par_shapeT_get)
				.add_property("clumped",par_clumped_get,par_clumped_set)
				.add_property("stateT",par_stateT_get)
				.add_property("dofs",par_dofs_get,par_dofs_set)
				.add_property("groups",par_groups_get,par_groups_set)
				.add_property("matId",par_matId_get,par_matId_set)
			;

			// from-python as list
			_custom_vector_from_seq<Particle>();
			// to-python as ParticleList proxy
			py::class_<std::vector<Particle>>("ParticleList").def(py::vector_indexing_suite<std::vector<Particle>>());
		}
	#endif

};

#endif

#endif
