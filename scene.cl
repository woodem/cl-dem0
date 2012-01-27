/* vim:set syntax=c: */
#ifndef _SCENE_CL_
#define _SCENE_CL_
#include<assert.h>
#include"../cl-math/basic-math.cl"

#ifdef __OPENCL_VERSION__
	#define cl_short2 short2
	#define cl_long2 long2
#endif

#include<stdbool.h>

typedef cl_short2 flagSpec; // offset and length, in bits
typedef long par_id_t;
typedef cl_long2 par_id2_t;

int flags_get(int flags, flagSpec spec){ return (flags>>spec.x)&((1<<spec.y)-1); }
void flags_set(int flags, flagSpec spec, int val){ flags|=(val&((1<<spec.y)-1))<<spec.x; }

struct Sphere{	Real radius; };
struct Sphere Sphere_new(){ struct Sphere ret; ret.radius=NAN; return ret; }

enum { Shape_Sphere=1, };

// http://gpu.doxos.eu/trac/wiki/OpenCLDataStructures
struct Particle{
	int flags;
	Vec3 pos, vel, angVel;
	Quat ori;
	Vec3 inertia;
	Real mass;
	int semaphore;
	Vec3 force, torque;
	union{
		struct Sphere sphere;
	} shape;
};
#define PAR_LEN_shapeT  0
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

const flagSpec par_flag_shapeT ={PAR_OFF_shapeT,PAR_LEN_shapeT};
const flagSpec par_flag_clumped={PAR_OFF_clumped,PAR_LEN_clumped};
const flagSpec par_flag_stateT ={PAR_OFF_stateT,PAR_LEN_stateT};
const flagSpec par_flag_dofs   ={PAR_OFF_dofs,PAR_LEN_dofs};
const flagSpec par_flag_groups ={PAR_OFF_groups,PAR_LEN_groups};
const flagSpec par_flag_matId  ={PAR_OFF_matId,PAR_LEN_matId};
/* END */

// static_assert(par_flag_groups.x+par_flag_groups.y<=32); // don't overflow int

#define PARTICLE_FLAG_GET_SET(what) \
	int par_##what##_get(struct Particle *p){ return flags_get(p->flags,par_flag_##what); } \
	void par_##what##_set(struct Particle *p, int val){ flags_set(p->flags,par_flag_##what,val); }
PARTICLE_FLAG_GET_SET(shapeT);
PARTICLE_FLAG_GET_SET(clumped);
PARTICLE_FLAG_GET_SET(stateT);
PARTICLE_FLAG_GET_SET(dofs);
PARTICLE_FLAG_GET_SET(groups);
PARTICLE_FLAG_GET_SET(matId);

int dof_axis(int axis, bool rot){ return 1<<(axis+(rot?3:0)); }
const int par_dofs_trans=7; // 0b0000111
const int par_dofs_rot=56;  // 0b0111000
const int par_dofs_all=63;  // 0b0111111


struct Particle Particle_new(){
	struct Particle p;
	p.flags=0;
	p.pos={NAN,NAN,NAN};
	p.inertia={1,1,1};
	p.mass=1.
	p.vel=p.angVel={0,0,0};
	p.semaphore=0;
	p.force=p.torque={0,0,0};

	par_shapeT_set(&p,0);
	par_clumped_set(&p,0);
	par_stateT_set(&p,0);
	par_dofs_set(&p,par_dofs_all);
	par_groups_set(&p,0);
	par_matId_set(&p,0);
}



struct L1Geom{ Real uN; };
struct L1Geom L1Geom_new(){ struct L1Geom ret; ret.uN=NAN; return ret; }

enum { Geom_L1Geom=1, };

struct NormPhys{ Real kN; };
struct NormPhys NormPhys_new(){ struct NormPhys ret; ret.kN=NAN; return ret; }

enum { Phys_NormPhys=1, };

struct Contact{
	int flags;
	par_id2_t ids;
	Vec3 pos;
	Quat ori;
	Vec3 force, torque;
	union {
		struct L1Geom l1g;
	} geom;
	union {
		struct NormPhys normPhys;
	} phys;
};

#define CON_LEN_shapesT 2*(PAR_LEN_shapeT) // currently unused in the actual code
#define CON_LEN_geomT   2
#define CON_LEN_physT   2

#define CON_OFF_shapesT 0
#define CON_OFF_geomT   CON_OFF_shapesT + CON_LEN_shapesT
#define CON_OFF_physT   CON_OFF_geomT + CON_LEN_geomT

const flagSpec con_flag_shapesT={CON_OFF_shapesT,CON_LEN_shapesT};
const flagSpec con_flag_geomT  ={CON_OFF_geomT,CON_LEN_geomT};
const flagSpec con_flag_physT  ={CON_OFF_physT,CON_LEN_physT};
#define CONTACT_FLAG_GET_SET(what) \
	int con_##what##_get(struct Contact *c){ return flags_get(c->flags,con_flag_##what); } \
	void con_##what##_set(struct Contact *c, int val){ flags_set(c->flags,con_flag_##what,val); }
CONTACT_FLAG_GET_SET(shapesT);
CONTACT_FLAG_GET_SET(geomT);
CONTACT_FLAG_GET_SET(physT);

struct Contact Contact_new(){
	struct Contact c;
	c.pos={0,0,0};
	c.ori={1,0,0,0};
	c.force=c.torque={0,0,0};
	con_shapesT_set(&c,0);
	con_geomT_set(&c,0);
	con_physT_set(&c,0);
	return c;
}

int shapeT_combine2(int sh1, int sh2){ return sh1 & sh2<<(par_flag_shapeT.y); }
int geomT_physT_combine2(int g, int p){ return g & p<<(con_flag_physT.y); }


/* possibly we won't need different materials in one single simulation; stay general for now */
struct ElastMat{ Real young; };
struct ElastMat ElastMat_new(){ struct ElastMat ret; ret.young=NAN; return ret; }

enum { Mat_ElastMat=1, };

struct Material{
	int flags;
	union{
		struct ElastMat elast;
	} mat;
};

const flagSpec mat_flag_matT={0,2};
#define MATERIAL_FLAG_GET_SET(what) \
	int mat_##what##_get(struct Material *m){ return flags_get(m->flags,mat_flag_##what); } \
	void mat_##what##_set(struct Material *m, int val){ flags_set(m->flags,mat_flag_##what,val); }
MATERIAL_FLAG_GET_SET(matT);

int matT_combine2(int m1, int m2){ return m1 & m2<<(mat_flag_matT.y); }


struct Scene{
	Real t;
	long step;
	Real dt;
	Vec3 gravity;
	Real damping;
	struct Material materials[8];
};

struct Scene Scene_new(){
	struct Scene s;
	s.t=0;
	s.dt=1e-8;
	s.step=0;
	s.gravity={0,0,0};
	s.damping=0.;
}

#ifdef __OPENCL_VERSION__

kernel void nextTimestep(global Scene* scene);
kernel void forcesToParticles(global const Scene* scene, global Particle* par, global const Contact* con);
kernel void integrator (global const Scene* scene, const Particles* par);
kernel void contCompute(global const Scene*, global const Partices* par, global const Contact* con);

kernel void nextTimestep(global Scene* scene){
	scene->t+=scene->dt;
}

#define PSEUDO_SEM_CRITICAL_BEGIN(sem) \
	while(true){ \
		while((sem)!=0); \
		int _sem_old=atomic_inc(&(sem)); \
		/* semaphore taken between while and atomic_inc, release and try again */ \
		if(atomic_inc(&(sem))!=1){ atomic_dec(&(sem)); continue }
#define PSEUDO_SEM_CRITICAL_END(sem) \
		atomic_dec(&(sem)); /* release the semaphore*/ \
		break; /* end the infinite loop */ \
	}


kernel void forcesToParticles(global Scene* scene, global Particle* par, global const Contact* con){
	if(scene->step==0) return;
	/* how to synchronize access to particles? */
	Contact* con=scene->con[get_global_id()];
	Particle* p1=scene->par[con->ids.x], p2=scene->par[con->ids.y];
	Mat3 R_T=Quat_toMat3(Quat_conjugate(con->ori));
	Vec3 Fp1=+Mat3_multV(R_T,con->force);
	Vec3 Fp2=-Mat3_multV(R_T,con->force);
	Vec3 Tp1=+cross(ori->pos-p1->pos,Mat3_multV(R_T,con->force))+Mat3_multV(R_T,con->torque);
	Vec3 Tp2=-cross(ori->pos-p2->pos,Mat3_multV(R_T,con->force))-Mat3_multV(R_T,con->torque);
	// write to particles, "protect" via semaphore
	for(int i=0; i<2; i++){
		p=(i==0?p1:p2);
		PSEUDO_SEM_CRITICAL_BEGIN(p->semaphore);
			p->force+=(i==0?Fp1:Fp2);
			p->torque+=(i==0?Tp1:Tp2);
		PSEUDO_SEM_CRITICAL_END(p->semaphore);
	}
}

kernel void integrator(global const Scene* scene, global Particle* par){
	if(scene->step==0) return;
	Particle* p=scene->par[get_global_id()];
	int dofs=par_dofs_get(p);
	Vec3 accel=0., angAccel=0.;
	// optimize for particles with all translations/rotations (vector ops)
	// aspherical integration (if p->inertia has the same components)
	for(int ax=0; ax<3; ax++){
		if(dof_axis(dofs,ax,false)) (*Real(&accel))[ax]+=p->force/p->mass+(*Real(&scene->gravity))[ax];
		if(dof_axis(dofs,ax,true )) (*Real(&angAccel))[ax]+=p->torque/p->inertia.x;
	}
	if(scene->damping!=0){
		accel   =accel   *(1-damping*sgn(p->force *(p->vel   +accel   *scene->dt/2)));
		angAccel=angAccel*(1-damping*sgn(p->torque*(p->angVel+angAccel*scene->dt/2)));
	}
	p->vel+=accel*scene->dt;
	p->angVel+=angAccel*scene->dt;
	p->pos+=p->vel*scene->dt;
	p->ori=Quat_multQ(Quat_fromRotVec(p->angVel*scene->dt),p->ori);
}

kernel void contCompute(globa const Scene*, global const Particles* par, global const Contact* con){
	Contact* c=scene->con[get_global_id()];
	Particle* p1=scene->par[c->ids.s0];
	Particle* p2=scene->par[c->ids.s1];
	/* create geometry for new contacts */
	if(con_geomT_get(c)==0){
		con_geomT_set(c,Geom_L1Geom); c->geom.l1g=L1Geom_new();
	}
	/* update geometry */
	switch(shapeT_combine2(par_shapeT_get(p1),par_shapeT_get(p2))){
		case shapeT_combine2(Shape_Sphere,Shape_Sphere):
			assert(con_geomT_get(c)==Geom_L1Geom);
			const L1G* l1g=c->geom.l1g;
			l1g.uN=distance(p1->pos,p2->pos)-p1->shape.sphere.radius-p2->shape.sphere.radius;
		default: /* signal error */
	}
	/* update physical params: only if there are no physical params yet */
	if(con_physT_get(c)==0){
		Material* m1=scene->materials[p1->matId];
		Material* m2=scene->materials[p2->matId];
		switch(matT_combine2(m1->matT,m2->matT)){
			case matT_combine(Mat_ElastMat,Mat_ElastMat):
				c->phys.normPhys=NormPhys_new();
				// ratio should depend on particle radius (TODO)
				c->phys.normPhys.kN=.5*(m1.mat.young+m2.mat.young) /* multiply by cross-section times length etc */;
			default: /* error */
		}
	}
	/* contact law */
	switch(geomT_physT_combine2(c->geomT,c->physT)){
		case geomT_physT_combine2(Geom_L1Geom,Phys_NormPhys):
			c->force=c->geom->l1g.uN*c->phys.normPhys.kN;
			c->torque=(Vec3)0.
		default: /* error */
	}
}

#endif /* __OPENCL_VERSION__ */
#endif /* _SCENE_CL_ */
