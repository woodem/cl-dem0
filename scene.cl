/* vim:set syntax=c: */

#include"../cl-math/basic-math.cl"

typedef short2 flagSpec; // offset and length, in bits
typedef long par_id_t;
typedef long2 par_id2_t;

int flags_get(int flags, flagSpec spec){ return (flags>>spec.s0)&((1<<spec.s1)-1); }
void flags_set(int flags, int val, flagSpec spec){ flags|=(val&((1<<spec.s1)-1))<<spec.s0; }

struct Sphere{	Real radius; };
Sphere Sphere_new(){ Sphere ret; ret.radius=nan(); return ret; }

enum { Shape_Sphere=1, };

// http://gpu.doxos.eu/trac/wiki/OpenCLDataStructures
struct Particle{
	int flags;
	Vec3 pos, vel, angVel;
	Quat ori;
	Vec3 inertia;
	Real mass;
	int semaphore; // HACK!
	Vec3 force, torque;
	union{
		Sphere sphere;
	} shape;
};

const flagSpec par_flag_shapeT =(0,2);
const flagSpec par_flag_clumped=(par_flag_shapeT.x+par_flag_shapeT.y,1);
const flagSpec par_flag_stateT =(par_flag_clumped.x+par_flag.clumped.y,2);
const flagSpec par_flag_dofs   =(par_flag_stateT.x+par_flag_stateT.t,6);
const flagSpec par_flag_groups =(par_flag_dofs.x+par_flag_dofx.y,4);
const flagSpec par_flag_matId  =(par_flag_groups.x+par_flag_groups.y,8);
assert(par_flag.groups.x+par_flag_groups.y<=32); // don't overflow int

#define PARTICLE_FLAG_GET_SET(what) \
	int par_##what##_get(Particle *p){ return flags_get(p->flags,par_flag_##what); } \
	void par_##what##_set(Particle *p, int what){ flags_set(p->flags,par_flag_##what,what); }
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


struct L1Geom{ Real uN; };
L1Geom L1Geom_new(){ L1Geom ret; ret.uN=nan(); return L1Geom; }

enum { Geom_L1Geom=1, };

struct NormPhys{ Real kN; }
NormPhys NormPhys_new{ NormPhys ret; ret.kN=nan(); return NormPhys; }

enum { Phys_NormPhys=1, };

struct Contact{
	int flags;
	par_id2_t ids;
	Vec3 pos;
	Quat ori;
	Vec3 force, torque;
	union {
		L1Geom l1g;
	} geom;
	union {
		NormPhys normPhys;
	} phys;
};

const flagSpec con_flag_shapesT=(0,2*par_flag_shapeT.y);
const flagSpec con_flag_geomT  =(con_flag_shapesT.x+con_flag_shapesT.y,2);
const flagSpec con_flag_physT  =(con_flag_geomT.x+con_flag_geomT.y,2);
#define PARTICLE_FLAG_GET_SET(what) \
	int con_##what##_get(Contact *c){ return flags_get(p->flags,con_flag_##what); } \
	void con_##what##_set(Contact *c, int what){ flags_set(p->flags,con_flag_##what,what); }
PARTICLE_FLAG_GET_SET(shapesT);
PARTICLE_FLAG_GET_SET(geomT);
PARTICLE_FLAG_GET_SET(physT);

int shapeT_combine2(int sh1, int sh2){ return sh1 & sh2<<(par_flag_shapeT.y); }
int geomT_physT_combine2(int g, int p){ return g & p<<(con_flag_physT.y); }


/* possibly we won't need different materials in one single simulation; stay general for now */
struct ElastMat{ Real young; }
ElastMat ElastMat_new(){ ElastMat ret; ret.young=nan(); }

enum { Mat_ElastMat=1, };

struct Material{
	int flags;
	union{
		ElastMat elast;
	} mat;
};

const flagSpec mat_flag_matT=(0,2);
#define MATERIAL_FLAG_GET_SET(what) \
	int mat_##what##_get(Material *m){ return flags_get(m->flags,mat_flag_##what); } \
	void mat_##what##_set(Material *m, int what){ flags_set(m->flags,mat_flag_##what,what); }
MATERIAL_FLAG_GET_SET(matT);

int matT_combine2(int m1, int m2){ return m1 & m2<<(mat_flag_matT.y); }


struct Scene{
	Real t;
	long step;
	Real dt;
	Vec3 gravity;
	Real damping;
	Material[8] materials;
};


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
		if(_sem_old!=0){ atomic_dec(&(sem)); continue }
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
		PSEUDO_SEM_CRITICAL_BEGIN(p->sem);
			p->force+=(i==0?Fp1:Fp2);
			p->torque+=(i==0?Tp1:Tp2);
		PSEUDO_SEM_CRITICAL_END(p->sem);
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
	switch(con_shapesT_get(c)){
		case shapeT_combine2(par_shapeT_get(p1),par_shapeT_get(p2)):
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
				c->phys.normPhys.kN=.5*(m1.mat.young+m2.mat.young);
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

