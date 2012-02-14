/* vim:set syntax=c: */
#ifndef _SCENE_CL_
#define _SCENE_CL_

#include"../cl-math/basic-math.cl"


// AMD does not align unions correctly, help it a bit here
// the 128 is the biggest alignment (on cl_double16), which is very
// wasteful
// optionally, alignment could be specified on each union
// to accomodate the biggest alignment of contained structs
#define AMD_UNION_ALIGN_BUG_WORKAROUND() __attribute__((aligned(128)))


#ifdef __OPENCL_VERSION__
	#define cl_short2 short2
	#define cl_long2 long2
	// printf in OpenCL code
	#ifdef cl_intel_printf
		#pragma OPENCL EXTENSION cl_intel_printf: enable
	#elif defined(cl_amd_printf)
		#pragma OPENCL EXTENSION cl_amd_printf: enable
	#else
		#define printf(...)
	#endif
#else
	#define global
	#define constant const
#endif

//#ifdef cl_khr_global_int32_base_atomics
	#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
//#else
//	#error cl_khr_global_int32_base_atomics extension not supported
//#endif

#ifdef __cplusplus
	#include<cstring>
	#define UNION_BITWISE_CTORS(a) a(){}; a(const a& _b){ memcpy(this,&_b,sizeof(a));}; a& operator=(const a& _b){ memcpy(this,&_b,sizeof(a)); return *this; };
#else
	#define UNION_BITWISE_CTORS(a)
#endif



typedef cl_short2 flagSpec; // offset and length, in bits
typedef long par_id_t;
typedef cl_long2 par_id2_t;

int flags_get(const int flags, const flagSpec spec){ return (flags>>spec.x)&((1<<spec.y)-1); }
void flags_set(global int *flags, const flagSpec spec, int val){
	(*flags)&=~(((1<<spec.y)-1)<<spec.x); /* zero field */
	val&=((1<<spec.y)-1); /* zero excess bits */
	(*flags)|=val<<spec.x;  /* set field */
}
void flags_set_local(int *flags, const flagSpec spec, int val){
	(*flags)&=~(((1<<spec.y)-1)<<spec.x); /* zero field */
	val&=((1<<spec.y)-1); /* zero excess bits */
	(*flags)|=val<<spec.x;  /* set field */
}

struct Sphere{	Real radius; };
struct Sphere Sphere_new(){ struct Sphere ret; ret.radius=NAN; return ret; }

enum _shape_enum { Shape_Sphere=1, };

// http://gpu.doxos.eu/trac/wiki/OpenCLDataStructures
struct Particle{
	int flags;
	Vec3 pos, vel, angVel;
	Quat ori;
	Vec3 inertia;
	Real mass;
	Vec3 force, torque;
	union {
		struct Sphere sphere;
	} AMD_UNION_ALIGN_BUG_WORKAROUND() shape;
	int mutex;
};
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
	int par_##what##_get(global const struct Particle *p){ return flags_get(p->flags,par_flag_##what); } \
	void par_##what##_set(global struct Particle *p, int val){ flags_set(&(p->flags),par_flag_##what,val); } \
	void par_##what##_set_local(struct Particle *p, int val){ flags_set_local(&(p->flags),par_flag_##what,val); }
PARTICLE_FLAG_GET_SET(shapeT);
PARTICLE_FLAG_GET_SET(clumped);
PARTICLE_FLAG_GET_SET(stateT);
PARTICLE_FLAG_GET_SET(dofs);
PARTICLE_FLAG_GET_SET(groups);
PARTICLE_FLAG_GET_SET(matId);

int dof_axis(int axis, int rot){ return 1<<(axis+(rot?3:0)); }
constant int par_dofs_trans=7; // 0b0000111
constant int par_dofs_rot=56;  // 0b0111000
constant int par_dofs_all=63;  // 0b0111111


struct Particle Particle_new(){
	struct Particle p;
	p.flags=0;
	p.pos=Vec3_set(NAN,NAN,NAN);
	p.ori=Quat_identity();
	p.inertia=Vec3_set(1,1,1);
	p.mass=1.;
	p.vel=p.angVel=Vec3_set(0,0,0);
	p.mutex=0;
	p.force=p.torque=Vec3_set(0,0,0);

	par_shapeT_set_local(&p,0);
	par_clumped_set_local(&p,0);
	par_stateT_set_local(&p,0);
	par_dofs_set_local(&p,par_dofs_all);
	par_groups_set_local(&p,0);
	par_matId_set_local(&p,0);

	return p;
}



struct L1Geom{ Real uN; };
struct L1Geom L1Geom_new(){ struct L1Geom ret; ret.uN=NAN; return ret; }

struct L6Geom{ Real uN; Vec3 vel; Vec3 angVel; };
struct L6Geom L6Geom_new(){ struct L6Geom ret; ret.uN=NAN; ret.vel=ret.angVel=Vec3_set(0.,0.,0.); return ret; }

enum _geom_enum { Geom_L1Geom=1, Geom_L6Geom };

struct NormPhys{ Real kN; };
struct NormPhys NormPhys_new(){ struct NormPhys ret; ret.kN=NAN; return ret; }

enum _phys_enum { Phys_NormPhys=1, };

struct Contact{
	int flags;
	par_id2_t ids;
	Vec3 pos;
	Mat3 ori;
	Vec3 force, torque;
	union _gg {
		UNION_BITWISE_CTORS(_gg)
		struct L1Geom l1g;
		struct L6Geom l6g;
	}  AMD_UNION_ALIGN_BUG_WORKAROUND() geom;
	union {
		struct NormPhys normPhys;
	}  AMD_UNION_ALIGN_BUG_WORKAROUND() phys;
};

#define CON_LEN_shapesT 2*(PAR_LEN_shapeT) // currently unused in the actual code
#define CON_LEN_geomT   2
#define CON_LEN_physT   2

#define CON_OFF_shapesT 0
#define CON_OFF_geomT   CON_OFF_shapesT + CON_LEN_shapesT
#define CON_OFF_physT   CON_OFF_geomT + CON_LEN_geomT

constant flagSpec con_flag_shapesT={CON_OFF_shapesT,CON_LEN_shapesT};
constant flagSpec con_flag_geomT  ={CON_OFF_geomT,CON_LEN_geomT};
constant flagSpec con_flag_physT  ={CON_OFF_physT,CON_LEN_physT};
#define CONTACT_FLAG_GET_SET(what) \
	int con_##what##_get(global const struct Contact *c){ return flags_get(c->flags,con_flag_##what); } \
	void con_##what##_set(global struct Contact *c, int val){ flags_set(&c->flags,con_flag_##what,val); } \
	int con_##what##_get_local(const struct Contact *c){ return flags_get(c->flags,con_flag_##what); } \
	void con_##what##_set_local(struct Contact *c, int val){ flags_set_local(&c->flags,con_flag_##what,val); } 
CONTACT_FLAG_GET_SET(shapesT);
CONTACT_FLAG_GET_SET(geomT);
CONTACT_FLAG_GET_SET(physT);

struct Contact Contact_new(){
	struct Contact c;
	c.pos=Vec3_set(0,0,0);
	c.ori=Mat3_identity();
	c.force=c.torque=Vec3_set(0,0,0);
	con_shapesT_set_local(&c,0);
	con_geomT_set_local(&c,0);
	con_physT_set_local(&c,0);
	return c;
}

//int shapeT_combine2(int sh1, int sh2){ return sh1 & sh2<<(par_flag_shapeT.y); }
//int geomT_physT_combine2(int g, int p){ return g & p<<(con_flag_physT.y); }
#define SHAPET2_COMBINE(s1,s2) ((s1) | (s2)<<(PAR_LEN_shapeT))
#define GEOMT_PHYST_COMBINE(g,p) ((g) | (p)<<(CON_LEN_geomT))


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

// int matT_combine2(int m1, int m2){ return m1 & m2<<(mat_flag_matT.y); }
#define MATT2_COMBINE(m1,m2) ((m1) | (m2)<<(MAT_LEN_matT))

#define SCENE_NUM_MAT 8

struct Scene{
	Real t;
	Real dt;
	long step;
	Vec3 gravity;
	Real damping;
	struct Material materials[SCENE_NUM_MAT];
};

struct Scene Scene_new(){
	struct Scene s;
	s.t=0;
	s.dt=1e-8;
	s.step=-1;
	s.gravity=Vec3_set(0,0,0);
	s.damping=0.;
	// no materials
	for(int i=0; i<SCENE_NUM_MAT; i++) mat_matT_set_local(&s.materials[i],0); 
	return s;
}

#ifdef __OPENCL_VERSION__

kernel void nextTimestep(global struct Scene* scene);
kernel void forcesToParticles(global const struct Scene* scene, global struct Particle* par, global const struct Contact* con);
kernel void integrator (global struct Scene* scene, global struct Particle* par);
kernel void contCompute(global const struct Scene*, global const struct Particle* par, global struct Contact* con);

kernel void nextTimestep(global struct Scene* scene){
	// if step is -1, we are at the very beginning
	// keep time at 0, only increase step number
	if(scene->step>0) scene->t+=scene->dt;
	scene->step++;
}

#define TRYLOCK(a) atom_cmpxchg(a,0,1)
#define LOCK(a) while(TRYLOCK(a))
#define UNLOCK(a) atom_xchg(a,0)

kernel void forcesToParticles(global const struct Scene* scene, global struct Particle* par, global const struct Contact* con){
	//if(scene->step==0) return;
	/* how to synchronize access to particles? */
	global const struct Contact* c=&(con[get_global_id(0)]);
	global struct Particle *p1=&(par[c->ids.x]), *p2=&(par[c->ids.y]);
	Mat3 R_T=c->ori;
	Vec3 Fp1=+Mat3_multV(R_T,c->force);
	Vec3 Fp2=-Mat3_multV(R_T,c->force);
	Vec3 Tp1=+cross(c->pos-p1->pos,Mat3_multV(R_T,c->force))+Mat3_multV(R_T,c->torque);
	Vec3 Tp2=-cross(c->pos-p2->pos,Mat3_multV(R_T,c->force))-Mat3_multV(R_T,c->torque);
	LOCK(&p1->mutex);
		p1->force+=Fp1; p1->torque+=Tp1;
	UNLOCK(&p1->mutex);
	LOCK(&p2->mutex);
		p2->force+=Fp2; p2->torque+=Tp2;
	UNLOCK(&p2->mutex);
}

kernel void integrator(global struct Scene* scene, global struct Particle* par){
	if(scene->step==0) return;
	global struct Particle *p=&(par[get_global_id(0)]);
	int dofs=par_dofs_get(p);
	Vec3 accel=0., angAccel=0.;
	// optimize for particles with all translations/rotations (vector ops)
	// aspherical integration (if p->inertia has the same components)
#if 0
	for(int ax=0; ax<3; ax++){
		if(dofs & dof_axis(ax,false)) ((Real*)(&accel))[ax]+=((Real*)(&p->force))[ax]/p->mass+((Real*)(&scene->gravity))[ax]*p->mass;
		if(dofs & dof_axis(ax,true ))	((Real*)(&angAccel))[ax]+=((Real*)(&p->torque))[ax]/p->inertia.x;
	}
#else
	// use vector ops for particles free in all 3 dofs, and on-eby-one only for those which are partially fixed
	p->force+=scene->gravity*p->mass; // put this up here so that grav appears in the force term
	if((dofs&par_dofs_trans)==par_dofs_trans){ // all translations possible, use vector expr
		accel+=p->force/p->mass;
	} else {
		for(int ax=0; ax<3; ax++){
			if(dofs & dof_axis(ax,false)) ((Real*)(&accel))[ax]+=((global Real*)(&(p->force)))[ax]/p->mass+((global Real*)(&(scene->gravity)))[ax];
			else ((Real*)(&accel))[ax]=0.;
		}
	}
	if((dofs&par_dofs_rot)==par_dofs_rot){
		angAccel+=p->torque/p->inertia.x;
	} else {
		for(int ax=0; ax<3; ax++){
			if(dofs & dof_axis(ax,true )) ((Real*)(&angAccel))[ax]+=((global Real*)(&p->torque))[ax]/p->inertia.x;
			else ((Real*)(&angAccel))[ax]=0.;
		}
	}
#endif
	#ifdef DUMP_INTEGRATOR
		if(dofs!=0) printf("$%ld/%d: v=%v3g, ω=%v3g, F=%v3g, T=%v3g, a=%v3g",scene->step,get_global_id(0),p->vel,p->angVel,p->force,p->torque,accel);
	#endif
	if(scene->damping!=0){
		accel   =accel   *(1-scene->damping*sign(p->force *(p->vel   +accel   *scene->dt/2)));
		angAccel=angAccel*(1-scene->damping*sign(p->torque*(p->angVel+angAccel*scene->dt/2)));
	}
	p->vel+=accel*scene->dt;
	p->angVel+=angAccel*scene->dt;
	p->pos+=p->vel*scene->dt;
	#ifdef DUMP_INTEGRATOR
		if(dofs!=0) printf(", aDamp=%v3g, vNew=%v3g, posNew=%v3g\n",accel,p->vel,p->pos);
	#endif
	p->ori=Quat_multQ(Quat_fromRotVec(p->angVel*scene->dt),p->ori);
	//
	// reset forces on particles
	p->force=p->torque=(Vec3)0.;
}

void computeL6GeomGeneric(global struct Contact* c, const Vec3 pos1, const Vec3 vel1, const Vec3 angVel1, const Vec3 pos2, const Vec3 vel2, const Vec3 angVel2, const Vec3 normal, const Vec3 contPt, const Real uN, const Real r1, const Real r2, Real dt){
	// new contact
	if(con_geomT_get(c)==0){
		con_geomT_set(c,Geom_L6Geom); c->geom.l6g=L6Geom_new();
		c->ori=Mat3_rot_setYZ(normal);
		c->pos=contPt;
		c->geom.l6g.uN=uN;
		return;
	}
	Vec3 currNormal=normal;
	Vec3 prevNormal=Mat3_row(c->ori,0);
	Vec3 prevContPt=c->pos;
	Vec3 normRotVec=cross(prevNormal,currNormal);
	Vec3 midNormal=normalize(.5*(prevNormal+currNormal));
	Vec3 normTwistVec=midNormal*dt*.5*dot(midNormal,angVel1+angVel2);
	Mat3 prevOri=c->ori;
	Vec3 midOri1=Mat3_row(prevOri,1)-.5*cross(Mat3_row(prevOri,1),normRotVec+normTwistVec);
	Mat3 midOri=Mat3_setRows(midNormal,midOri1,cross(midNormal,midOri1));
	Vec3 currOri1=Mat3_row(prevOri,1)-cross(midOri1,normRotVec+normTwistVec);
	Mat3 currOri=Mat3_setRows(currNormal,currOri1,cross(currNormal,currOri1));
	currOri=Mat3_orthonorm_r0(currOri);
	Vec3 midContPt=.5*(prevContPt+contPt), midPos1=pos1-vel1*dt/2., midPos2=pos2-vel2*dt/2.;
	Vec3 c1x=midContPt-midPos1, c2x=midContPt-midPos2;
	Vec3 relVel=(vel2+cross(angVel2,c2x))-(vel1+cross(angVel1,c1x));
	// update contact
	c->pos=contPt;
	c->ori=currOri;
	c->geom.l6g.vel=Mat3_multV(midOri,relVel);
	c->geom.l6g.angVel=Mat3_multV(midOri,angVel2-angVel1);
	c->geom.l6g.uN=uN;
}

// #define GEOM_L1GEOM

kernel void contCompute(global const struct Scene* scene, global const struct Particle* par, global struct Contact* con){
	global struct Contact *c=&(con[get_global_id(0)]);
	global const struct Particle* p1=&(par[c->ids.s0]);
	global const struct Particle* p2=&(par[c->ids.s1]);
	/* create geometry for new contacts */
	#ifdef GEOM_L1GEOM
		if(con_geomT_get(c)==0){
			con_geomT_set(c,Geom_L1Geom); c->geom.l1g=L1Geom_new();
		}
		switch (SHAPET2_COMBINE(par_shapeT_get(p1),par_shapeT_get(p2))){
			case SHAPET2_COMBINE(Shape_Sphere,Shape_Sphere): {
				//assert(con_geomT_get(c)==Geom_L1Geom);
				global struct L1Geom* l1g=&(c->geom.l1g);
				Real r1=p1->shape.sphere.radius, r2=p2->shape.sphere.radius;
				Real dist=distance(p1->pos,p2->pos);
				Vec3 normal=(p2->pos-p1->pos)/dist;
				l1g->uN=dist-(r1+r2);
				c->pos=p1->pos+(r1+.5*l1g->uN)*normal;
				c->ori=Mat3_rot_setYZ(normal);
			}
			default: /* error */ ;
		}
	#else
		switch(SHAPET2_COMBINE(par_shapeT_get(p1),par_shapeT_get(p2))){
			case SHAPET2_COMBINE(Shape_Sphere,Shape_Sphere):{
				Real r1=p1->shape.sphere.radius, r2=p2->shape.sphere.radius;
				Real dist=distance(p1->pos,p2->pos);
				Vec3 normal=(p2->pos-p1->pos)/dist;
				Real uN=dist-(r1+r2);
				Vec3 contPt=p1->pos+(r1+.5*uN)*normal;
				computeL6GeomGeneric(c,p1->pos,p1->vel,p1->angVel,p2->pos,p2->vel,p2->angVel,normal,contPt,uN,r1,r2,scene->dt);
			}
			default: /* signal error */ ;
		}
	#endif
	/* update physical params: only if there are no physical params yet */
	if(con_physT_get(c)==0){
		int matId1=par_matId_get(p1), matId2=par_matId_get(p2);
		global const struct Material* m1=&(scene->materials[matId1]);
		global const struct Material* m2=&(scene->materials[matId2]);
		int matT1=mat_matT_get(m1), matT2=mat_matT_get(m2);
		switch(MATT2_COMBINE(matT1,matT2)){
			case MATT2_COMBINE(Mat_ElastMat,Mat_ElastMat):{
				c->phys.normPhys=NormPhys_new();
				con_physT_set(c,Phys_NormPhys);
				/* fixme: d1 should depend on current distance, only indirectly on radius */
				Real r1=(par_shapeT_get(p1)==Shape_Sphere?p1->shape.sphere.radius:INFINITY);
				Real r2=(par_shapeT_get(p2)==Shape_Sphere?p2->shape.sphere.radius:INFINITY);
				Real d1=(par_shapeT_get(p1)==Shape_Sphere?p1->shape.sphere.radius:0);
				Real d2=(par_shapeT_get(p1)==Shape_Sphere?p2->shape.sphere.radius:0);
				Real A=M_PI*min(r1,r2)*min(r1,r2);
				c->phys.normPhys.kN=1./(d1/(A*m1->mat.elast.young)+d2/(A*m2->mat.elast.young));
				// printf("d1=%f, d2=%f, A=%f, E1=%f, E2=%f, kN=%f\n",d1,d2,A,m1->mat.elast.young,m2->mat.elast.young,c->phys.normPhys.kN);
			}
			default: /* error */ ;
		}
	}
	/* contact law */
	int geomT=con_geomT_get(c), physT=con_physT_get(c);
	switch(GEOMT_PHYST_COMBINE(geomT,physT)){
		case GEOMT_PHYST_COMBINE(Geom_L1Geom,Phys_NormPhys):
			c->force=(Vec3)(c->geom.l1g.uN*c->phys.normPhys.kN,0,0);
			c->torque=(Vec3)0.;
		case GEOMT_PHYST_COMBINE(Geom_L6Geom,Phys_NormPhys):{
			// [ HACK: this will be removed once the params are in the material
			const Real ktDivKn=.2;
			// const Real l=(par_shapeT_get(p2)==Shape_Sphere?p2->shape.sphere.radius:1.);
			Real charLen=INFINITY;
			// ]
			Vec3 kntt=(Vec3)(c->phys.normPhys.kN,c->phys.normPhys.kN*ktDivKn,c->phys.normPhys.kN*ktDivKn);
			Vec3 ktbb=kntt/charLen;
			c->force+=scene->dt*c->geom.l6g.vel*kntt;
			c->torque+=scene->dt*c->geom.l6g.angVel*ktbb;
			// set this directly
			((global Real*)&(c->force))[0]=c->phys.normPhys.kN*c->geom.l6g.uN;
		}
		default: /* error */ ;
	}
}

#endif /*__OPENCL_VERSION__*/

#ifdef __OPENCL_VERSION__
	#undef constant
	#undef global
#endif

#endif /* _SCENE_CL */

// this padding is important since otherwise the file is read in a werid way givin errors at the end
