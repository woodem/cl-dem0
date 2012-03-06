/* vim:set syntax=c: */
#ifndef _KERNELS_CL_H
#define _KERNELS_CL_H

#include"Particle.cl.h"
#include"Contact.cl.h"
#include"Collider.cl.h"
#include"Scene.cl.h"

CLDEM_NAMESPACE_BEGIN()

// substep numbers here
enum _substeps{ SUB_nextTimestep=0, SUB_integrator, SUB_updateBboxes, SUB_checkPotCon, SUB_contCompute, SUB_forcesToParticles };

#ifdef __cplusplus
/* how is kernel parallellized: single (task), over particles, over contacts, over potential contacts */
enum _kargs{ KARGS_SINGLE, KARGS_PAR, KARGS_CON, KARGS_POT };
/* information about kernel necessary for the queueing and interrupt logic */
struct KernelInfo {
	int substep;
	const char* name;
	int argsType;
};
/* this is the loop which will be running */
static struct KernelInfo clDemKernels[]={
	{ SUB_nextTimestep     , "nextTimestep_1"    , KARGS_SINGLE },
	{ SUB_integrator       , "integrator_P"      , KARGS_PAR },
	{ SUB_updateBboxes     , "updateBboxes_P"    , KARGS_PAR },
	{ SUB_checkPotCon      , "checkPotCon_PC"    , KARGS_POT },
	{ SUB_contCompute      , "contCompute_C"     , KARGS_CON },
	{ SUB_forcesToParticles, "forcesToParticles_C", KARGS_CON },
	{ }, /*sentinel*/
};
#endif

CLDEM_NAMESPACE_END()

/* the rest of the code is only used on the device */
#ifdef __OPENCL_VERSION__



// all kernels take the same set of arguments, for simplicity in the host code
#define KERNEL_ARGUMENT_LIST global struct Scene* scene, global struct Particle* par, global struct Contact* con, global int *conFree, global long2* pot, global int *potFree, global struct CJournalItem* cJournal, global int* clumps, global float* bboxes

#ifdef cl_amd_printf
	#define PRINT_TRACE(name)
#else
	#define PRINT_TRACE(name) // { if(get_global_id(0)==0) printf("%s:%-4d %4ld/%d %s\n",__FILE__,__LINE__,scene->step,substep,name); }
#endif

#if 1
	#define TRYLOCK(a) atom_cmpxchg(a,0,1)
	#define LOCK(a) while(TRYLOCK(a))
	#define UNLOCK(a) atom_xchg(a,0)
#else
	#define TRYLOCK(a)
	#define LOCK(a)
	#define UNLOCK(a)
#endif

#ifdef TRACK_ENERGY
	#define ADD_ENERGY(scene,name,E) Scene_energyAdd(scene,ENERGY_##name,E);
#else
	#define ADD_ENERGY(scene,name,E)
#endif

// if defined, printf's show when contacts are created/deleted
#define CON_LOG

void Scene_energyAdd(global struct Scene* scene, int energyIndex, Real E){
	//atomic_add(&(scene->energy[energyIndex]),(float)E);
	global int* mutex=&(scene->energyMutex[energyIndex]);
	LOCK(mutex);
	scene->energy[energyIndex]+=(float)E;
	UNLOCK(mutex);
}
void Scene_energyZeroNonincremental(global struct Scene* scene){
	// reset non-incremental energies
	//printf("[");
	for(int i=0; i<SCENE_ENERGY_NUM; i++){
		//printf("%d ",(int)energyDefinitions[i].incremental);
		if(!energyDefinitions[i].incremental) scene->energy[i]=0.;
	}
	//printf("]\n");
}

/**** kernel code, will be in separate files ****/


kernel void nextTimestep_1(KERNEL_ARGUMENT_LIST){
	const int substep=SUB_nextTimestep;
	scene->step++;
	//PRINT_TRACE("** nextTimestep_1");
	if(scene->step>0 && Scene_skipKernel(/*use new value already, but without assigning it to scene*/scene,substep)) return;
	PRINT_TRACE("nextTimestep_1");

	// if step is -1, we are at the very beginning
	// keep time at 0, only increase step number
	if(scene->step>0) scene->t+=scene->dt;

	scene->updateBboxes=false;

	#ifdef TRACK_ENERGY
		Scene_energyZeroNonincremental(scene);
	#endif
}

kernel void integrator_P(KERNEL_ARGUMENT_LIST){
	const int substep=SUB_integrator;
	if(Scene_skipKernel(scene,substep)) return;
	PRINT_TRACE("integrator_P");

	global struct Particle *p=&(par[get_global_id(0)]);
	struct Particle pp=*p;
	Real dt=scene->dt;
	Vec3 gravity=scene->gravity;
	Real damping=scene->damping;
	Real verletDist=scene->verletDist;

	int dofs=par_dofs_get_local(&pp);
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
	pp.force+=gravity*pp.mass; // put this up here so that grav appears in the force term

	// some translations are allowed, compute gravity contribution
	if((dofs & par_dofs_trans)!=0) ADD_ENERGY(scene,grav,-dot(gravity,pp.vel)*pp.mass*dt);

	if((dofs&par_dofs_trans)==par_dofs_trans){ // all translations possible, use vector expr
		accel+=pp.force/pp.mass;
	} else {
		for(int ax=0; ax<3; ax++){
			if(dofs & dof_axis(ax,false)) ((Real*)(&accel))[ax]+=((Real*)(&(pp.force)))[ax]/pp.mass+((Real*)(&(gravity)))[ax];
			else ((Real*)(&accel))[ax]=0.;
		}
	}
	if((dofs&par_dofs_rot)==par_dofs_rot){
		angAccel+=pp.torque/pp.inertia.x;
	} else {
		for(int ax=0; ax<3; ax++){
			if(dofs & dof_axis(ax,true )) ((Real*)(&angAccel))[ax]+=((Real*)(&pp.torque))[ax]/pp.inertia.x;
			else ((Real*)(&angAccel))[ax]=0.;
		}
	}
#endif
	#ifdef DUMP_INTEGRATOR
		if(dofs!=0) printf("$%ld/%d: v=%v3g, ω=%v3g, F=%v3g, T=%v3g, a=%v3g",scene->step,get_global_id(0),pp.vel,pp.angVel,pp.force,pp.torque,accel);
	#endif
	if(damping!=0){
		ADD_ENERGY(scene,damp,(dot(fabs(pp.vel),fabs(pp.force))+dot(fabs(pp.angVel),fabs(pp.torque)))*damping*dt);
		accel   =accel   *(1-damping*sign(pp.force *(pp.vel   +accel   *dt/2)));
		angAccel=angAccel*(1-damping*sign(pp.torque*(pp.angVel+angAccel*dt/2)));
	}

	// this estimates velocity at next on-step using current accel; it follows how yade computes it currently
	ADD_ENERGY(scene,Ekt,.5*pp.mass*Vec3_sqNorm(pp.vel+.5*dt*accel));
	// valid only for spherical particles (fix later)
	ADD_ENERGY(scene,Ekr,.5*dot(pp.inertia*(pp.angVel+.5*dt*angAccel),pp.angVel+.5*dt*angAccel));

	pp.vel+=accel*dt;
	pp.angVel+=angAccel*dt;
	pp.pos+=pp.vel*dt;
	#ifdef DUMP_INTEGRATOR
		if(dofs!=0) printf(", aDamp=%v3g, vNew=%v3g, posNew=%v3g\n",accel,pp.vel,pp.pos);
	#endif
	pp.ori=Quat_multQ(Quat_fromRotVec(pp.angVel*dt),pp.ori); // checks automatically whether |rotVec|==0
	//
	// reset forces on particles
	// pp.force=pp.torque=(Vec3)0.;

	// is this safe for multi-threaded? is the new value always going to be written?
	// it is not read until in the next kernel
	#if 0
		printf("#%ld boxPos=%g\n",get_global_id(0),p->bboxPos.s0);
		printf("#%ld boxPos=%g\n",get_global_id(0),p->bboxPos.s1);
		printf("#%ld boxPos=%g\n",get_global_id(0),p->bboxPos.s2);
	#endif
	if(verletDist>=0 && (Vec3_sqNorm(pp.pos-pp.bboxPos)>pown(verletDist,2) || isnan(pp.bboxPos.s0))) scene->updateBboxes=true;

	/* write back */
	p->pos=pp.pos;
	p->ori=pp.ori;
	p->vel=pp.vel;
	p->angVel=pp.angVel;
	p->force=(Vec3)0.;
	p->torque=(Vec3)0.;
}

/** This kernel runs only when INT_OUT_OF_BBOX is set **/
kernel void updateBboxes_P(KERNEL_ARGUMENT_LIST){
	const int substep=SUB_updateBboxes;
	if(!scene->updateBboxes) return;
	if(Scene_skipKernel(scene,substep)) return;
	PRINT_TRACE("updateBboxes_P");

	par_id_t id=get_global_id(0);
	global struct Particle *p=&par[id];
	struct Particle pp=*p;
	Vec3 mn, mx;
	switch(par_shapeT_get_local(&pp)){
		case Shape_Sphere:
			mn=pp.pos-(Vec3)(pp.shape.sphere.radius);
			mx=pp.pos+(Vec3)(pp.shape.sphere.radius);
			break;
		case Shape_Wall:
			mn=(Vec3)-INFINITY; mx=(Vec3)INFINITY;
			if(pp.shape.wall.axis==0) mn.s0=mx.s0=p->pos.s0;
			else if(pp.shape.wall.axis==1) mn.s1=mx.s1=p->pos.s1;
			else if(pp.shape.wall.axis==2) mn.s2=mx.s2=p->pos.s2;
			else /* error */ printf("ERROR: #%ld is a wall with axis=%d!\n",id,pp.shape.wall.axis);
		default: /* */;
	}
	p->bboxPos=pp.pos;
	mn-=(Vec3)scene->verletDist;
	mx+=(Vec3)scene->verletDist;
	bboxes[id*6+0]=mn.x; bboxes[id*6+1]=mn.y; bboxes[id*6+2]=mn.z;
	bboxes[id*6+3]=mx.x; bboxes[id*6+4]=mx.y; bboxes[id*6+5]=mx.z;
	//printf("%ld: %v3g±(%g+%g) %v3g %v3g\n",id,p->pos,scene->verletDist,p->shape.sphere.radius,mn,mx);

	// not immediate since all threads must finish
	if(id==0) Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE | INT_NOT_DESTRUCTIVE | INT_BBOXES);
}

kernel void checkPotCon_PC(KERNEL_ARGUMENT_LIST){
	const int substep=SUB_checkPotCon;
	if(Scene_skipKernel(scene,substep)) return;
	size_t cid=get_global_id(0);
	if(cid>=scene->arrSize[ARR_POT]) return; // in case we are past real number of contacts
	par_id2_t ids=pot[cid];
	if(ids.s0<0 || ids.s1<0) return; // deleted contact
	PRINT_TRACE("checkPotCon_PC");
	
	// always make contacts such that shape index of the first particle <= shape index of the second particle
	int flip=par_shapeT_get(&par[ids.s0])>par_shapeT_get(&par[ids.s1]);
	if(flip) ids=(par_id2_t)(ids.s1,ids.s0);

	const struct Particle p1=par[ids.s0];
	const struct Particle p2=par[ids.s1];
	switch(SHAPET2_COMBINE(par_shapeT_get_local(&p1),par_shapeT_get_local(&p2))){
		case SHAPET2_COMBINE(Shape_Sphere,Shape_Sphere):{
			Real r1=p1.shape.sphere.radius, r2=p2.shape.sphere.radius;
			if(distance(p1.pos,p2.pos)>=r1+r2){
				// printf("pot ##%d+%d: distance %g.\n",(int)ids.s0,(int)ids.s1,distance(p1->pos,p2->pos)-(r1+r2));
				return;
			}
			break;
		}
		case SHAPET2_COMBINE(Shape_Sphere,Shape_Wall):{
			if(fabs(((Real*)(&(p1.pos)))[p2.shape.wall.axis]-((Real*)(&(p2.pos)))[p2.shape.wall.axis])>=p1.shape.sphere.radius) return;
			break;
		}
		default: printf("ERROR: Invalid shape indices in pot ##%ld+%ld: %d+%d (%d; sphere+sphere=%d, sphere+wall=%d)!\n",ids.s0,ids.s1,par_shapeT_get_local(&p1),par_shapeT_get_local(&p2),SHAPET2_COMBINE(par_shapeT_get_local(&p1),par_shapeT_get_local(&p2)),SHAPET2_COMBINE(Shape_Sphere,Shape_Sphere),SHAPET2_COMBINE(Shape_Sphere,Shape_Wall)); return;
	};
	#ifdef CON_LOG
		printf("Creating ##%ld+%ld\n",ids.s0,ids.s1);
	#endif
	// the contact is real, create it
	// this will be moved to a separate function later

	// delete from pot
	pot[cid]=(par_id2_t)(-1,-1);

	// alocate indices; group together so that the interrupt handler reallocates all of them at once
	int ixPotFree=Scene_arr_findNegative_or_append(scene,ARR_POTFREE,potFree); // index where to write free slot
	int ixCon=Scene_arr_fromFreeArr_or_append(scene,conFree,ARR_CONFREE,ARR_CON,/*shrink*/true);
	int ixCJournal=Scene_arr_append(scene,ARR_CJOURNAL);
	// printf("Indices %d %d %d\n",ixPotFree,ixCon,ixCJournal);

	// add to potFree
	// not immediate, since other work-items might increase the required capacity yet more
	if(ixPotFree<0){
		//printf("!! potfree insufficient, size=%d, alloc=%d\n",scene->arrSize[ARR_POTFREE],scene->arrAlloc[ARR_POTFREE]);
		Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE|INT_DESTRUCTIVE|INT_ARRAYS); return; }
	potFree[ixPotFree]=cid;
	
	// add to con
	if(ixCon<0){ Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE|INT_DESTRUCTIVE|INT_ARRAYS); return; }
	// actually create the new contact here
	#if 1
		// debugging only (the race conditions is fixed now, so it should not happen anymore)
		if(con[ixCon].ids.s0>0){ printf("ERROR: con[%d] reported as free, but contains ##%ld+%ld\n",ixCon,con[ixCon].ids.s0,con[ixCon].ids.s1); }
	#endif
	Contact_init(&(con[ixCon]));
	con[ixCon].ids=ids;

	// add the change to cJournal
	if(ixCJournal<0){ Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE|INT_DESTRUCTIVE|INT_ARRAYS); return; }

	struct CJournalItem i; i.ids=ids; i.index=ixCon; i.what=CJOURNAL_POT2CON;
	cJournal[ixCJournal]=i;

	#if 0
		CJournalItem_init(&(cJournal[ixCJournal]));
		cJournal[ixCJournal].ids=ids; cJournal[ixCJournal].index=ixCon; cJournal[ixCJournal].what=CJOURNAL_POT2CON;
	#endif
};

bool Bbox_overlap(global float* _A, global float* _B){
	float8 A=vload8(0,_A), B=vload8(0,_B);
	return
		A.s0<B.s3 && B.s0<A.s3 && // xMinA<xMaxB && xMinB<xMaxA
		A.s1<B.s4 && B.s1<A.s4 &&
		A.s2<B.s5 && B.s2<A.s5;
}

void computeL6GeomGeneric(struct Contact* c, const Vec3 pos1, const Vec3 vel1, const Vec3 angVel1, const Vec3 pos2, const Vec3 vel2, const Vec3 angVel2, const Vec3 normal, const Vec3 contPt, const Real uN, Real dt){
	// new contact
	if(con_geomT_get_local(c)==0){
		con_geomT_set_local(c,Geom_L6Geom); c->geom.l6g=L6Geom_new();
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
	Vec3 prevOri1=Mat3_row(c->ori,1);
	Vec3 midOri1=prevOri1-.5*cross(prevOri1,normRotVec+normTwistVec);
	Mat3 midOri=Mat3_setRows(midNormal,midOri1,cross(midNormal,midOri1));
	//midOri=Mat3_orthonorm_row0(midOri); // not clear how much difference this one makes
	Vec3 currOri1=prevOri1-cross(midOri1,normRotVec+normTwistVec);
	Mat3 currOri=Mat3_setRows(currNormal,currOri1,cross(currNormal,currOri1));
	currOri=Mat3_orthonorm_row0(currOri);
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


kernel void contCompute_C(KERNEL_ARGUMENT_LIST){
	const int substep=SUB_contCompute;
	if(Scene_skipKernel(scene,substep)) return;
	size_t cid=get_global_id(0);
	if(cid>=scene->arrSize[ARR_CON]) return; // in case we are past real number of contacts
	struct Contact c=con[cid];
	if(c.ids.s0<0 || c.ids.s1<0) return; // deleted contact
	PRINT_TRACE("contCompute_C");

	const Real dt=scene->dt;

	bool contactBroken=false;
	const struct Particle p1=par[c.ids.s0];
	const struct Particle p2=par[c.ids.s1];
	if(par_shapeT_get_local(&p1)>par_shapeT_get_local(&p2)) printf("ERROR: ##%ld+%ld is not ordered by shape indices: %d+%d\n",c.ids.s0,c.ids.s1,par_shapeT_get_local(&p1),par_shapeT_get_local(&p2));

	/* create geometry for new contacts */
	switch(SHAPET2_COMBINE(par_shapeT_get_local(&p1),par_shapeT_get_local(&p2))){
		case SHAPET2_COMBINE(Shape_Sphere,Shape_Sphere):{
			Real r1=p1.shape.sphere.radius, r2=p2.shape.sphere.radius;
			Real dist=distance(p1.pos,p2.pos);
			Vec3 normal=(p2.pos-p1.pos)/dist;
			Real uN=dist-(r1+r2);
			#ifdef L6GEOM_BREAK_TENSION
				if(uN>0) contactBroken=true;
			#endif
			Vec3 contPt=p1.pos+(r1+.5*uN)*normal;
			computeL6GeomGeneric(&c,p1.pos,p1.vel,p1.angVel,p2.pos,p2.vel,p2.angVel,normal,contPt,uN,dt);
			break;
		}
		case SHAPET2_COMBINE(Shape_Sphere,Shape_Wall):{
			short axis=p2.shape.wall.axis, sense=p2.shape.wall.sense;
			// coordinates along wall normal axis
			Real cS=((Real*)&p1.pos)[axis], cW=((Real*)&p2.pos)[axis]; 
			Real signedDist=cW-cS;
			Vec3 normal=(Vec3)0.;
			// if sense==0, the normal is oriented from the wall towards the sphere's center
			if(sense==0) ((Real*)&normal)[axis]=signedDist>0?1.:-1.;
			// else it is always oriented either along +axis or -axis
			else ((Real*)&normal)[axis]=(sense>0?-1.:1.);
			Real uN=((Real*)&normal)[axis]*signedDist-p1.shape.sphere.radius;
			//printf("uN=%g, normal=%v3g, signedDist=%g, radius=%g\n",uN,normal,signedDist,p1->shape.sphere.radius);
			#ifdef L6GEOM_BREAK_TENSION
				if(uN>0) contactBroken=true;
			#endif
			// project sphere onto the wall
			Vec3 contPt=p1.pos; ((Real*)(&contPt))[axis]=((Real*)&p2.pos)[axis];
			computeL6GeomGeneric(&c,p1.pos,p1.vel,p1.angVel,p2.pos,p2.vel,p2.angVel,normal,contPt,uN,dt);
			break;
		}
		default: printf("ERROR: ##%ld+%ld has unknown shape index combination %d+%d",c.ids.s0,c.ids.s1,par_shapeT_get_local(&p1),par_shapeT_get_local(&p2));
	}
	/* update physical params: only if there are no physical params yet */
	if(!contactBroken && con_physT_get_local(&c)==0){
		int matId1=par_matId_get_local(&p1), matId2=par_matId_get_local(&p2);
		const struct Material m1=scene->materials[matId1];
		const struct Material m2=scene->materials[matId2];
		int matT1=mat_matT_get_local(&m1), matT2=mat_matT_get_local(&m2);
		switch(MATT2_COMBINE(matT1,matT2)){
			case MATT2_COMBINE(Mat_ElastMat,Mat_ElastMat):{
				c.phys.normPhys=NormPhys_new();
				con_physT_set_local(&c,Phys_NormPhys);
				/* fixme: d1 should depend on current distance, only indirectly on radius */
				Real r1=(par_shapeT_get_local(&p1)==Shape_Sphere?p1.shape.sphere.radius:INFINITY);
				Real r2=(par_shapeT_get_local(&p2)==Shape_Sphere?p2.shape.sphere.radius:INFINITY);
				Real d1=(par_shapeT_get_local(&p1)==Shape_Sphere?p1.shape.sphere.radius:0);
				Real d2=(par_shapeT_get_local(&p1)==Shape_Sphere?p2.shape.sphere.radius:0);
				Real A=M_PI*min(r1,r2)*min(r1,r2);
				c.phys.normPhys.kN=1./(d1/(A*m1.mat.elast.young)+d2/(A*m2.mat.elast.young));
				// printf("d1=%f, d2=%f, A=%f, E1=%f, E2=%f, kN=%f\n",d1,d2,A,m1->mat.elast.young,m2->mat.elast.young,c->phys.normPhys.kN);
				break;
			}
			default: /* error */ ;
		}
	}
	/* contact law */
	if(!contactBroken){
		int geomT=con_geomT_get_local(&c), physT=con_physT_get_local(&c);
		switch(GEOMT_PHYST_COMBINE(geomT,physT)){
			case GEOMT_PHYST_COMBINE(Geom_L1Geom,Phys_NormPhys):
				c.force=(Vec3)(c.geom.l1g.uN*c.phys.normPhys.kN,0,0);
				c.torque=(Vec3)0.;
				break;
			case GEOMT_PHYST_COMBINE(Geom_L6Geom,Phys_NormPhys):{
				#ifndef BEND_CHARLEN
					#define BEND_CHARLEN INFINITY
				#endif
				#ifndef SHEAR_KT_DIV_KN
					#define SHEAR_KT_DIV_KN 0.
				#endif
				// [ HACK: this will be removed once the params are in the material
				Real charLen=BEND_CHARLEN;
				const Real ktDivKn=SHEAR_KT_DIV_KN;
				// ]
				Vec3 kntt=(Vec3)(c.phys.normPhys.kN,c.phys.normPhys.kN*ktDivKn,c.phys.normPhys.kN*ktDivKn);
				Vec3 ktbb=kntt/charLen;
				c.force+=dt*c.geom.l6g.vel*kntt;
				c.torque+=dt*c.geom.l6g.angVel*ktbb;
				((Real*)&(c.force))[0]=c.phys.normPhys.kN*c.geom.l6g.uN; // set this one directly
				#ifdef TRACK_ENERGY
					Real E=.5*pown(c.force.s0,2)/kntt.s0;  // normal stiffness is always non-zero
					if(kntt.s1!=0.) E+=.5*dot(c.force.s12,c.force.s12/kntt.s12);
					if(ktbb.s0!=0.) E+=.5*pown(c.torque.s0,2)/ktbb.s0;
					if(ktbb.s1!=0.) E+=.5*dot(c.torque.s12,c.torque.s12/ktbb.s12);
					ADD_ENERGY(scene,elast,E);
					// gives NaN when some stiffness is 0
					// ADD_ENERGY(scene,elast,.5*dot(c.force,c.force/kntt)+.5*dot(c.torque,c.torque/ktbb));
				#endif
				break;
			}
			default: /* error */ ;
		}
		// wrice contact back to global mem
		con[cid]=c;
	} 

	if(contactBroken){
		#ifdef CON_LOG
			printf("Breaking ##%ld+%ld\n",c.ids.s0,c.ids.s1);
		#endif
		// append contact to conFree (with possible allocation failure)
		// if there is still bbox overlap, append to pot
		// delete contact from con
		int ixConFree=Scene_arr_findNegative_or_append(scene,ARR_CONFREE,conFree); // index where to write free slot
		int ixCJournal=Scene_arr_append(scene,ARR_CJOURNAL);
		//printf("ix=%d\n",ix);
		if(ixConFree<0){ Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE|INT_DESTRUCTIVE|INT_ARRAYS); return; }
		conFree[ixConFree]=cid;
		if(ixCJournal<0){ Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE|INT_DESTRUCTIVE|INT_ARRAYS); return; }
		// if there is overlap, put back to potential contacts
		//printf("bboxes (%v3g)--(%v3g)  (%v3g)--(%v3g)\n",(Vec3)(bboxes[6*c.ids.s0],bboxes[6*c.ids.s0+1],bboxes[6*c.ids.s0+2]),(Vec3)(bboxes[6*c.ids.s0+3],bboxes[6*c.ids.s0+4],bboxes[6*c.ids.s0+5]),(Vec3)(bboxes[6*c.ids.s1],bboxes[6*c.ids.s1+1],bboxes[6*c.ids.s1+2]),(Vec3)(bboxes[6*c.ids.s1+3],bboxes[6*c.ids.s1+4],bboxes[6*c.ids.s1+5]));
		if(Bbox_overlap(&(bboxes[6*c.ids.s0]),&(bboxes[6*c.ids.s1]))){
			int ixPot=Scene_arr_fromFreeArr_or_append(scene,potFree,ARR_POTFREE,ARR_POT,/*shrink*/true);
			if(ixPot<0){ Scene_interrupt_set(scene,substep,INT_NOT_IMMEDIATE|INT_DESTRUCTIVE|INT_ARRAYS); return; }
			pot[ixPot]=c.ids;
			// write to the log
			cJournal[ixCJournal].ids=c.ids; cJournal[ixCJournal].index=ixPot; cJournal[ixCJournal].what=CJOURNAL_CON2POT;
		} else {
			// write to the log
			cJournal[ixCJournal].ids=c.ids; cJournal[ixCJournal].index=-1; cJournal[ixCJournal].what=CJOURNAL_CON_DEL; 
		}
		// delete
		Contact_init(&con[cid]); // or just set ids=(-1,-1)?
	}
}


kernel void forcesToParticles_C(KERNEL_ARGUMENT_LIST){
	const int substep=SUB_forcesToParticles;
	if(Scene_skipKernel(scene,substep)) return;
	PRINT_TRACE("forcesToParticle_C");

	/* how to synchronize access to particles? */
	const struct Contact c=con[get_global_id(0)];
	if(con_geomT_get_local(&c)==0) return;
	global struct Particle *p1=&(par[c.ids.x]), *p2=&(par[c.ids.y]);
	Mat3 R_T=Mat3_transpose(c.ori);
	Vec3 Fp1=+Mat3_multV(R_T,c.force);
	Vec3 Fp2=-Mat3_multV(R_T,c.force);
	Vec3 Tp1=+cross(c.pos-p1->pos,Mat3_multV(R_T,c.force))+Mat3_multV(R_T,c.torque);
	Vec3 Tp2=-cross(c.pos-p2->pos,Mat3_multV(R_T,c.force))-Mat3_multV(R_T,c.torque);
	LOCK(&p1->mutex);
		p1->force+=Fp1; p1->torque+=Tp1;
	UNLOCK(&p1->mutex);
	LOCK(&p2->mutex);
		p2->force+=Fp2; p2->torque+=Tp2;
	UNLOCK(&p2->mutex);
}
#endif
#endif
