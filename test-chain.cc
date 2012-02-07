
#include<vector>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cassert>

#include<CL/cl.h>
#include"scene.cl"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;


const Real r=0.005; // radius
const Real E=1e9;   // young's modulus
const Real rho=1e4; // density


int main(int argc, char* argv[]){
	if(argc<=1){
		cerr<<"Usage: "<<argv[0]<<" platformNumber"<<endl;
		return -1;
	}
	int platformNum=atoi(argv[1]);

   // initialize OpenCL
   cl_int err;
   cl_platform_id platforms[8];
   cl_device_id device;
   cl_context context;
   cl_command_queue queue;
	cl_uint nPlatforms;
   err=clGetPlatformIDs(8,platforms,&nPlatforms); assert(!err);
	//cerr<<"Number of platforms available: "<<nPlatforms<<endl;
	if(nPlatforms==0){ cerr<<"No OpenCL platform available"<<endl; return -1; }
	if(platformNum>=nPlatforms){ cerr<<"Only "<<nPlatforms<<" available, using 0th platform."<<endl; platformNum=0; }
	cl_platform_id platform=platforms[platformNum];
   err=clGetDeviceIDs(platform,CL_DEVICE_TYPE_GPU,1,&device,NULL);
	if(err){ cerr<<"No GPU, trying CPU"<<endl; err=clGetDeviceIDs(platform,CL_DEVICE_TYPE_CPU,1,&device,NULL); }
	assert(!err);
	char pName[128],dName[128];
	clGetPlatformInfo(platform,CL_PLATFORM_NAME,sizeof(pName),pName,NULL);
	clGetDeviceInfo(device,CL_DEVICE_NAME,sizeof(dName),dName,NULL);
	cerr<<"** OpenCL ready: platform \""<<pName<<"\", device \""<<dName<<"\"."<<endl;
	context=clCreateContext(NULL,1,&device,NULL,NULL,&err); assert(!err);
   queue=clCreateCommandQueue(context,device,0,&err); assert(!err);

	size_t N=10;
	size_t Ncon=N-1;
	vector<size_t> supports={0,N-1};

	struct Scene scene=Scene_new();
	// scene setup
	scene.dt=.2*(r/sqrt(E/rho)); // .5 × p-Wave critical timestep
	scene.gravity={0,0,-10};
	scene.damping=.7;

	struct ElastMat em=ElastMat_new();
	em.young=E;
	mat_matT_set(&scene.materials[0],Mat_ElastMat);
	scene.materials[0].mat.elast=em;

	// create chain of particles connected with each other
	vector<Particle> par(N);
	vector<Contact> con(Ncon);
	for(int i=0; i<N; i++){
		par[i]=Particle_new();
		Particle& p(par[i]);
		p.pos={2*r*i,0,0};
		p.mass=rho*(4/3.)*M_PI*pow(r,3);
		Real inert=p.mass*(2/5.)*pow(r,2);
		p.inertia={inert,inert,inert};
		bool isSupport=(std::find(supports.begin(),supports.end(),i)!=supports.end());
		par_dofs_set(&p,isSupport?0:par_dofs_all);
		p.shape.sphere=Sphere_new();
		p.shape.sphere.radius=r;
		par_shapeT_set(&p,Shape_Sphere);
		cerr<<"#"<<i<<", flags="<<p.flags<<endl;
		if(i>0){ // first particle has no contact with the previous one
			con[i-1]=Contact_new();
			con[i-1].ids.s0=i-1;
			con[i-1].ids.s1=i;
		}
	}
	// compile source
	const char* src="#include\"scene.cl\"";
	cl_program prog=clCreateProgramWithSource(context,1,&src,NULL,NULL);
	err=clBuildProgram(prog,0,NULL," -I.",NULL,NULL);
	if(err!=CL_SUCCESS){
		switch(err){
			#define CASE_ERR(e) case e: std::cerr<<#e<<std::endl; break;
			CASE_ERR(CL_INVALID_PROGRAM);
			CASE_ERR(CL_INVALID_VALUE);
			CASE_ERR(CL_INVALID_DEVICE);
			CASE_ERR(CL_INVALID_BUILD_OPTIONS);
			CASE_ERR(CL_COMPILER_NOT_AVAILABLE);
			CASE_ERR(CL_BUILD_PROGRAM_FAILURE);
			CASE_ERR(CL_INVALID_OPERATION);
			CASE_ERR(CL_OUT_OF_HOST_MEMORY);
			default: std::cerr<<"Unknown error from clBuildProgram: "<<err<<std::endl;
		}
		char buildLog[1<<16];
		clGetProgramBuildInfo(prog,device,CL_PROGRAM_BUILD_LOG,sizeof(buildLog),buildLog,NULL);
		std::cerr<<buildLog;
		return -1;
	}

	// create buffers, enqueue copies to the device
	cl_mem sceneBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(scene),NULL,&err);
	clEnqueueWriteBuffer(queue,sceneBuf,CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL);
	cl_mem parBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(Particle)*par.size(),NULL,&err);
	clEnqueueWriteBuffer(queue,parBuf,CL_TRUE,0,sizeof(Particle)*par.size(),&(par[0]),0,NULL,NULL);
	cl_mem conBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(Contact)*con.size(),NULL,&err);
	clEnqueueWriteBuffer(queue,conBuf,CL_TRUE,0,sizeof(Contact)*con.size(),&(con[0]),0,NULL,NULL);

	// create kernels, set their arguments
	cl_kernel stepK=clCreateKernel(prog,"nextTimestep",&err);
	clSetKernelArg(stepK,0,sizeof(cl_mem),&sceneBuf);

	cl_kernel integratorK=clCreateKernel(prog,"integrator",&err);
	clSetKernelArg(integratorK,0,sizeof(cl_mem),&sceneBuf);
	clSetKernelArg(integratorK,1,sizeof(cl_mem),&parBuf);

	cl_kernel contactsK=clCreateKernel(prog,"contCompute",&err);
	clSetKernelArg(contactsK,0,sizeof(cl_mem),&sceneBuf);
	clSetKernelArg(contactsK,1,sizeof(cl_mem),&parBuf);
	clSetKernelArg(contactsK,2,sizeof(cl_mem),&conBuf);

	cl_kernel forcesK=clCreateKernel(prog,"forcesToParticles",&err);
	clSetKernelArg(forcesK,0,sizeof(cl_mem),&sceneBuf);
	clSetKernelArg(forcesK,1,sizeof(cl_mem),&parBuf);
	clSetKernelArg(forcesK,2,sizeof(cl_mem),&conBuf);
	
	for(int bigStep=0; bigStep<1; bigStep++){
		size_t one={1,};
		// put kernels in the queue
		for(int step=0; step<10000; step++){
			clEnqueueTask(queue,stepK,0,NULL,NULL);
			clEnqueueNDRangeKernel(queue,integratorK,1,NULL,&N,   &one,0,NULL,NULL);
			clEnqueueNDRangeKernel(queue,contactsK,  1,NULL,&Ncon,&one,0,NULL,NULL);
			clEnqueueNDRangeKernel(queue,forcesK,    1,NULL,&Ncon,&one,0,NULL,NULL);
		}
		/* blocking reads wait for kernels to finish */
		clEnqueueReadBuffer(queue,sceneBuf,CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL);
		clEnqueueReadBuffer(queue,parBuf,CL_TRUE,0,sizeof(Particle)*par.size(),&(par[0]),0,NULL,NULL);
		clEnqueueReadBuffer(queue,conBuf,CL_TRUE,0,sizeof(Contact)*con.size(),&(con[0]),0,NULL,NULL);
		// write out particle's positions
		cout<<"At step "<<scene.step<<" (t="<<scene.t<<"), Δt="<<scene.dt<<endl;
		for(size_t i=0; i<par.size(); i++){
			cout<<"#"<<i<<" x="<<par[i].pos<<"; v="<<par[i].vel<<endl;
		}
		for(const Contact& c: con){
			cout<<"#"<<c.ids.x<<"+"<<c.ids.y<<": pt="<<c.pos<<endl;
			cout<<"\tgeomT="<<con_geomT_get(&c)<<", uN="<<c.geom.l1g.uN<<", locX="<<c.ori<<endl;
			cout<<"\tphysT="<<con_physT_get(&c)<<", kN="<<c.phys.normPhys.kN<<", F="<<c.force<<endl;
			//cout<<"\tcomb="<<GEOMT_PHYST_COMBINE(con_geomT_get(&c),con_physT_get(&c))<<", handled "<<GEOMT_PHYST_COMBINE(Geom_L1Geom,Phys_NormPhys)<<endl;
		}
	}
	return 0;
}
