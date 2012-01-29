#include<CL/cl.h>

#include"scene.cl"
#include<vector>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
using std::vector;
using std::cerr;
using std::cout;
using std::endl;

const Real r=0.005; // radius
const Real E=1e6;   // young's modulus
const Real rho=1e3; // density


int main(int argc, char* argv[]){
	if(argc<=1){
		cerr<<"Usage: "<<argv[0]<<" N support1 support2 ..."<<endl;
		return 1;
	}
	size_t N=atoi(argv[1]);
	size_t Ncon=N-1;
	vector<int> supports;
	for(int i=2; i<argc; i++) supports.push_back(atoi(argv[i]));

	struct Scene scene=Scene_new();
	// scene setup
	scene.dt=.5*(r/sqrt(E/rho)); // .5 × p-Wave critical timestep
	scene.gravity=Vec3_set(0,0,-10);
	scene.damping=.3;

	struct ElastMat em=ElastMat_new();
	em.young=E;
	mat_matT_set(&scene.materials[0],Mat_ElastMat);
	scene.materials[0].mat.elast=em;

	// create chain of particles connected with each other
	vector<Particle> par(N);
	vector<Contact> con(Ncon);
	for(int i=0; i<N; i++){
		Particle& p(par[i]);
		p=Particle_new();
		p.pos=Vec3_set(2*r*i,0,0);
		p.mass=rho*(4/3.)*M_PI*pow(r,3);
		Real inert=p.mass*(2/5.)*pow(r,2);
		p.inertia=Vec3_set(inert,inert,inert);
		par_dofs_set(&p,(i==0||std::find(supports.begin(),supports.end(),i)!=supports.end()) ? 0 : par_dofs_all);
		if(i>0){ // first particle has no contact with the previous one
			con[i-1]=Contact_new();
			con[i-1].ids.s0=i-1;
			con[i-1].ids.s1=i;
		}
	}

	/* setup OpenCL */
   cl_int err;
   cl_platform_id platform;
   cl_device_id device;
   cl_context context;
   cl_command_queue queue;
   //Initialization
   err=clGetPlatformIDs(1,&platform,NULL);
   err=clGetDeviceIDs(platform,CL_DEVICE_TYPE_GPU,1,&device,NULL);
	context=clCreateContext(NULL,1,&device,NULL,NULL,&err);
   queue=clCreateCommandQueue(context,device,0,&err);

	// create buffers, enqueue copies to the device
	cl_mem sceneBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(scene),NULL,&err);
	clEnqueueWriteBuffer(queue,sceneBuf,CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL);
	cl_mem parBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(Particle)*par.size(),NULL,&err);
	clEnqueueWriteBuffer(queue,parBuf,CL_TRUE,0,sizeof(Particle)*par.size(),&(par[0]),0,NULL,NULL);
	cl_mem conBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(Contact)*con.size(),NULL,&err);
	clEnqueueWriteBuffer(queue,conBuf,CL_TRUE,0,sizeof(Contact)*con.size(),&(con[0]),0,NULL,NULL);
	//
	std::ifstream sceneCl("scene.cl",std::ios::in); 
	std::ostringstream oss; oss<<sceneCl.rdbuf();
	const char* srcStr(oss.str().c_str());
	cl_program prog=clCreateProgramWithSource(context,1,(const char**)&srcStr,NULL,NULL);
	err=clBuildProgram(prog,0,NULL,NULL,NULL,NULL);
	/* enqueue kernels (args: scene, par, con) */
	vector<cl_kernel> kernels;
	for(int step=0; step<100; step++){
		kernels.push_back(clCreateKernel(prog,"nextTimestep",&err));
		clSetKernelArg(*kernels.rbegin(),0,sizeof(cl_mem),&sceneBuf);
		clEnqueueTask(queue,*kernels.rbegin(),0,NULL,NULL);

		kernels.push_back(clCreateKernel(prog,"forcesToParticles",&err));
		clSetKernelArg(*kernels.rbegin(),0,sizeof(cl_mem),&sceneBuf);
		clSetKernelArg(*kernels.rbegin(),1,sizeof(cl_mem),&parBuf);
		clSetKernelArg(*kernels.rbegin(),2,sizeof(cl_mem),&conBuf);
		clEnqueueNDRangeKernel(queue,*kernels.rbegin(),1,NULL,&Ncon,NULL,0,NULL,NULL);

		kernels.push_back(clCreateKernel(prog,"integrator",&err));
		clSetKernelArg(*kernels.rbegin(),0,sizeof(cl_mem),sceneBuf);
		clSetKernelArg(*kernels.rbegin(),1,sizeof(cl_mem),&parBuf);
		clEnqueueNDRangeKernel(queue,*kernels.rbegin(),1,NULL,&N,NULL,0,NULL,NULL);

		kernels.push_back(clCreateKernel(prog,"contCompute",&err));
		clSetKernelArg(*kernels.rbegin(),0,sizeof(cl_mem),&sceneBuf);
		clSetKernelArg(*kernels.rbegin(),1,sizeof(cl_mem),&parBuf);
		clSetKernelArg(*kernels.rbegin(),2,sizeof(cl_mem),&conBuf);
		clEnqueueNDRangeKernel(queue,*kernels.rbegin(),1,NULL,&Ncon,NULL,0,NULL,NULL);
	}
	/* blocking reads wiat for kernels to finish */
	clEnqueueReadBuffer(queue,sceneBuf,CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL);
	clEnqueueReadBuffer(queue,parBuf,CL_TRUE,0,sizeof(Particle)*par.size(),&(par[0]),0,NULL,NULL);
	clEnqueueReadBuffer(queue,conBuf,CL_TRUE,0,sizeof(Contact)*con.size(),&(con[0]),0,NULL,NULL);
	// write out particle's positions
	cout<<"Finished at step "<<scene.step<<" (t="<<scene.t<<"), Δt="<<scene.dt<<endl;
	for(size_t i=0; i<par.size(); i++){
		cout<<"#"<<i<<": "<<par[i].pos.x<<", "<<par[i].pos.y<<", "<<par[i].pos.z<<endl;
	}
	return 0;
}
