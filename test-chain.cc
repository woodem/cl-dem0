#include<CL/cl.h>

#include"scene.cl"
#include<vector>
#include<iostream>
#include<algorithm>
using std::vector;
using std::cerr;
using std::endl;

const Real r=0.005; // radius
const Real E=1e6;   // young's modulus
const Real rho=1e3; // density

cl_command_queue getOpenCLContext(){
    cl_int err;
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue commandQueue;
    //Initialization
    err=clGetPlatformIDs(1,&platform,NULL);
    err=clGetDeviceIDs(platform,CL_DEVICE_TYPE_GPU,1,&device,NULL);
	 return clCreateContext(NULL,1,&device,NULL,NULL,&err);
	 return commandQueue;
}


int main(int argc, char* argv[]){
	if(argc<=1){
		cerr<<"Usage: "<<argv[0]<<" N support1 support2 ..."<<endl;
		return 1;
	}
	int N=atoi(argv[1]);
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
	vector<Contact> con(N-1);
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
	sceneBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(scene),NULL,&err);
	parBug=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(Particle)*par.size(),NULL,&err);
	conBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(Contact)*con.size(),NULL,&err);
	/* enqueue kernels (args: scene, par, con) */
	/* gather results */
}
