#include<iostream>
#include<fstream>
#include<sstream>
#include<cassert>

#include<CL/cl.h>
#include"scene.cl"


int main(int argc, char* argv[]){
   // initialize OpenCL
   cl_int err;
   cl_platform_id platform;
   cl_device_id device;
   cl_context context;
   cl_command_queue queue;
   err=clGetPlatformIDs(1,&platform,NULL); assert(!err);
   err=clGetDeviceIDs(platform,CL_DEVICE_TYPE_GPU,1,&device,NULL); assert(!err);
	context=clCreateContext(NULL,1,&device,NULL,NULL,&err); assert(!err);
   queue=clCreateCommandQueue(context,device,0,&err); assert(!err);
	// compile source
	std::ifstream sceneCl("scene.cl",std::ios::in); 
	std::ostringstream oss; oss<<"#line 1 \"scene.cl\"\n"<<sceneCl.rdbuf();
	const char* srcStr(oss.str().c_str());
	cl_program prog=clCreateProgramWithSource(context,1,(const char**)&srcStr,NULL,NULL);
	err=clBuildProgram(prog,0,NULL,"-Werror ",NULL,NULL);
	if(err!=CL_SUCCESS){
		char buildLog[1<<16];
		clGetProgramBuildInfo(prog,device,CL_PROGRAM_BUILD_LOG,sizeof(buildLog),buildLog,NULL);
		std::cerr<<buildLog;
		return -1;
	}

	// scene setup
	struct Scene scene=Scene_new();
	scene.dt=1e-8;
	std::cout<<"Starting at step "<<scene.step<<" (t="<<scene.t<<"), Δt="<<scene.dt<<std::endl;

	// create buffers, enqueue copies to the device
	cl_mem sceneBuf=clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(scene),NULL,&err); assert(!err);
	clEnqueueWriteBuffer(queue,sceneBuf,/*blocking*/CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL); 
	// create kernels, set their arguments
	cl_kernel stepK=clCreateKernel(prog,"nextTimestep",&err);
	clSetKernelArg(stepK,0,sizeof(cl_mem),&sceneBuf);
	// enque kernels
	for(int step=0; step<20; step++){
		clEnqueueTask(queue,stepK,0,NULL,NULL);
	}
	clEnqueueReadBuffer(queue,sceneBuf,/*block until done*/CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL);
	std::cout<<"Finished at step "<<scene.step<<" (t="<<scene.t<<"), Δt="<<scene.dt<<std::endl;
}
