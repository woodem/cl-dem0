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
   err=clGetDeviceIDs(platform,CL_DEVICE_TYPE_CPU,1,&device,NULL); assert(!err);
	context=clCreateContext(NULL,1,&device,NULL,NULL,&err); assert(!err);
   queue=clCreateCommandQueue(context,device,0,&err); assert(!err);
	// compile source
	const char* src="#include\"scene.cl\"\n";
	cl_program prog=clCreateProgramWithSource(context,1,(const char**)&src,NULL,&err); assert(!err);
	err=clBuildProgram(prog,0,NULL,"-Werror -I.",NULL,NULL);
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
		}
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
	for(int step=0; step<20000; step++){
		clEnqueueTask(queue,stepK,0,NULL,NULL);
	}
	clEnqueueReadBuffer(queue,sceneBuf,/*block until done*/CL_TRUE,0,sizeof(scene),&scene,0,NULL,NULL);
	std::cout<<"Finished at step "<<scene.step<<" (t="<<scene.t<<"), Δt="<<scene.dt<<std::endl;
}
