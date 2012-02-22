#pragma once

#include"Particle.cl.h"
#include"Contact.cl.h"
#include"Scene.cl.h"

namespace clDem{

	struct Simulation{
		Scene scene;
		vector<Particle> par;
		vector<Contact> con;
		vector<par_id_t> clumps;
		vector<cl_float> bboxes;

		cl::Platform platform;
		cl::Device device;
		cl::Context context;
		cl::CommandQueue queue;
		cl::Program program;

		long lastWriteStep;
		int lastWriteSubstep;

		struct Buffers{
			size_t parSize,conSize,clumpsSize,bboxesSize;
			cl::Buffer scene,par,con,clumps,bboxes;
		} bufs;

		/* read and write fixed-size objects (all async, use queue.finish() to wait) */
		template<typename T> cl::Buffer writeBuf(const T& obj){
			cl::Buffer buf(context,CL_MEM_READ_WRITE,sizeof(T),NULL);
			queue.enqueueWriteBuffer(buf,CL_FALSE,0,sizeof(T),&obj);
			return buf;
		}
		template<typename T> void readBuf(cl::Buffer& buf, T& obj){
			queue.enqueueReadBuffer(buf,CL_FALSE,0,sizeof(T),&obj);
		}
		/* read and write std::vector containers (all async) */
		template<typename T> cl::Buffer writeVecBuf(const std::vector<T>& obj, size_t& writtenSize){
			if(obj.empty()) throw std::runtime_error("Buffer created from std::vector<T> may not be empty.");
			cl::Buffer buf(context,CL_MEM_READ_WRITE,sizeof(T)*obj.size(),NULL);
			queue.enqueueWriteBuffer(buf,CL_FALSE,0,sizeof(T)*obj.size(),&(obj[0]));
			writtenSize=obj.size();
			return buf;
		}
		template<typename T> void readVecBuf(cl::Buffer& buf, std::vector<T>& obj, const size_t& writtenSize){
			obj.resize(writtenSize); // in case size was changed meanwhile
			queue.enqueueReadBuffer(buf,CL_FALSE,0,sizeof(T)*obj.size(),&(obj[0]));
		}

		void writeBufs(bool wait=true);
		void readBufs(bool wait=true);
		void runKernels(int nSteps, int substepStart=0);
		cl::Kernel makeKernel(const char* name);

		Simulation(int pNum=-1,int dNum=-1, const string& opts=""){
			scene=Scene_new();
			initCl(pNum,dNum,opts);
		}
		void initCl(int pNum, int dNum, const string& opts);
		void run(int nSteps);
		Real pWaveDt();
		py::list saveVtk(string prefix, bool compress=true, bool ascii=false);
	};

};

void Simulation_hpp_expose(){
	py::class_<Simulation>("Simulation",py::init<int,int,string>((py::arg("platformNum")=-1,py::arg("deviceNum")=-1,py::arg("opts")="")))
		.PY_RW(Simulation,scene)
		.PY_RW(Simulation,par)
		.PY_RW(Simulation,con)
		.def("run",&Simulation::run)
		#ifdef CLDEM_VTK
		.def("saveVtk",&Simulation::saveVtk,(py::arg("prefix"),py::arg("compress")=true,py::arg("ascii")=false))
		#endif
		.def("pWaveDt",&Simulation::pWaveDt)
	;
}
