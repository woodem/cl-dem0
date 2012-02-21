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

		struct Buffers{
			size_t parSize,conSize;
			cl::Buffer scene,par,con,clumps,bboxes;
		} bufs;

		void writeBufs(bool wait=true);
		void readBufs(bool wait=true);
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
