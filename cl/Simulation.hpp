#pragma once

#include"common.cl.h"
#include"serialization.cl.h"
#include"ObjectIO.hpp"

#include"Particle.cl.h"
#include"Contact.cl.h"
#include"Scene.cl.h"
#include"Collider.cl.h"

// if included from yade
#ifdef YADE_VTK
	#define CLDEM_VTK
#endif

namespace clDem{
	struct Simulation{
		Scene scene;
		cl::Buffer sceneBuf;

		enum { _par=0,_con,_conFree,_pot,_potFree,_cJournal,_clumps,_bboxes,_arrMax };
		struct BufSize{ cl::Buffer buf; size_t size; };
		BufSize bufSize[_arrMax];

		vector<Particle> par;
		vector<Contact> con;
		vector<cl_int> conFree; // free slots in con
		vector<par_id2_t> pot; // potential contacts (only the id1,id2-tuple)
		vector<cl_int> potFree; // free slots in pot
		vector<CJournalItem> cJournal; // logging changes in contact arrays so that the collider's internal structures can be updated accordingly
		vector<ClumpMember> clumps;

		vector<cl_float> bboxes;

		// OopenCL device
		int pNum, dNum;
			
		// numerical parameters
		bool trackEnergy;
		Real ktDivKn;
		bool breakTension;
		Real charLen;
		bool collideGpu;

		string extraOpts;
		int maxScheduledSteps;

		CLDEM_SERIALIZE_ATTRS((pNum)(dNum)/*binary*/(scene)(trackEnergy)(ktDivKn)(breakTension)(charLen)(collideGpu)(extraOpts)(maxScheduledSteps)/*other*/(par)(con)(conFree)(pot)(potFree)(cJournal)(clumps)(bboxes)(cpuCollider),/*otherCode*/);

		Simulation(int pNum=-1,int dNum=-1, bool _trackEnergy=false, Real _ktDivKn=NAN, bool _breakTension=false, Real _charLen=NAN, bool _collideGpu=false, const string& _opts=""): trackEnergy(_trackEnergy), ktDivKn(_ktDivKn), breakTension(_breakTension), charLen(_charLen), collideGpu(_collideGpu), extraOpts(_opts), maxScheduledSteps(-1) {}

		void save(const string& s){ ObjectIO::save(s,"cldem",*this); }
		void load(const string& s){ ObjectIO::load(s,"cldem",*this); }

		shared_ptr<cl::Platform> platform;
		shared_ptr<cl::Device> device;
		shared_ptr<cl::Context> context;
		shared_ptr<cl::CommandQueue> queue;
		shared_ptr<cl::Program> program;

		shared_ptr<CpuCollider> cpuCollider;

		/* read and write fixed-size objects (all async, use queue.finish() to wait) */
		template<typename T> cl::Buffer writeBuf(const T& obj,bool wait=false){
			cl::Buffer buf(*context,CL_MEM_READ_WRITE,sizeof(T),NULL);
			queue->enqueueWriteBuffer(buf,wait?CL_TRUE:CL_FALSE,0,sizeof(T),&obj);
			return buf;
		}
		template<typename T> void readBuf(cl::Buffer& buf, T& obj,bool wait=false){
			queue->enqueueReadBuffer(buf,wait?CL_TRUE:CL_FALSE,0,sizeof(T),&obj);
		}
		/* read and write std::vector containers (all async) */
		template<typename T> void writeVecBuf(std::vector<T>& obj, BufSize& bs){
			if(obj.empty()) throw std::runtime_error("Buffer created from std::vector<"+string(typeid(T).name())+"> may not be empty.");
			bs.buf=cl::Buffer(*context,CL_MEM_READ_WRITE,sizeof(T)*obj.size(),NULL);
			queue->enqueueWriteBuffer(bs.buf,CL_FALSE,0,sizeof(T)*obj.size(),obj.data());
			bs.size=obj.size();
		}
		template<typename T> void readVecBuf(std::vector<T>& obj, BufSize& bs){
			//cerr<<typeid(T).name()<<"["<<bs.size;
			obj.resize(bs.size); // in case size was changed meanwhile
			queue->enqueueReadBuffer(bs.buf,CL_FALSE,0,sizeof(T)*obj.size(),obj.data());
			//cerr<<"]";
		}

		size_t arrNewSize(const Scene&, int arrIx);


		void writeBufs(bool wait=true);
		void readBufs(bool wait=true);
		void runKernels(int nSteps, int substepStart=0);
		cl::Kernel makeKernel(const char* name);

		void initCl();
		void run(int nSteps);
		Real pWaveDt();
		py::tuple getBbox(par_id_t id);
		py::list saveVtk(string prefix, bool compress=true, bool ascii=false);
		int addClump(vector<Particle>&);

		Matrix3r inertiaTranslate(const Matrix3r& I, const Real m, const Vector3r& off){
			return I+m*(off.dot(off)*Matrix3r::Identity()-off*off.transpose());
		}
		Matrix3r inertiaRotate(const Matrix3r& I, const Matrix3r& T){
			return T.transpose()*I*T;
		}
	};
};

//static std::vector<Particle>* Simulation_par_get(Simulation* s){ return &(s->par); }
//static void Simulation_par_set(Simulation* s, const std::vector<Particle>& par){ s->par.vec=par; }

static
void Simulation_hpp_expose(){
	py::class_<clDem::Simulation,shared_ptr<clDem::Simulation>>("Simulation",py::init<int,int,bool,Real,bool,Real,bool,string>((py::arg("platformNum")=-1,py::arg("deviceNum")=-1,py::arg("trackEnergy")=false,py::arg("ktDivKn")=NAN,py::arg("breakTension")=false,py::arg("charLen")=NAN,py::arg("collideGpu")=false,py::arg("opts")="")))
		.PY_RW(Simulation,scene)
		.PY_RW(Simulation,par)
		.PY_RW(Simulation,con)
		.PY_RWV(Simulation,conFree)
		.PY_RWV(Simulation,pot)
		.PY_RWV(Simulation,potFree)
		.PY_RWV(Simulation,cJournal)
		.PY_RWV(Simulation,bboxes)
		.PY_RW(Simulation,maxScheduledSteps)
		.def_readonly("trackEnergy",&Simulation::trackEnergy)
		.def_readonly("ktDivKn",&Simulation::ktDivKn)
		.def_readonly("breakTension",&Simulation::breakTension)
		.def_readonly("charLen",&Simulation::charLen)
		.def_readwrite("collideGpu",&Simulation::collideGpu)
		.def("run",&Simulation::run,(py::arg("nSteps"),py::arg("resetArrays")=true))
		#ifdef CLDEM_VTK
		.def("saveVtk",&Simulation::saveVtk,(py::arg("prefix"),py::arg("compress")=true,py::arg("ascii")=false))
		#endif
		.def("pWaveDt",&Simulation::pWaveDt)
		.def("getBbox",&Simulation::getBbox)
		.def("addClump",&Simulation::addClump)
		.def("save",&Simulation::save)
		.def("load",&Simulation::load)
	;
}

