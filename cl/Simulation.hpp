#pragma once

#include"Particle.cl.h"
#include"Contact.cl.h"
#include"Scene.cl.h"
#include"Collider.cl.h"

// if included from yade
#ifdef YADE_VTK
	#define CLDEM_VTK
#endif

#if 0
#ifndef YADE_CLDEM
	#include<boost/archive/binary_oarchive.hpp>
	#include<boost/archive/binary_iarchive.hpp>
	// after supported archive types
	#include<boost/serialization/export.hpp>
	#include<boost/serialization/vector.hpp>
#endif

namespace boost {
	namespace align { struct __attribute__((__aligned__(128))) a128 {};}
	template<> class type_with_alignment<128> { public: typedef align::a128 type; };
};

#define _BITWISE_PRIMITIVE(klass) BOOST_IS_BITWISE_SERIALIZABLE(klass);  namespace boost{ namespace serialization{ template<class Archive> void serialize(Archive & ar, klass & obj, const unsigned int version){ throw std::runtime_error("Class " #klass " must be serialized bitwise, not by calling the serialization function (are you serializing to XML?)."); } }}
#define _BITWISE(klass) _BITWISE_PRIMITIVE(klass); BOOST_CLASS_TRACKING(klass,boost::serialization::track_never)
	_BITWISE(clDem::Scene)
	_BITWISE(clDem::Particle)
	_BITWISE(clDem::Contact)
	_BITWISE(clDem::CJournalItem)
	_BITWISE_PRIMITIVE(cl_long2)
	_BITWISE_PRIMITIVE(clDem::par_id_t)
#undef _BITWISE
#undef _BITWISE_PRIMITIVE
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
		vector<cl_long2> pot; // potential contacts (only the id1,id2-tuple)
		vector<cl_int> potFree; // free slots in pot
		vector<CJournalItem> cJournal; // logging changes in contact arrays so that the collider's internal structures can be updated accordingly
		vector<par_id_t> clumps;
		vector<cl_float> bboxes;
			
		#if 0
		friend class boost::serialization::access;
		template<class ArchiveT> void serialize(ArchiveT & ar, unsigned int version){
			ar & BOOST_SERIALIZATION_NVP(scene);
			ar & BOOST_SERIALIZATION_NVP(par);
			ar & BOOST_SERIALIZATION_NVP(con);
			ar & BOOST_SERIALIZATION_NVP(conFree);
			ar & BOOST_SERIALIZATION_NVP(pot);
			ar & BOOST_SERIALIZATION_NVP(potFree);
			ar & BOOST_SERIALIZATION_NVP(clumps);
			ar & BOOST_SERIALIZATION_NVP(bboxes);
		}
		#endif

		Simulation(int pNum=-1,int dNum=-1, const string& opts=""){
			scene=Scene();
			maxScheduledSteps=-1;
			initCl(pNum,dNum,opts);
			cpuCollider=make_shared<CpuCollider>();
		}

		cl::Platform platform;
		cl::Device device;
		cl::Context context;
		cl::CommandQueue queue;
		cl::Program program;

		shared_ptr<CpuCollider> cpuCollider;

		/* read and write fixed-size objects (all async, use queue.finish() to wait) */
		template<typename T> cl::Buffer writeBuf(const T& obj,bool wait=false){
			cl::Buffer buf(context,CL_MEM_READ_WRITE,sizeof(T),NULL);
			queue.enqueueWriteBuffer(buf,wait?CL_TRUE:CL_FALSE,0,sizeof(T),&obj);
			return buf;
		}
		template<typename T> void readBuf(cl::Buffer& buf, T& obj,bool wait=false){
			queue.enqueueReadBuffer(buf,wait?CL_TRUE:CL_FALSE,0,sizeof(T),&obj);
		}
		/* read and write std::vector containers (all async) */
		template<typename T> void writeVecBuf(std::vector<T>& obj, BufSize& bs){
			if(obj.empty()) throw std::runtime_error("Buffer created from std::vector<"+string(typeid(T).name())+"> may not be empty.");
			bs.buf=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(T)*obj.size(),NULL);
			queue.enqueueWriteBuffer(bs.buf,CL_FALSE,0,sizeof(T)*obj.size(),obj.data());
			bs.size=obj.size();
		}
		template<typename T> void readVecBuf(std::vector<T>& obj, BufSize& bs){
			//cerr<<typeid(T).name()<<"["<<bs.size;
			obj.resize(bs.size); // in case size was changed meanwhile
			queue.enqueueReadBuffer(bs.buf,CL_FALSE,0,sizeof(T)*obj.size(),obj.data());
			//cerr<<"]";
		}

		size_t arrNewSize(const Scene&, int arrIx);


		void writeBufs(bool wait=true);
		void readBufs(bool wait=true);
		void runKernels(int nSteps, int substepStart=0);
		cl::Kernel makeKernel(const char* name);

		int maxScheduledSteps;


		void initCl(int pNum, int dNum, const string& opts);
		void run(int nSteps);
		Real pWaveDt();
		py::tuple getBbox(par_id_t id);
		py::list saveVtk(string prefix, bool compress=true, bool ascii=false);
	};
};

//static std::vector<Particle>* Simulation_par_get(Simulation* s){ return &(s->par); }
//static void Simulation_par_set(Simulation* s, const std::vector<Particle>& par){ s->par.vec=par; }

static
void Simulation_hpp_expose(){
	py::class_<clDem::Simulation,shared_ptr<clDem::Simulation>>("Simulation",py::init<int,int,string>((py::arg("platformNum")=-1,py::arg("deviceNum")=-1,py::arg("opts")="")))
		.PY_RW(Simulation,scene)
		.PY_RW(Simulation,par)
		.PY_RW(Simulation,con)
		.PY_RWV(Simulation,conFree)
		.PY_RWV(Simulation,pot)
		.PY_RWV(Simulation,potFree)
		.PY_RWV(Simulation,cJournal)
		.PY_RWV(Simulation,bboxes)
		.PY_RW(Simulation,maxScheduledSteps)
		.def("run",&Simulation::run,(py::arg("nSteps"),py::arg("resetArrays")=true))
		#ifdef CLDEM_VTK
		.def("saveVtk",&Simulation::saveVtk,(py::arg("prefix"),py::arg("compress")=true,py::arg("ascii")=false))
		#endif
		.def("pWaveDt",&Simulation::pWaveDt)
		.def("getBbox",&Simulation::getBbox)
	;
}

