#ifndef _COLLIDER_CL_H_
#define _COLLIDER_CL_H_

#include"common.cl.h"

CLDEM_NAMESPACE_BEGIN()

// values of CJournalItem::what
enum _cjournal {
	CJOURNAL_POT2CON=0, // potential contact becomes real, deleted from pot and added to con[index]
	CJOURNAL_CON2POT,   // real contact broken & bboxes overlap; deleted from con and added to pot[index]
	CJOURNAL_CON_DEL    // real contact broken & no bboxes overlap; deleted from con (index is unused)
};

/* structure for logging changes in contact arrays */
struct CJournalItem {
	// type of change; (-1 means an invalid item which should be skipped)
	cl_int what;
	// new index within array, if applicable
	cl_int index;
	// particles in contact
	par_id2_t ids;
	#ifdef __cplusplus
		std::string __str__(){
			if(ids.s0<0) return "{}";
			return "{##"+lexical_cast<string>(ids.s0)+"+"+lexical_cast<string>(ids.s1)+" "+(what==CJOURNAL_POT2CON?("→ con["+lexical_cast<string>(index)+"]"):(what==CJOURNAL_CON2POT?("⇒ pot["+lexical_cast<string>(index)+"]"):"×"))+"}";
		}
	#endif
};

static struct CJournalItem CJournalItem_new(){
	struct CJournalItem ret;
	ret.what=ret.index=ret.ids.s0=ret.ids.s1=-1;
	return ret;
};


CLDEM_NAMESPACE_END()

// define to make collider check that changes in cLog
// are consistent with what is actually in con & pot arrays
#define CJOURNAL_CHECK_CONSISTENCY

#ifdef __cplusplus
namespace clDem{
	struct Simulation;
	struct Scene;
	struct CpuCollider{
		// index in either con (if isReal) or pot (if !isReal)
		struct ConLoc{
			ConLoc(size_t _ix=-1, bool _isReal=false): ix(_ix), isReal(_isReal) {};
			//ConLoc(const ConLoc& b){ *this=b; }
			//void operator=(const ConLoc& b){ ix=b.ix; isReal=b.isReal; } 
			long ix; bool isReal;
		};
		typedef std::vector<std::map<par_id_t,ConLoc>> cMapT;
		cMapT cMap;

		struct AxBound{
			// simply sort by coord; isThin takes care of mn==mx
			bool operator<(const AxBound& b) const { return coord<b.coord; }
			par_id_t id;
			float coord;
			bool isMin: 1;
			bool isThin:1;
			// AxBound(par_id_t _id, float _coord, bool _isMin, bool _isThin);
		};
		std::vector<AxBound> bounds[3];

		void add(par_id_t id1, par_id_t id2, const ConLoc&);
		void remove(par_id_t id1, par_id_t id2);
		ConLoc* find(par_id_t id1, par_id_t id2);

		bool bboxOverlap(par_id_t id1, par_id_t id2) const;

		void clearSimPot();
		void addPot(par_id_t id1, par_id_t id2, bool useFree /*use potFree array */);
		void delPot(par_id_t id1, par_id_t id2, ConLoc* cl);

		std::string ids2str(par_id_t id1, par_id_t id2){ return "##"+lexical_cast<string>(id1)+"+"+lexical_cast<string>(id2); }
		std::string ids2str(const cl_long2& ii){ return ids2str(ii.s0,ii.s1); }

		void run(Simulation*);
		Simulation* sim;
		Scene* scene;

		void initialStep();
		void incrementalStep();

		void replayJournal();
		void checkConsistency();
		void updateBounds();
		void insertionSort();
	};
};

static void Collider_cl_h_expose(){
	VECTOR_SEQ_CONV(CJournalItem);
	py::class_<CJournalItem>("CJournalItem")
		.def("__str__",&CJournalItem::__str__)
		.def("__repr__",&CJournalItem::__str__)
	;
	//PY_RW(CJournalItem,what).PY_RW(CJournalItem,index).PY_RW(CJournalItem,index).PY

};

#endif

#endif
