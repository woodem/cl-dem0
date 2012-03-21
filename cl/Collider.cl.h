#ifndef _COLLIDER_CL_H_
#define _COLLIDER_CL_H_

#include"common.cl.h"
#include"serialization.cl.h"

CLDEM_NAMESPACE_BEGIN()

// values of CJournalItem::what
enum _cjournal {
	CJOURNAL_POT2CON=0, // potential contact becomes real, deleted from pot and added to con[index]
	CJOURNAL_CON2POT,   // real contact broken & bboxes overlap; deleted from con and added to pot[index]
	CJOURNAL_CON_DEL    // real contact broken & no bboxes overlap; deleted from con (index is unused)
};

struct CJournalItem;
inline void CJournalItem_init(global struct CJournalItem*);

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
		CJournalItem(){ CJournalItem_init(this); }
	#endif
	CLDEM_SERIALIZE_ATTRS((what)(index)(ids),/**/);
};

inline void CJournalItem_init(global struct CJournalItem* i){
	i->what=i->index=i->ids.s0=i->ids.s1=-1;
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
			CLDEM_SERIALIZE_ATTRS((ix)(isReal),/**/);
		};
		typedef std::vector<std::map<par_id_t,ConLoc>> cMapT;
		cMapT cMap;

		struct AxBound{
			//data
			float coord;
			par_id_t idMinThin;
			AxBound(){}
			AxBound(par_id_t id, float _coord, bool _isMin, bool _isThin): idMinThin(id<<2 | (_isMin?1:0) | (_isThin?2:0)), coord(_coord) {}
			//AxBound(const AxBound& b): coord(b.coord), idMinThin(b.idMinThin){}
			//AxBound& operator=(const AxBound& b){ coord=b.coord; idMinThin=b.idMinThin; return *this; }
			bool operator<(const AxBound& b) const { return coord<b.coord; } // simply sort by coord; isThin takes care of mn==mx
			bool isMin() const { return idMinThin&1; }
			bool isThin() const{ return idMinThin&2; }
			par_id_t getId() const { return idMinThin>>2; }
			CLDEM_SERIALIZE_ATTRS((idMinThin)(coord),);
		};
		std::vector<AxBound> bounds[3];

		CLDEM_SERIALIZE_ATTRS((cMap)(bounds),/**/);

		void add(par_id_t id1, par_id_t id2, const ConLoc&);
		void remove(par_id_t id1, par_id_t id2);
		ConLoc* find(par_id_t id1, par_id_t id2);

		bool bboxOverlap(par_id_t id1, par_id_t id2) const;

		void clearSimPot();
		void addPot(par_id_t id1, par_id_t id2, bool useFree /*use potFree array */);
		void delPot(par_id_t id1, par_id_t id2, ConLoc* cl);
		void compactFree();

		std::string ids2str(par_id_t id1, par_id_t id2){ return "##"+lexical_cast<string>(id1)+"+"+lexical_cast<string>(id2); }
		std::string ids2str(const par_id2_t& ii){ return ids2str(ii.s0,ii.s1); }

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
