#include"Collider.cl.h"
#include"Simulation.hpp"

namespace clDem{

CpuCollider::ConLoc* CpuCollider::find(par_id_t id1, par_id_t id2){
	// input out of range
	if(id1<0 || id1>=cMap.size() || id2<0 || id2>=cMap.size()) return NULL;
	if(id1>id2) std::swap(id1,id2);
	auto iter=cMap[id1].find(id2);
	if(iter!=cMap[id1].end()) return &(iter->second);
	return NULL;
};

void CpuCollider::add(par_id_t id1, par_id_t id2, const ConLoc& cl){
	if(id1>id2) std::swap(id1,id2);
	if(id1<0 || id1>cMap.size() || id2<0 || id2>=cMap.size()) throw std::runtime_error("Adding ConLoc with out-of-range ids ("+lexical_cast<string>(id1)+","+lexical_cast<string>(id2)+").");
	if(cMap[id1].count(id2)>0) throw std::runtime_error("Adding ConLoc with already-existing ids ("+lexical_cast<string>(id1)+","+lexical_cast<string>(id2)+").");
	cMap[id1][id2]=cl;
};

void CpuCollider::remove(par_id_t id1, par_id_t id2){
	if(id1>id2) std::swap(id1,id2);
	if(find(id1,id2)==NULL) throw std::runtime_error("Deleting ConLon with non-existent ids ("+lexical_cast<string>(id1)+","+lexical_cast<string>(id2)+").");
	cMap[id1].erase(id2);
};

bool CpuCollider::bboxOverlap(par_id_t id1, par_id_t id2) const {
	return
		sim->bboxes[6*id1+0]<=sim->bboxes[6*id2+3] && sim->bboxes[6*id2+0]<=sim->bboxes[6*id1+3] &&
		sim->bboxes[6*id1+1]<=sim->bboxes[6*id2+4] && sim->bboxes[6*id2+1]<=sim->bboxes[6*id1+4] &&
		sim->bboxes[6*id1+2]<=sim->bboxes[6*id2+5] && sim->bboxes[6*id2+2]<=sim->bboxes[6*id1+5];
}

// clear sim->pot (all potential contacts)
void CpuCollider::clearSimPot(){
	sim->pot.clear();
	// potFree is filled with invalid items, to avoid reallocation during simulation later
	std::fill(sim->potFree.begin(),sim->potFree.end(),-1);
	/* scene->arrSize & scene->arrAlloc are handled automatically by writeBufs(setArrays=true) */
}

// add new potential contact
// NB: does not check if the contact is not already present
size_t CpuCollider::addPot(par_id_t id1, par_id_t id2){
	cl_long2 ids={id1,id2};
	sim->pot.push_back(ids);
	return sim->pot.size()-1; // return index of the new item
}


void CpuCollider::run(Simulation* _sim){
	sim=_sim;
	scene=&sim->scene;
	bool init=false; // run the initial (full) contact detection, or only the incremental one?
	if(bounds[0].size()!=sim->par.size()*2) init=true;
	if(init) initialStep();
	else incrementalStep();
};

void CpuCollider::initialStep(){
	size_t N=sim->par.size();
	cMap.clear(); cMap.resize(N);

	// copy real CL contact to cMap (if any)
	for(size_t i=0; i<scene->arrSize[ARR_CON]; i++){
		const Contact& c=sim->con[i];
		if(c.ids.s0<0) continue;
		add(c.ids.s0,c.ids.s1,ConLoc(i,/*isReal*/true));
	}
	clearSimPot();

	// copy bboxes to unsorted AxBound arrays
	for(int ax:{0,1,2}){
		bounds[ax].resize(2*N);
		for(size_t id=0; id<N; id++){
			float mn=sim->bboxes[6*id+ax], mx=sim->bboxes[6*id+3+ax];
			// add NaN as thin bbox 0--0
			if(isnan(mn) || isnan(mx)) mn=mx=0;
			if(!(mn<=mx)) throw std::runtime_error("#"+lexical_cast<string>(id)+", axis "+lexical_cast<string>(ax)+": min>max, "+lexical_cast<string>(mn)+">"+lexical_cast<string>(mx));
			AxBound bMin={id,mn,true,mn==mx}, bMax={id,mx,false,mn==mx};
			bounds[ax][2*id]=bMin; bounds[ax][2*id+1]=bMax;
		}
	}
	// sort initial bounds
	for(int ax:{0,1,2}){
		std::sort(bounds[ax].begin(),bounds[ax].end());
		// for(const auto& b: bounds[ax]) cerr<<ax<<" "<<b.id<<" "<<(b.isMin?"<-":"->")<<" "<<b.coord<<" "<<(b.isThin?"THIN":"")<<endl;
	}
	// create potential contacts
	const int ax0=0; // traverse along x, for example
	for(size_t i=0; i<2*N; i++){
		const AxBound& b0=bounds[ax0][i];
		if(!b0.isMin || b0.isThin) continue;
		for(size_t j=i+1; bounds[ax0][j].id!=b0.id; j++){
			par_id_t id1=b0.id,id2=bounds[ax0][j].id;
			// cerr<<"##"<<id1<<"+"<<id2<<endl;
			if(!bboxOverlap(id1,id2)) continue;
			if(find(id1,id2)) continue;
			if(id1>id2) std::swap(id1,id2);
			size_t index=addPot(id1,id2);
			add(id1,id2,ConLoc(index,/*isReal*/false));
			cerr<<"* new pot. ##"<<id1<<"+"<<id2<<" @ pot["<<index<<"]"<<endl;
		}
	}
}

void CpuCollider::incrementalStep(){
	throw std::runtime_error("CpuCollider::incrementalStep not yet implemented!");
}



};
