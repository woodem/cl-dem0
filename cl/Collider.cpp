#include"Collider.cl.h"
#include"Simulation.hpp"

#include<set>

//#define COLL_DBG(a) cerr<<a<<endl
#define COLL_DBG(a)

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
	scene->arrSize[ARR_POTFREE]=0;
	/* scene->arrAlloc is handled automatically by writeBufs */
}

// add new potential contact
// NB: does not check if the contact is not already present
void CpuCollider::addPot(par_id_t id1, par_id_t id2, bool useFree){
	// use the same logic as kernels to find a free slot in pot
	long ix=-1;
	if(useFree) ix=Scene_arr_fromFreeArr_or_append(scene,sim->potFree.data(),ARR_POTFREE,ARR_POT,/*shrink*/true);
	par_id2_t ids={id1,id2};
	if(ix<0){
		sim->pot.push_back(ids);
		scene->arrSize[ARR_POT]=sim->pot.size();
		ix=sim->pot.size()-1;
	}
	else{
		sim->pot[ix]=par_id2_t(ids);
	}
	//
	add(id1,id2,ConLoc(ix,/*isReal*/false));
	COLL_DBG("* new pot. ##"<<id1<<"+"<<id2<<" @ pot["<<ix<<"] ("<<scene->arrSize[ARR_POT]<<"/"<<scene->arrAlloc[ARR_POT]<<")");
}

void CpuCollider::delPot(par_id_t id1, par_id_t id2, ConLoc* cl){
	assert(find(id1,id2)==cl);
	assert(!cl->isReal);
	const par_id2_t& ids=sim->pot[cl->ix];
	assert(min(id1,id2)==min(ids.s0,ids.s1) && max(id1,id2)==max(ids.s0,ids.s1));
	// remove from pot
	par_id2_t no2={-1,-1};
	sim->pot[cl->ix]=no2;
	// add to potFree; this changes the element atomically to 0
	long ixPotFree=Scene_arr_findNegative_or_append(scene,ARR_POTFREE,sim->potFree.data());
	if(ixPotFree<0){
		sim->potFree.push_back(cl->ix); 
		scene->arrSize[ARR_POTFREE]=sim->potFree.size();
	} else {
		// 0 is when the slot was one which was previously unused (to avoid races)
		// -1 is 
		assert(sim->potFree[ixPotFree]==0 || sim->potFree[ixPotFree]==-1); 
		sim->potFree[ixPotFree]=cl->ix;
	}
	// remove from cMap
	remove(id1,id2);
	//
	COLL_DBG("* del pot ##"<<id1<<"+"<<id2<<" @ pot["<<cl->ix<<"] ("<<scene->arrSize[ARR_POT]<<"/"<<scene->arrAlloc[ARR_POT]<<")");
}

/* compact free arrays: move valid items to the beginning and
trailing -1's to the end
*/
void CpuCollider::compactFree(){
	for(int what=0; what<2; what++){
		std::vector<cl_int>& FF(what==0?sim->potFree:sim->conFree);
		#if 0
			cerr<<"$$$$$$ "<<what<<": ";
			for(int a: FF) cerr<<a<<" ";
			cerr<<endl;
		#endif
		long lastUsed=-1;
		for(long i=0; i<(long)FF.size(); i++){
			if(FF[i]<0) continue;
			// FF[i]>=0 now
			if(lastUsed<(i-1)){ // this item can be moved towards the front
				lastUsed++;
				FF[lastUsed]=FF[i];
				FF[i]=-1;
			} else { // if it cannot be moved, this is the last used index
				lastUsed=i;
			}
		}
		// set used array size
		sim->scene.arrSize[what==0?ARR_POTFREE:ARR_CONFREE]=lastUsed+1;
		#if 0
			cerr<<"@@@@@@ "<<what<<"["<<lastUsed+1<<"]: ";
			for(int a: FF) cerr<<a<<" ";
			cerr<<endl;
		#endif
	}
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
		cerr<<"Previous contact ##"<<c.ids.s0<<"+"<<c.ids.s1<<endl;
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
			AxBound bMin(id,mn,true,mn==mx), bMax(id,mx,false,mn==mx);
			bounds[ax][2*id]=bMin; bounds[ax][2*id+1]=bMax;
	if(2*id+1 == 1 || 2*id+1 ==5){
			std::cout << "L/R : " << 2*id << " / " << 2*id+1 << std::endl;
			std::cout << "bmax: " << bMax.coord << std::endl;
			std::cout << "mx: " << mx << std::endl;
	}
		}
	}


for (int i = 0; i < 2*N; i++){
	if(bounds[0][i].coord < 0.6){

	} else {
		std::cout << "tu: " << i << " coord: " << bounds[0][i].coord << std::endl;
	}
}

#if 1
	std::cout << "GPU" << std::endl;
	cl::Buffer boundBufs[3];

	int powerOfTwo = (pow(2, trunc(log2(2*N))) == 2*N) ?
		pow(2, trunc(log2(2*N))) : pow(2, trunc(log2(2*N)) + 1);

	cl_uint bits = 0;
	for (int temp = powerOfTwo; temp > 1; temp >>= 1) {
		++bits;
	}

	int less = powerOfTwo - 2*N;
	cerr << "=======================" << endl;
	cerr << "powerOfTwo : " << powerOfTwo << endl;
	cerr << "less : " << less << endl;
	cerr << "=======================" << endl;

	int local_size = 16;
	int global_size = (trunc(trunc(sqrt(powerOfTwo / 2)) / local_size) +  1) * local_size;
	//prepare size of bounds for bittonic-sort
	std::cout << "priprava bufferu" << std::endl;
	for (int ax:{0, 1, 2}){
		std::cout << "1 ax: " << ax << std::endl;
		boundBufs[ax] = cl::Buffer(sim->context, CL_MEM_READ_WRITE,	powerOfTwo * sizeof (AxBound), NULL);
			std::cout << "2" << std::endl;
		bounds[ax].resize(powerOfTwo);
			std::cout << "3" << std::endl;	
		sim->queue.enqueueWriteBuffer(boundBufs[ax], CL_TRUE, 0, powerOfTwo * sizeof (AxBound), bounds[ax].data());
			std::cout << "4" << std::endl;
	}


	std::cout << "create kernel" << std::endl;
	try {
		cl::Kernel sortKernel(sim->program, "sortBitonic");
		cl::Event eve;
		sortKernel.setArg(3, powerOfTwo);
		sortKernel.setArg(4, 1);
		
		for (int ax:{0, 1, 2}){
			sortKernel.setArg(0, boundBufs[ax]);
			for (cl_uint stage = 0; stage < bits; stage++) {
				sortKernel.setArg(1, stage);
				for (cl_uint passOfStage = 0; passOfStage < stage + 1; passOfStage++) {
					sortKernel.setArg(2, passOfStage);
					sim->queue.enqueueNDRangeKernel(sortKernel, cl::NullRange,
						cl::NDRange(global_size, global_size),
						cl::NDRange(local_size, local_size), NULL, &eve);
					eve.wait();
				}
			}
		}

	} catch (cl::Error& e){
		std::cerr << e.what() << e.err() << std::endl;			
		throw;
	}
	std::cout << "sort OK" << std::endl;

	for (int ax:{0, 1, 2}){
		bounds[ax].resize(2*N);
		sim->queue.enqueueReadBuffer(boundBufs[ax], CL_TRUE, 0, 2*N*sizeof (AxBound), bounds[ax].data());
		boundBufs[ax] = cl::Buffer(sim->context, CL_MEM_READ_WRITE, 2*N*sizeof (AxBound), NULL);
		sim->queue.enqueueWriteBuffer(boundBufs[ax], CL_TRUE, 0, 2*N*sizeof (AxBound), bounds[ax].data());
	}
	std::cout << "shrink buffer OK" << std::endl;
	
	/*for (int i = 0; i < 2*N; i++){
		if(bounds[0][i].coord > bounds[0][i+1].coord){
			std::cout << bounds[0][i].coord << std::endl;
		}
	}*/

	std::vector<AxBound> bounds1[3];
	
	bounds1[0] = bounds[0];
	bounds1[1] = bounds[1];
	bounds1[2] = bounds[2];


//	for(int ax:{0,1,2}){
//		std::sort(bounds[ax].begin(),bounds[ax].end());
//	}

	std::cout << "test sorted " << std::endl;


	for(int ax:{0,1,2}){
	for(int i = 0;  i < 2*N; i++) {
	//std::cout << "id: " << i << " : " << bounds[0][i].coord << " x " <<bounds[1][i].coord << " x " << bounds[2][i].coord << std::endl;
	/*if(bounds[ax][i].getId() != bounds1[ax][i].getId()){
		std::cout << "error: " << ax << "coord cpu/gpu" << bounds[ax][i].coord << "/" 
			<< bounds1[ax][i].coord << std::endl;
	}*/
	}
	}
	std::cout << "test sorted OK : " << N << std::endl;

	const int ax0=0; // traverse along x, for example
	for(size_t i=0; i<2*N; i++){
		const AxBound& b0=bounds[ax0][i];
	//std::cout << "A" << std::endl;
		if(!b0.isMin() || b0.isThin()) continue;
		for(size_t j=i+1; bounds[ax0][j].getId()!=b0.getId(); j++){
	//std::cout << "B" << std::endl;
	//std::cout << j << std::endl;
			par_id_t id1=b0.getId(),id2=bounds[ax0][j].getId();
			// cerr<<"##"<<id1<<"+"<<id2<<endl;
			if(!bboxOverlap(id1,id2)) continue;
	//std::cout << "C" << std::endl;
			if(!Scene_particles_may_collide(scene,&(sim->par[id1]),&(sim->par[id2]))) continue;
	//std::cout << "D" << std::endl;
			if(find(id1,id2)) continue;
	//std::cout << "E" << std::endl;
			if(id1>id2) std::swap(id1,id2);
	//std::cout << "F" << std::endl;
			addPot(id1,id2, /* useFree */ false);

		}
	}

#else

	// sort initial bounds
	for(int ax:{0,1,2}){
		std::sort(bounds[ax].begin(),bounds[ax].end());
		// for(const auto& b: bounds[ax]) cerr<<ax<<" "<<b.id<<" "<<(b.isMin?"<-":"->")<<" "<<b.coord<<" "<<(b.isThin?"THIN":"")<<endl;
	}

	std::cout << bounds[0][2*N-2].coord << "/" << bounds[0][2*N-1].coord << "/" << bounds[0][2*N].coord << std::endl;

	// create potential contacts
	const int ax0=0; // traverse along x, for example
	for(size_t i=0; i<2*N; i++){
		const AxBound& b0=bounds[ax0][i];
		if(!b0.isMin() || b0.isThin()) continue;
		for(size_t j=i+1; bounds[ax0][j].getId()!=b0.getId(); j++){
			//std::cout << j << std::endl;
			par_id_t id1=b0.getId(),id2=bounds[ax0][j].getId();
			// cerr<<"##"<<id1<<"+"<<id2<<endl;
			if(!bboxOverlap(id1,id2)) continue;
			if(!Scene_particles_may_collide(scene,&(sim->par[id1]),&(sim->par[id2]))) continue;
			if(find(id1,id2)) continue;
			if(id1>id2) std::swap(id1,id2);
			addPot(id1,id2, /* useFree */ false);
		}
	}

#endif

	std::cout << "GPU code is OK" << std::endl;

	if(sim->pot.empty()){
		par_id2_t no2={-1,-1};
		sim->pot.push_back(no2);
		scene->arrSize[ARR_POT]=0;
		scene->arrAlloc[ARR_POT]=1;
	} else {
		scene->arrSize[ARR_POT]=scene->arrAlloc[ARR_POT]=sim->pot.size();
	}
	/* potFree was filled with -1, but not deallocated */
}

void CpuCollider::incrementalStep(){
	//throw std::runtime_error("CpuCollider::incrementalStep not yet implemented!");
	replayJournal();
	checkConsistency();
	updateBounds();
	insertionSort();
	compactFree();
}

/* reply all contact changes and update cMap accordingly */
void CpuCollider::replayJournal(){
	for(size_t i=0; i<scene->arrSize[ARR_CJOURNAL]; i++){
		assert(i<sim->cJournal.size());
		CJournalItem& J(sim->cJournal[i]);
		std::string cStr=ids2str(J.ids);
		ConLoc* cl=find(J.ids.s0,J.ids.s1);
		if(!cl) throw std::runtime_error("Journal: referencing non-existent contact "+cStr);
		switch(J.what){
			case CJOURNAL_POT2CON:
				if(cl->isReal) throw std::runtime_error("Journal: making already real "+cStr+" real");
				cl->ix=J.index;
				cl->isReal=true;
				break;
			case CJOURNAL_CON2POT:
				if(!cl->isReal) throw std::runtime_error("Journal: making already porential contact "+cStr+" potential");
				cl->ix=J.index;
				cl->isReal=false;
				break;
			case CJOURNAL_CON_DEL:
				if(!cl->isReal) throw std::runtime_error("Journal: potential contact "+cStr+" being deleted outside the collider");
				remove(J.ids.s0,J.ids.s1);
				break;
			default: abort();
		}
		J=CJournalItem();
	}
	scene->arrSize[ARR_CJOURNAL]=0;
};
/* once the journal has been updated, we can check that cMap corresponds to arrays within the simulation:

1. all contacts in cMap are in con/pot
2. all valid contacts in con/pot are in cMap
3. all indices in conFree/potFree are really free (invalid) in con/pot
4. all invalid items in con/pot are in conFree/potFree (possibly expensive)

TODO: check array sizes before dereferencing their elements

*/
void CpuCollider::checkConsistency(){
	// all contacts in cMap are in con/pot
	for(long id1=0; id1<(long)cMap.size(); id1++){
		for(const auto& id2_cl: cMap[id1]){
			const long& id2=id2_cl.first;
			const ConLoc& cl=id2_cl.second;
			par_id2_t s2=(cl.isReal?sim->con[cl.ix].ids:sim->pot[cl.ix]);
			if(min(id1,id2)!=min(s2.s0,s2.s1) || max(id1,id2)!=max(s2.s0,s2.s1)) throw std::runtime_error("Inconsistency: "+string(cl.isReal?"con":"pot")+"["+lexical_cast<string>(cl.ix)+"] should be "+ids2str(id1,id2)+", is "+ids2str(s2));
		}
	}
	// all valid contacts in con/pot are in cMap
	// also check that unused but allocated space is occupied by invalid contacts
	for(size_t i=0; i<sim->con.size(); i++){
		const Contact& c=sim->con[i];
		string conStr="con["+lexical_cast<string>(i)+"]";
		if(i>=scene->arrSize[ARR_CON]){ // allocated but unused contact
			if(c.ids.s0>=0 || c.ids.s1>=0) throw std::runtime_error("Inconsistency: "+conStr+" is unused and should have ids (-1,-1), has "+ids2str(c.ids));
		} else { // contact in the used part of the array
			if(min(c.ids.s0,c.ids.s1)<0 && max(c.ids.s0,c.ids.s1)>=0) throw std::runtime_error("Inconsistency: "+conStr+" has one id negative (must be both or none): "+ids2str(c.ids));
			// invalid contact: make sure it does not exist
			ConLoc* cl=find(c.ids.s0,c.ids.s1);
			if(!cl && c.ids.s0>=0) throw std::runtime_error("Inconsistency: "+conStr+" references "+ids2str(c.ids)+" which the collider knows nothing about");
			if(cl && !cl->isReal) throw std::runtime_error("Inconsistency: "+conStr+" is real, but the collider thinks it is potential");
		}
	}
	for(size_t i=0; i<sim->pot.size(); i++){
		const par_id2_t& ids=sim->pot[i];
		string potStr="pot["+lexical_cast<string>(i)+"]";
		if(i>=scene->arrSize[ARR_POT]){
			if(ids.s0>=0 || ids.s1>=0) throw std::runtime_error("Inconsistency: "+potStr+" is unused and should have ids (-1,-1), has "+ids2str(ids));
		} else {
			if(min(ids.s0,ids.s1)<0 && max(ids.s0,ids.s1)>=0) throw std::runtime_error("Inconsistency: "+potStr+" has one id negative (must be both or none): "+ids2str(ids));
			ConLoc* cl=find(ids.s0,ids.s1);
			if(!cl && ids.s0>=0) throw std::runtime_error("Inconsistency: "+potStr+" references "+ids2str(ids)+" which the collider knows nothing about");
			if(cl && cl->isReal) throw std::runtime_error("Inconsistency: "+potStr+" is potential, but the collider things it is real");
		}
	}
	
	// all indices in conFree/potFree are really free (invalid) in con/pot
	for(size_t i=0; i<sim->conFree.size(); i++){
		string cfStr("conFree["+lexical_cast<string>(i)+"]");
		long cf=sim->conFree[i];
		if(i>=scene->arrSize[ARR_CONFREE]){
			if(cf>=0) throw std::runtime_error("Inconsistency: "+cfStr+" is unused ("+lexical_cast<string>(scene->arrSize[ARR_CONFREE])+"/"+lexical_cast<string>(scene->arrAlloc[ARR_CONFREE])+") and should be negative, not "+lexical_cast<string>(cf));
		} else {
			if(cf>=0 && sim->con[cf].ids.s0>=0) throw std::runtime_error("Inconsistency: "+cfStr+" says con["+lexical_cast<string>(cf)+"] is free, but is "+ids2str(sim->con[cf].ids));
		}
	}
	for(size_t i=0; i<sim->potFree.size(); i++){
		string pfStr("potFree["+lexical_cast<string>(i)+"]");
		long pf=sim->potFree[i];
		if(i>=scene->arrSize[ARR_POTFREE]){
			if(pf>=0) throw std::runtime_error("Inconsistency: "+pfStr+" is unused ("+lexical_cast<string>(scene->arrSize[ARR_POTFREE])+"/"+lexical_cast<string>(scene->arrAlloc[ARR_POTFREE])+") and should be negative, not "+lexical_cast<string>(pf));
		} else {
			if(pf>=0 && sim->pot[pf].s0>=0) throw std::runtime_error("Inconsistency: "+pfStr+" says pot["+lexical_cast<string>(pf)+"] is free, but is "+ids2str(sim->pot[pf]));
		}
	}

	// all invalid items in con/pot are in conFree/potFree (possibly expensive)
	// create sets for fast lookup
	std::set<cl_long> conFree(sim->conFree.begin(),sim->conFree.end()), potFree(sim->potFree.begin(),sim->potFree.end());
	for(size_t i=0; i<scene->arrSize[ARR_CON]; i++){
		const Contact& c=sim->con[i];
		if(c.ids.s0<0 && conFree.count(i)==0) throw std::runtime_error("Inconsistency: con["+lexical_cast<string>(i)+"] is free, but not in conFree");
		if(c.ids.s0>=0 && conFree.count(i)>0) throw std::runtime_error("Inconsistency: con["+lexical_cast<string>(i)+"] is "+ids2str(c.ids)+", but is in conFree");
	}
	for(size_t i=0; i<scene->arrSize[ARR_POT]; i++){
		const par_id2_t& ids=sim->pot[i];
		if(ids.s0<0 && potFree.count(i)==0) throw std::runtime_error("Inconsistency: pot["+lexical_cast<string>(i)+"] is free, but not in potFree");
		if(ids.s0>=0 && potFree.count(i)>0) throw std::runtime_error("Inconsistency: pot["+lexical_cast<string>(i)+"] is "+ids2str(ids)+", but is in potFree");
	}
}

void CpuCollider::updateBounds(){
	for(int ax:{0,1,2}){
		for(AxBound& ab: bounds[ax]){
			ab.coord=sim->bboxes[6*ab.getId()+ax+(ab.isMin()?0:3)];
			//ab.isThin=(sim->bboxes[6*ab.id+ax]==sim->bboxes[6*ab.id+ax+3]); // this is perhaps not needed anymore
		}
	}
}

void CpuCollider::insertionSort(){
	for(int ax:{0,1,2}){
		vector<AxBound>& bb(bounds[ax]);
		long iMax=bb.size();
		for(long i=0; i<iMax; i++){
			const AxBound bbInit=bb[i]; // copy, so that it is const
			if(isnan(bbInit.coord)) continue;
			long j=i-1;
			while(j>=0 && (bb[j].coord>bbInit.coord || isnan(bb[j].coord))){
				// bbInit is bb[j+1] which is travelling downwards and swaps with b[j]
				// do not store min-min, max-max, nan inversions
				if(bbInit.isMin()!=bb[j].isMin() && !isnan(bb[j].coord)){
					int minId=min(bbInit.getId(),bb[j].getId()), maxId=max(bbInit.getId(),bb[j].getId());
					#if 0
						// min going below max
						if(bbInit.isMin) inv.push_back(Vector2i(minId,maxId));
						// max going below min
						else inv.push_back(Vector2i(maxId,minId));
					#else
						// instead of writing to the inversions array, handle contacts right here
						ConLoc* cl=find(minId,maxId);

						// min going below max; if the contact is potential, delete it
						// (no need to check bboxes, they cannot possibly overlap anymore)
						if(!bbInit.isMin() && cl && !cl->isReal){
							assert(!bboxOverlap(minId,maxId));
							delPot(minId,maxId,cl);
						}
						// min going below max: the contact might be created, if there is overlap
						if(bbInit.isMin() && !cl && bboxOverlap(minId,maxId)){
							if(Scene_particles_may_collide(scene,&(sim->par[minId]),&(sim->par[maxId]))) addPot(minId,maxId,/*useFree*/true);
						}
						// otherwise, the contact is real or non-existent, do nothing
					#endif
					// COLL_DBG("["<<maxId<<"↔ "<<minId<<"]");
				}
				bb[j+1]=bb[j];
				j--;
			}
			bb[j+1]=bbInit;
		}
	}
}

};
