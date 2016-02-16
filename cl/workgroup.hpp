#include"common.cl.h"

cl::NDRange makeLinearNDRange(size_t i, const shared_ptr<cl::Device> dev){
	assert(dev->getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>()>=3);
	vector<size_t> wis=dev->getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
	assert(wis.size()>=3);
	wis.resize(3); // in case it is greater than 3
	size_t iMax=wis[0]*wis[1]*wis[2];
	if(i>iMax) throw std::runtime_error("makeLinearNDRange: required number of items ("+lexical_cast<string>(i)+" is bigger than the device's maximum "+lexical_cast<string>(iMax)+" ("+lexical_cast<string>(wis[0])+"×"+lexical_cast<string>(wis[1])+"×"+lexical_cast<string>(wis[2])+")");
	if(i<=wis[2]) return cl::NDRange(1,1,i);
	if(i<=wis[1]*wis[2]) return cl::NDRange(1,i/wis[2],wis[2]);
	return cl::NDRange(i/(wis[1]*wis[2]),wis[1],wis[2]);
}
