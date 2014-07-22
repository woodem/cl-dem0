#ifndef _COMMON_CL_H
#define _COMMON_CL_H

#include"../../cl-math/basic-math.cl"
#ifdef __cplusplus
	#define __CL_ENABLE_EXCEPTIONS
	#include"common.hpp"
#endif

// some things, such as unrestricted unions, don't work with gcc4.5
#ifdef __cplusplus
	#if __GNUC__ >= 4 && __GNUC_MINOR__ > 5
		#define GCC46
	#endif
#endif


// AMD does not align unions correctly, help it a bit here
// the 128 is the biggest alignment (on cl_double16), which is very
// wasteful
// optionally, alignment could be specified on each union
// to accomodate the biggest alignment of contained structs
#define AMD_UNION_ALIGN_BUG_WORKAROUND() __attribute__((aligned(128)))


#ifdef __OPENCL_VERSION__
	#define static
	typedef double3 cl_double3;
	typedef short2 cl_short2;
	typedef short cl_short;
	typedef long2 cl_long2;
	typedef bool cl_bool;
	typedef int cl_int;
	typedef ulong cl_ulong;
	typedef long cl_long;
	// printf in OpenCL code
	#ifdef cl_intel_printf
		#pragma OPENCL EXTENSION cl_intel_printf: enable
	#elif defined(cl_amd_printf)
		#pragma OPENCL EXTENSION cl_amd_printf: enable
	#else
		// hope for the best anyway
		// #define printf()
	#endif
	//#ifdef cl_khr_global_int32_base_atomics
		#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
	//#else
	//	#error cl_khr_global_int32_base_atomics extension not supported
	//#endif
#else
	#define global
	#define constant const
	// some funcs for source-level compat
	// NB: there are not real atomics, and are safe only single-threaded!
	template<typename T> T atom_cmpxchg(T* p, const T& cmp, const T& val){ T old=*p; *p=(old==cmp?val:old); return old; }
	template<typename T> T atom_xchg(T* p, const T& val){ T old=*p; *p=val; return old; }
	template<typename T> T atom_inc(T* p){ T old=*p; (*p)++; return old; }
#endif

// these functions belong together
#ifdef __cplusplus

	static boost::tuple<cl::NDRange,cl::NDRange> makeGlobLocNDRanges(size_t globReqSize, bool onePerGroup, cl::Device& dev, cl::Kernel kernel, const string& log=""){
		cl::NDRange globRange, locRange;
		size_t loc[3], glob[3];
		assert(dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>()>=3);
		// https://www.khronos.org/message_boards/viewtopic.php?f=28&t=4423
		size_t maxLocSize=kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(dev);
		assert(dev.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()>=maxLocSize);
		size_t prefGroupSizeMult=kernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(dev);
		vector<size_t> maxItems=dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
		size_t cores=dev.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();

		if(maxItems[0]<prefGroupSizeMult) throw std::logic_error("Maximum items in work-group "+lexical_cast<string>(maxItems[0])+" are smaller than preferred multiple "+lexical_cast<string>(prefGroupSizeMult));

		size_t itemsPerCore=globReqSize/cores;
	
		/*
			1. work-group should have total size a multiple of CL_KERNEL_PREFERRED_WORK_GROUPSIZE_MULTIPLE
			2. work-group dimensions MUST NOT exceed CL_DEVICE_MAX_WORK_ITEM_SIZES

		 */
		if(globReqSize==1 || onePerGroup){ loc[0]=1; loc[1]=1; loc[2]=1; }
		// less items than make all cores busy; use smaller local group, but spread it accross cores better
		else if(itemsPerCore<prefGroupSizeMult){
			loc[0]=itemsPerCore; loc[1]=1; loc[2]=1;
		}
		else {
			// make a the multiple itself, this guarantees b and c can be any numbers
			int a=prefGroupSizeMult;
			int b=min(maxItems[1],maxLocSize/a+(maxLocSize%a?1:0));
			int c=maxLocSize/(a*b)+(maxLocSize%(a*b)?1:0);  // leftovers for c
			loc[0]=a; loc[1]=b; loc[2]=c;
			if((loc[0]*loc[1]*loc[2])%prefGroupSizeMult!=0){ cerr<<"WARN: Local size "<<loc[0]*loc[1]*loc[2]<<" not a multiple of "<<prefGroupSizeMult<<endl; }
		}
		size_t locSize=(loc[0])*(loc[1])*(loc[2]);

		/* check constraints once again */
		if(loc[0]>maxItems[0] || loc[1]>maxItems[1] || loc[2]>maxItems[2]){
			cerr<<"Work item sizes "<<loc[0]<<","<<loc[1]<<","<<loc[2]<<" exceed device maximum "<<maxItems[0]<<","<<maxItems[1]<<","<<maxItems[2]<<endl; 
			throw std::logic_error("CL_DEVICE_MAX_WORK_ITEM_SIZES exceeded. (bug)");
		}
		if(locSize>maxLocSize){ cerr<<"local size "<<loc[0]<<"×"<<loc[1]<<"×"<<loc[2]<<"="<<locSize<<" is greater than maximum work-group size "<<maxLocSize<<endl; throw std::logic_error("CL_KERNEL_WORK_GROUP_SIZE exceeded (bug)"); }
		assert(locSize<=maxLocSize);
		assert(!onePerGroup || locSize==1);

		// glob[i] should be multiples of loc[i]
		glob[0]=((globReqSize%locSize?1:0)+globReqSize/locSize)*loc[0]; glob[1]=loc[1]; glob[2]=loc[2];
		size_t globSize=glob[0]*glob[1]*glob[2];

		if(!log.empty()) cerr<<log<<" "<<globReqSize<<" → "<<glob[0]<<"×"<<glob[1]<<"×"<<glob[2]<<" = ("<<glob[0]/loc[0]<<"×"<<glob[1]/loc[1]<<"×"<<glob[2]/loc[2]<<")×grp("<<loc[0]<<"×"<<loc[1]<<"×"<<loc[2]<<"="<<locSize<<"/"<<maxLocSize<<") = "<<globSize<<" (waste "<<int((globSize-globReqSize)*100./globReqSize)<<"%); "<<(onePerGroup?"[single] ":"")<<"grp max ("<<maxItems[0]<<"×"<<maxItems[1]<<"×"<<maxItems[2]<<")"<<endl;

		if(globSize<globReqSize) throw std::logic_error("Total number of work-items "+lexical_cast<string>(globSize)+" smaller than the required number "+lexical_cast<string>(globReqSize));

		return boost::make_tuple(cl::NDRange(glob[0],glob[1],glob[2]),cl::NDRange(loc[0],loc[1],loc[2]));
	}

	static cl::NDRange makeLinear3DRange(size_t i, const shared_ptr<cl::Device> dev){
		assert(dev->getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>()>=3);
		vector<size_t> wis=dev->getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
		assert(wis.size()>=3);
		wis.resize(3); // in case it is greater than 3
		size_t iMax=wis[0]*wis[1]*wis[2];
		if(i>iMax) throw std::runtime_error("makeLinearNDRange: required number of items ("+lexical_cast<string>(i)+" is bigger than the device's maximum "+lexical_cast<string>(iMax)+" ("+lexical_cast<string>(wis[0])+"×"+lexical_cast<string>(wis[1])+"×"+lexical_cast<string>(wis[2])+")");
		cl::NDRange ret;
		if(i<=wis[2]) ret=cl::NDRange(1,1,i);
		else if(i<=wis[1]*wis[2]) ret=cl::NDRange(1,i/wis[2]+(i%wis[2]==0?0:1),wis[2]);
		else ret = cl::NDRange(i/(wis[1]*wis[2])+(i%(wis[1]*wis[2])==0?0:1),wis[1],wis[2]);
		// const size_t* rr(ret); std::cout<<"NDRange ("<<i<<"≤"<<rr[0]*rr[1]*rr[2]<<"): "<<rr[0]<<","<<rr[1]<<","<<rr[2]<<endl;

		return ret;
	}

#if 0


    static cl::NDRange makeGlobal3DRange(size_t i ,const shared_ptr<cl::Device> dev){
		uint cores = dev->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
		uint tmpCores = trunc(pow(cores, 1.0/3.0)) + 1;
		uint oneCore = i/tmpCores + 1;
		size_t maxItemInGroup = dev->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		uint min = 64;
		if (maxItemInGroup >= 512) {
			if(pow(tmpCores, 3)*min < i) {
				uint multi = trunc(i / tmpCores / (8*8*8)) + 1;		
				uint multi1 = multi*8*tmpCores;
				return cl::NDRange(multi1,multi1,multi1);
			} else {
				return cl::NDRange(4*tmpCores, 4*tmpCores, 4*tmpCores);
			}
		}

		uint multiple = 1;
		if (oneCore > maxItemInGroup) {
			multiple = oneCore / maxItemInGroup + 1;
			oneCore = maxItemInGroup;
		}
		multiple += tmpCores;
		uint oneCoreOneDim = trunc(pow(oneCore/cores, 1.0/3.0)) + 1;
		return cl::NDRange(multiple*oneCoreOneDim, multiple*oneCoreOneDim, multiple*oneCoreOneDim);
	}

    
	static cl::NDRange makeLocal3DRange(size_t i, const shared_ptr<cl::Device> dev){
		uint cores = dev->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
		uint tmpCores = trunc(pow(cores, 1.0/3.0)) + 1;
		uint oneCore = i/tmpCores + 1;
		size_t maxItemInGroup = dev->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		uint min = 64;
		if (maxItemInGroup >= 512){
			if(pow(tmpCores, 3)*min < i) {
				return cl::NDRange(8,8,8);
			} else {
				return cl::NDRange(4,4,4);
			}
		}
	
		if (oneCore > maxItemInGroup) {
			oneCore = maxItemInGroup;
		}
		uint oneCoreOneDim = trunc(pow(oneCore/cores, 1.0/3.0)) + 1;
		return cl::NDRange(oneCoreOneDim, oneCoreOneDim, oneCoreOneDim);
	}
#endif
	
#endif
#ifdef __OPENCL_VERSION__
	static size_t getLinearWorkItem(){ return get_global_id(0)*(get_global_size(1)*get_global_size(2))+get_global_id(1)*get_global_size(2)+get_global_id(2); }
#endif


#define CL_ASSERT2(what,explanation) if(!(what)) printf(__FILE__ ":%d: ERROR: assertion " #what " failed. " explanation "\n",__LINE__)
#define CL_ASSERT(what) CL_ASSERT2(what,)

/* to avoid symbol clashes in shared libs, define our namespace clDem */
#ifdef __cplusplus
	#define CLDEM_NAMESPACE_BEGIN() namespace clDem{
	#define CLDEM_NAMESPACE_END() }; using namespace clDem;
#else
	#define CLDEM_NAMESPACE_BEGIN()
	#define CLDEM_NAMESPACE_END()
#endif

/*
In c++11, unions can contain datatypes with non-trivial ctor, but the union itself
must provide corresponding ctor/copyctor/operator= as well. Math vector types have
initializer_list ctor in c++, hence all unions in which they are must use this macro
to define memcpy-based ctors. In OpenCL code, initializer_list ctor is not defined,
and this macro expands to nothing.
*/
#ifdef __cplusplus
	#include<cstring>
	#define UNION_BITWISE_CTORS(a) a(){}; a(const a& _b){ memcpy(this,&_b,sizeof(a));}; a& operator=(const a& _b){ memcpy(this,&_b,sizeof(a)); return *this; };
#else
	#define UNION_BITWISE_CTORS(a)
#endif

CLDEM_NAMESPACE_BEGIN()

typedef cl_short2 flagSpec; // offset and length, in bits
typedef long par_id_t;

#if 0
	/* strong typedef for cl_long2, so that it can be serialized easily */
	struct par_id2_t{
		typedef cl_long Scalar;
		cl_long2 t;
		cl_long &x, &y, &s0, &s1;
		par_id2_t(): x(t.x), y(t.y), s0(t.x), s1(t.y){}
		par_id2_t(const par_id2_t& other): t(other.t), x(t.x), y(t.y), s0(t.x), s1(t.y){}
		par_id2_t(std::initializer_list<Scalar> l): x(t.x), y(t.y), s0(t.x), s1(t.y){ assert(l.size()==2); int i=0; for(auto n: l) ((Scalar*)this)[i++]=n;  };
		par_id2_t(cl_long xx, cl_long yy): x(t.x), y(t.y), s0(t.x), s1(t.y){ t.x=xx; t.y=yy; }
		par_id2_t& operator=(const par_id2_t& other){ t=other.t; return *this; }
		par_id2_t& operator=(const cl_long2& other){ t=other; return *this; }
		CLDEM_SERIALIZE_ATTRS(/*no regular attrs*/,
			// references don't get saved properly somehow, so work it around here:
			ar&boost::serialization::make_nvp("x",t.x);
			ar&boost::serialization::make_nvp("y",t.y);
		);
	};
#else
	typedef cl_long2 par_id2_t;
#endif

inline int flags_get(const int flags, const flagSpec spec){
	return (flags >> spec.x) & ((1 << spec.y) - 1); 
}

inline void flags_set_global(global int *flags, const flagSpec spec, int val){
	(*flags)&=~(((1<<spec.y)-1)<<spec.x); /* zero field */
	#ifdef __cplusplus
		if(val>=(1<<spec.y)) throw std::runtime_error("Flag field value overflows its size: value "+lexical_cast<string>(val)+", field width "+lexical_cast<string>(spec.y)+", offset "+lexical_cast<string>(spec.x));
	#endif
	val&=((1<<spec.y)-1); /* zero excess bits */
	(*flags)|=val<<spec.x;  /* set field */
}
inline void flags_set_local(int *flags, const flagSpec spec, int val){
	(*flags)&=~(((1<<spec.y)-1)<<spec.x); /* zero field */
	val&=((1<<spec.y)-1); /* zero excess bits */
	(*flags)|=val<<spec.x;  /* set field */
}

CLDEM_NAMESPACE_END()


#endif
