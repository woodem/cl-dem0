#ifndef _COMMON_CL_H
#define _COMMON_CL_H

#include"../../cl-math/basic-math.cl"
#ifdef __cplusplus
	#define __CL_ENABLE_EXCEPTIONS
	#include"common.hpp"
#endif

// AMD does not align unions correctly, help it a bit here
// the 128 is the biggest alignment (on cl_double16), which is very
// wasteful
// optionally, alignment could be specified on each union
// to accomodate the biggest alignment of contained structs
#define AMD_UNION_ALIGN_BUG_WORKAROUND() __attribute__((aligned(128)))


#ifdef __OPENCL_VERSION__
	#define static
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
		#define printf(...)
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


    static cl::NDRange makeGlobal3DRange(size_t i ,const shared_ptr<cl::Device> dev){
		uint cores = dev->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
		int oneCore = i / cores;
		size_t maxItemInGroup = dev->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
	
		uint multiple = 1;
		if (oneCore > maxItemInGroup) {
			multiple = oneCore / maxItemInGroup;
			oneCore = maxItemInGroup;
		}
		uint oneCoreOneDim = trunc(pow(oneCore, 1.0/3.0)) + 1;

	    return cl::NDRange(multiple * oneCoreOneDim * 2, multiple * oneCoreOneDim * 2,
				multiple * oneCoreOneDim * 2);
	}

    
	static cl::NDRange makeLocal3DRange(size_t i, const shared_ptr<cl::Device> dev){
		uint cores = dev->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
		int oneCore = i / cores;
		size_t maxItemInGroup = dev->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
	
		uint multiple = 1;
		if (oneCore > maxItemInGroup) {
			multiple = oneCore / maxItemInGroup;
			oneCore = maxItemInGroup;
		}
		uint oneCoreOneDim = trunc(pow(oneCore, 1.0/3.0)) + 1;

		return cl::NDRange(oneCoreOneDim, oneCoreOneDim, oneCoreOneDim);
	}

	
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

inline int flags_get(const int flags, const flagSpec spec){ return (flags>>spec.x)&((1<<spec.y)-1); }
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
