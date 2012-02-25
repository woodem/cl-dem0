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
	#define cl_short2 short2
	#define cl_long2 long2
	#define cl_bool bool
	#define cl_int int
	#define cl_ulong ulong
	// printf in OpenCL code
	#ifdef cl_intel_printf
		#pragma OPENCL EXTENSION cl_intel_printf: enable
	#elif defined(cl_amd_printf)
		#pragma OPENCL EXTENSION cl_amd_printf: enable
	#else
		#define printf(...)
	#endif
#else
	#define global
	#define constant const
#endif

//#ifdef cl_khr_global_int32_base_atomics
	#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
//#else
//	#error cl_khr_global_int32_base_atomics extension not supported
//#endif

#ifdef __cplusplus
	#include<Eigen/Core>
	#include<Eigen/Geometry>
	typedef Eigen::Matrix<Real,3,1> Vector3r;
#endif


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
typedef cl_long2 par_id2_t;

inline int flags_get(const int flags, const flagSpec spec){ return (flags>>spec.x)&((1<<spec.y)-1); }
inline void flags_set(global int *flags, const flagSpec spec, int val){
	(*flags)&=~(((1<<spec.y)-1)<<spec.x); /* zero field */
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
