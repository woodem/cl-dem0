#ifndef _CLDEM_SERIALIZATION_HPP_
#define _CLDEM_SERIALIZATION_HPP_

#ifndef __cplusplus
	// for OpenCL
	#define CLDEM_SERIALIZE_ATTRS(attrs,otherCode)
#else
	// include archive types
	#include<boost/archive/binary_oarchive.hpp>
	#include<boost/archive/binary_iarchive.hpp>
	#include<boost/archive/xml_oarchive.hpp>
	#include<boost/archive/xml_iarchive.hpp>

	#include<boost/serialization/vector.hpp>
	#include<boost/serialization/shared_ptr.hpp>
	#include<boost/serialization/list.hpp>
	#include<boost/serialization/map.hpp>
	// what we need below
	#include<boost/serialization/nvp.hpp>
	#include<boost/serialization/is_bitwise_serializable.hpp>

	// perhaps not needed?
	//#include<boost/serialization/export.hpp>

	#include<boost/preprocessor.hpp>

	// macros for attribute serialization
	#define _CLDEM_SER_ATTR(x,y,z) /* std::cerr<<"["<<BOOST_PP_STRINGIZE(z)<<"]"; */ ar & BOOST_SERIALIZATION_NVP(z);
	#define CLDEM_SERIALIZE_ATTRS(attrs,otherCode) \
		friend class boost::serialization::access; \
		template<class ArchiveT> void serialize(ArchiveT & ar, unsigned int version){ \
			BOOST_PP_SEQ_FOR_EACH(_CLDEM_SER_ATTR,~,attrs) \
			otherCode \
		;}


	#include<CL/cl.h>
	/* serialization code for (some) OpenCL vector types */
	//BOOST_IS_BITWISE_SERIALIZABLE(cl_int2);
	BOOST_IS_BITWISE_SERIALIZABLE(cl_long2);
	BOOST_CLASS_IMPLEMENTATION(cl_long2,boost::serialization::object_serializable);

	BOOST_IS_BITWISE_SERIALIZABLE(cl_double4);
	BOOST_CLASS_IMPLEMENTATION(cl_double4,boost::serialization::object_serializable);

	namespace boost{ namespace serialization {

	template<class Archive> void serialize(Archive &ar, cl_long2 &i, const unsigned version){
		ar & make_nvp("x",i.x);
		ar & make_nvp("y",i.y);
	}
	/* hack for gcc 4.5 */
	template<class Archive> void serialize(Archive &ar, cl_double4& v, const unsigned version){ ar & boost::serialization::make_nvp("x",v.x); ar & boost::serialization::make_nvp("y",v.y); ar & boost::serialization::make_nvp("z",v.z); ar & boost::serialization::make_nvp("w",v.w); }
	/* end hack */


	}};

#endif

#endif
