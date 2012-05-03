#pragma once

#define __CL_ENABLE_EXCEPTIONS
#include"cl.hpp"

#include"serialization.cl.h"

#include<boost/shared_ptr.hpp>
#include<boost/make_shared.hpp>
using boost::shared_ptr;
using boost::make_shared;
#include<boost/tuple/tuple.hpp>

#include<boost/python.hpp>
namespace py=boost::python;

#include<boost/foreach.hpp>
#define FOREACH BOOST_FOREACH

#define PY_RWV(clss,attr) add_property(BOOST_PP_STRINGIZE(attr),/*read access*/py::make_getter(&clss::attr,py::return_value_policy<py::return_by_value>()),/*write access*/make_setter(&clss::attr,py::return_value_policy<py::return_by_value>()))
#define PY_RW(clss,attr) def_readwrite(BOOST_PP_STRINGIZE(attr),&clss::attr)

// when included from yade, this is already in the yade's header
#if BOOST_VERSION<=104900 && !defined(YADE_CLDEM)
	// workaround for Contact alignment issues http://thread.gmane.org/gmane.comp.python.c++/15639
	// it was fixed in the SVN http://svn.boost.org/svn/boost/trunk/boost/type_traits/type_with_alignment.hpp
	// before boost 1.50
	namespace boost {
		namespace align { struct __attribute__((__aligned__(128))) a128 {};}
		template<> class type_with_alignment<128> { public: typedef align::a128 type; };
	};
#endif

/* self-stolen from Yade */

/*** c++-list to python-list */
template<typename containedType>
struct custom_vector_to_list{
	static PyObject* convert(const std::vector<containedType>& v){
		py::list ret; /*for(const containedType& e: v)*/ FOREACH(const containedType& e, v) ret.append(e);
		return py::incref(ret.ptr());
	}
};;

template<typename containedType>
struct custom_vector_from_seq{
	custom_vector_from_seq(){ py::converter::registry::push_back(&convertible,&construct,py::type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		 void* storage=((py::converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(py::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};

#define VECTOR_SEQ_CONV(Type) custom_vector_from_seq<Type>();  py::to_python_converter<std::vector<Type>, custom_vector_to_list<Type> >();



#include<boost/lexical_cast.hpp>
#include<boost/static_assert.hpp>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cassert>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>

using std::vector;
using boost::lexical_cast;
using std::string;
using std::cerr;
using std::endl;
using std::max;
using std::min;


/* eigen stuff */
#define EIGEN_DONT_ALIGN
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Eigenvalues>

namespace clDem{
	typedef Eigen::Matrix<int,2,1> Vector2i;
	typedef Eigen::Matrix<Real,3,3> Matrix3r;
	typedef Eigen::Matrix<Real,3,1> Vector3r;
	typedef Eigen::Quaternion<Real> Quaternionr;

	inline Vector3r toEigen(const Vec3& v){ return Vector3r(v[0],v[1],v[2]); }
	inline Quaternionr toEigen(const Quat& q){ return Quaternionr(q[3],q[0],q[1],q[2]); }
	inline Matrix3r toEigen(const Mat3& m){ Matrix3r ret; ret<<m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8]; return ret; }

	inline Vec3 fromEigen(const Vector3r& v){ return Vec3_set(v.x(),v.y(),v.z()); }
	inline Quat fromEigen(const Quaternionr& q){ return Quat_set_wxyz(q.w(),q.x(),q.y(),q.z()); }
	inline Mat3 fromEigen(const Matrix3r& m){ return Mat3_set(m(0,0),m(0,1),m(0,2),m(1,0),m(1,1),m(1,2),m(2,0),m(2,1),m(2,2)); }
};
