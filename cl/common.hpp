#pragma once

#define __CL_ENABLE_EXCEPTIONS
#include"cl.hpp"

#include<boost/python.hpp>
namespace py=boost::python;

#define PY_RWV(clss,attr) add_property(BOOST_PP_STRINGIZE(attr),/*read access*/py::make_getter(&clss::attr,py::return_value_policy<py::return_by_value>()),/*write access*/make_setter(&clss::attr,py::return_value_policy<py::return_by_value>()))
#define PY_RW(clss,attr) def_readwrite(BOOST_PP_STRINGIZE(attr),&clss::attr)

// workaround for Contact alignment issues
// http://thread.gmane.org/gmane.comp.python.c++/15639
namespace boost {
	namespace align { struct __attribute__((__aligned__(128))) a128 {};}
	template<> class type_with_alignment<128> { public: typedef align::a128 type; };
	// namespace detail{ BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(is_pod,::boost::align::a128,true) }
};

/* self-stolen from Yade */

/*** c++-list to python-list */
template<typename containedType>
struct custom_vector_to_list{
	static PyObject* convert(const std::vector<containedType>& v){
		py::list ret; for(const containedType& e: v) ret.append(e);
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

