#include"Simulation.hpp"

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Geometry>
/* miniEigen types */
typedef double Real;
typedef Eigen::Matrix<int,2,1> Vector2i;
typedef Eigen::Matrix<Real,3,1> Vector3r;
typedef Eigen::Matrix<Real,3,3> Matrix3r;
typedef Eigen::Quaternion<Real> Quaternionr;




template<typename T> struct get_scalar{ typedef typename T::Scalar type; };
template<> struct get_scalar<cl_long2>{ typedef cl_long type; };

/* fixed-length vector types */
template<typename clType, typename eigType, int len>
struct custom_clType_to_eigType{
	typedef typename get_scalar<clType>::type clType_Scalar;
	static PyObject* convert(const clType& a){
		eigType ret;
		for(int i=0; i<len; i++) ((typename eigType::Scalar*)&ret)[i]=((clType_Scalar*)&a)[i];
		return py::incref(py::object(ret).ptr());
	}
};

template<typename clType, typename eigType, int len>
struct custom_clType_from_eigType{
	typedef typename get_scalar<clType>::type Scalar;
	custom_clType_from_eigType(){ py::converter::registry::push_back(&convertible,&construct,py::type_id<clType>()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		if(PySequence_Length(obj_ptr)!=len) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		 void* storage=((py::converter::rvalue_from_python_storage<clType>*)(data))->storage.bytes;
		 new (storage) clType();
		 clType* v=(clType*)(storage);
		 int l=PySequence_Size(obj_ptr); assert(l>0);
		 for(int i=0; i<l; i++) { (((Scalar*)v)[i]=py::extract<Scalar>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};




BOOST_PYTHON_MODULE(_clDem){
	/*
	NOTE:
	* types which are wrapped directly (such as scalars) can be returned by reference, with PY_RW
	* types which are not directly wrapped (such as Vec3, which is converted to and from eigen's Vector3r etc) must be returned by value, with PY_RWV
	*/

	py::to_python_converter<cl_long2,custom_clType_to_eigType<cl_long2,Vector2i,2>>();
	py::to_python_converter<Vec3    ,custom_clType_to_eigType<Vec3,Vector3r,3>>();
	py::to_python_converter<Quat    ,custom_clType_to_eigType<Quat,Quaternionr,4>>();
	py::to_python_converter<Mat3    ,custom_clType_to_eigType<Mat3,Matrix3r,9>>();

	custom_clType_from_eigType<cl_long2,Vector2i,2>();
	custom_clType_from_eigType<Vec3,Vector3r,3>();
	custom_clType_from_eigType<Mat3,Matrix3r,9>();
	custom_clType_from_eigType<Quat,Quaternionr,4>();

	VECTOR_SEQ_CONV(cl_long2);
	VECTOR_SEQ_CONV(int);

	py::class_<std::vector<cl_float>>("FloatList").def(py::vector_indexing_suite<std::vector<cl_float>>());
	custom_vector_from_seq<cl_float>();
	py::class_<std::vector<cl_double>>("DoubleList").def(py::vector_indexing_suite<std::vector<cl_double>>());
	custom_vector_from_seq<cl_double>();

	Particle_cl_h_expose();
	Contact_cl_h_expose();
	Scene_cl_h_expose();
	Simulation_hpp_expose();
}