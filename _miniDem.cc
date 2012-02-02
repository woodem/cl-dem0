/* :vim:makeprg=make: */
#include<boost/python.hpp>
namespace py=boost::python;

// workaround for Contact alignment issues
// have to check on a mailing list whether it is OK
namespace boost {
	namespace align { struct __attribute__((__aligned__(128))) a128 {};}
	template<> class type_with_alignment<128> { public: typedef align::a128 type; };
	// namespace detail{ BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(is_pod,::boost::align::a128,true) }
};

#include<boost/lexical_cast.hpp>
#include<boost/static_assert.hpp>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cassert>

using std::vector;
using boost::lexical_cast;
using std::string;


#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Geometry>
/* miniEigen types */
typedef double Real;
typedef Eigen::Matrix<int,2,1> Vector2i;
typedef Eigen::Matrix<Real,3,1> Vector3r;
typedef Eigen::Matrix<Real,3,3> Matrix3r;
typedef Eigen::Quaternion<Real> Quaternionr;


#define __CL_ENABLE_EXCEPTIONS
#include"cl.hpp"



#include"scene.cl"

py::object Material_mat_get(const Material* m){
	int matT=mat_matT_get(m);
	switch(matT){
		case 0: return py::object();
		case Mat_ElastMat: return py::object(m->mat.elast);
		default: throw std::runtime_error("Material has mat with unknown index "+lexical_cast<string>(matT));
	}
}

void Material_mat_set(Material *m, py::object mat){
	py::extract<ElastMat> elast(mat);
	if(elast.check()){ m->mat.elast=elast(); mat_matT_set(m,Mat_ElastMat); }
	else throw std::runtime_error("Unknown mat object.");
}

py::list Scene_mats_get(Scene* self){
	py::list ret;
	for(const Material& m: self->materials) ret.append(Material_mat_get(&m));
	return ret;
}

void Scene_mats_set(Scene* self, const py::list mm){
	if(py::len(mm)>=SCENE_NUM_MAT) throw std::runtime_error("Too many materials, to be defined (maximum "+lexical_cast<string>(SCENE_NUM_MAT)+", given "+lexical_cast<string>(py::len(mm)));
	for(int i=0; i<py::len(mm); i++){ Material_mat_set(&(self->materials[i]),mm[i]); }
};

py::object Particle_shape_get(Particle* p){
	int shapeT=par_shapeT_get(p);
	switch(shapeT){
		case 0: return py::object(); // None
		case Shape_Sphere: return py::object(p->shape.sphere);
		default:	throw std::runtime_error("Particle has shape with unknown index "+lexical_cast<string>(shapeT));
	}
}
void Particle_shape_set(Particle *p, py::object sh){
	py::extract<Sphere> sphere(sh);
	if(sphere.check()){
		p->shape.sphere=sphere();
		par_shapeT_set(p,Shape_Sphere);
	}
	else throw std::runtime_error("Unknown shape object.");
}

py::object Contact_geom_get(Contact* c){
	int geomT=con_geomT_get(c);
	switch(geomT){
		case 0: return py::object();
		case Geom_L1Geom: return py::object(c->geom.l1g);
		default: throw std::runtime_error("Contact has geom with unknown index "+lexical_cast<string>(geomT));
	}
}

void Contact_geom_set(Contact *c,py::object g){
	py::extract<L1Geom> l1g(g);
	if(l1g.check()){ c->geom.l1g=l1g(); con_geomT_set(c,Geom_L1Geom); }
	else throw std::runtime_error("Unknown geom object.");
}

py::object Contact_phys_get(Contact* c){
	int physT=con_physT_get(c);
	switch(physT){
		case 0: return py::object();
		case Phys_NormPhys: return py::object(c->phys.normPhys);
		default: throw std::runtime_error("Contact has phys with unknown index "+lexical_cast<string>(physT));
	}
}

void Contact_phys_set(Contact *c,py::object p){
	py::extract<NormPhys> np(p);
	if(np.check()){ c->phys.normPhys=np(); con_geomT_set(c,Phys_NormPhys); }
	else throw std::runtime_error("Unknown geom object.");
}

struct Simulation{
	Scene scene;
	vector<Particle> par;
	vector<Contact> con;

	cl::Platform platform;
	cl::Device device;
	cl::Context context;
	cl::CommandQueue queue;
	cl::Program program;

	Simulation(int pNum=-1,int dNum=-1){
		scene=Scene_new();
		initCl(pNum,dNum);
	}
	void initCl(int pNum, int dNum){
		std::vector<cl::Platform> platforms;
		std::vector<cl::Device> devices;
		cl::Platform::get(&platforms);
		if(pNum<0){
			std::cerr<<"==================================="<<std::endl;
			for(int i=0; i<platforms.size(); i++){
				std::cerr<<i<<". platform: "<<platforms[i].getInfo<CL_PLATFORM_NAME>()<<std::endl;
				platforms[i].getDevices(CL_DEVICE_TYPE_ALL,&devices);
				for(int j=0; j<devices.size(); j++){
					std::cerr<<"\t"<<j<<". device: "<<devices[j].getInfo<CL_DEVICE_NAME>()<<std::endl;
				}
			}
			std::cerr<<"==================================="<<std::endl;
		}
		cl::Platform::get(&platforms);
		if(pNum>=(int)platforms.size()){ std::cerr<<"Only "<<platforms.size()<<" platforms available, taking 0th platform."<<std::endl; pNum=0; }
		if(pNum<0) pNum=0;
		platform=platforms[pNum];
		platforms[pNum].getDevices(CL_DEVICE_TYPE_ALL,&devices);
		if(dNum>=(int)devices.size()){ std::cerr<<"Only "<<devices.size()<<" devices available, taking 0th platform."<<std::endl; dNum=0; }
		if(dNum<0) dNum=0;
		context=cl::Context(devices);
		device=devices[dNum];
		std::cerr<<"** OpenCL ready: platform \""<<platform.getInfo<CL_PLATFORM_NAME>()<<"\", device \""<<device.getInfo<CL_DEVICE_NAME>()<<"\"."<<std::endl;
		queue=cl::CommandQueue(context,device);
		// compile source
		const char* src="#include\"scene.cl\"\n";
		cl::Program::Sources source(1,std::make_pair(src,std::string(src).size()));
		program=cl::Program(context,source);
		try{
			program.build(std::vector<cl::Device>({device}),"-I.",NULL,NULL);
		}catch(cl::Error& e){
			std::cerr<<"Error building source. Build log:\n"<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)<<std::endl;
			throw std::runtime_error("Error compiling OpenCL code.");
		}
		auto log=program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		if(!log.empty()) std::cerr<<log<<std::endl;
		std::cerr<<"** Program compiled.\n";
	};
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


BOOST_PYTHON_MODULE(_miniDem){
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

	#define PY_RWV(clss,attr) add_property(BOOST_PP_STRINGIZE(attr),/*read access*/py::make_getter(&clss::attr,py::return_value_policy<py::return_by_value>()),/*write access*/make_setter(&clss::attr,py::return_value_policy<py::return_by_value>()))
	#define PY_RW(clss,attr) def_readwrite(BOOST_PP_STRINGIZE(attr),&clss::attr)


	py::class_<Simulation>("Simulation",py::init<int,int>((py::arg("platformNum")=-1,py::arg("deviceNum")=-1)))
		.PY_RW(Simulation,scene)
		.PY_RWV(Simulation,par)
		.PY_RWV(Simulation,con)
	;
	py::class_<Scene>("Scene").def("__init__",Scene_new)
		.PY_RW(Scene,t).PY_RW(Scene,dt).PY_RWV(Scene,gravity).PY_RWV(Scene,damping)
		//.PY_RW_BYVALUE(Scene,materials)
		.add_property("materials",Scene_mats_get,Scene_mats_set)
	;
	#if 0
		py::class_<Material>("Material").add_property("mat",Material_mat_get,Material_mat_set);
	#endif

	py::class_<Particle>("Particle").def("__init__",Particle_new)
		.PY_RWV(Particle,pos).PY_RWV(Particle,ori).PY_RWV(Particle,inertia).PY_RW(Particle,mass).PY_RWV(Particle,vel).PY_RWV(Particle,force).PY_RWV(Particle,torque)
		.add_property("shape",Particle_shape_get,Particle_shape_set)
		// flags
		.def_readonly("flags",&Particle::flags)
		.add_property("shapeT",par_shapeT_get)
		.add_property("clumped",par_clumped_get,par_clumped_set)
		.add_property("stateT",par_stateT_get)
		.add_property("dofs",par_dofs_get,par_dofs_set)
		.add_property("groups",par_groups_get,par_groups_set)
		.add_property("matId",par_matId_get,par_matId_set)
	;
	py::class_<Contact>("Contact").def("__init__",Contact_new)
		.PY_RWV(Contact,ids).PY_RWV(Contact,pos).PY_RWV(Contact,ori).PY_RWV(Contact,force).PY_RWV(Contact,torque)
		.def_readonly("flags",&Contact::flags)
		.add_property("geom",Contact_geom_get,Contact_geom_set)
		.add_property("phys",Contact_phys_get,Contact_phys_set)
		.add_property("shapesT",con_shapesT_get)
		.add_property("geomT",con_geomT_get)
		.add_property("physT",con_physT_get)
	;
	// "derived" classes (union members)
	py::class_<Sphere>("Sphere").def("__init__",Sphere_new).PY_RW(Sphere,radius);
	py::class_<ElastMat>("ElastMat").def("__init__",ElastMat_new).PY_RW(ElastMat,young);
	py::class_<L1Geom>("L1Geom").def("__init__",L1Geom_new).PY_RW(L1Geom,uN);
	py::class_<NormPhys>("NormPhys").def("__init__",NormPhys_new).PY_RW(NormPhys,kN);

	#define VECTOR_SEQ_CONV(Type) custom_vector_from_seq<Type>();  py::to_python_converter<std::vector<Type>, custom_vector_to_list<Type> >();
	VECTOR_SEQ_CONV(Particle);
	VECTOR_SEQ_CONV(Contact);
	VECTOR_SEQ_CONV(Material);
};
