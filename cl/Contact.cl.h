#ifndef _CONTACT_CL_H_
#define _CONTACT_CL_H_

#include"common.cl.h"

CLDEM_NAMESPACE_BEGIN();

struct L1Geom{ Real uN; };
static struct L1Geom L1Geom_new(){ struct L1Geom ret; ret.uN=NAN; return ret; }

struct L6Geom{ Real uN; Vec3 vel; Vec3 angVel; };
static struct L6Geom L6Geom_new(){ struct L6Geom ret; ret.uN=NAN; ret.vel=ret.angVel=Vec3_set(0.,0.,0.); return ret; }

enum _geom_enum { Geom_L1Geom=1, Geom_L6Geom };

struct NormPhys{ Real kN; };
static struct NormPhys NormPhys_new(){ struct NormPhys ret; ret.kN=NAN; return ret; }

enum _phys_enum { Phys_NormPhys=1, };

struct Contact{
	int flags;
	par_id2_t ids;
	Vec3 pos;
	Mat3 ori;
	Vec3 force, torque;
	union _gg {
		UNION_BITWISE_CTORS(_gg)
		struct L1Geom l1g;
		struct L6Geom l6g;
	}  AMD_UNION_ALIGN_BUG_WORKAROUND() geom;
	union {
		struct NormPhys normPhys;
	}  AMD_UNION_ALIGN_BUG_WORKAROUND() phys;
};

#define CON_LEN_shapesT 2*(PAR_LEN_shapeT) // currently unused in the actual code
#define CON_LEN_geomT   2
#define CON_LEN_physT   2

#define CON_OFF_shapesT 0
#define CON_OFF_geomT   CON_OFF_shapesT + CON_LEN_shapesT
#define CON_OFF_physT   CON_OFF_geomT + CON_LEN_geomT

constant flagSpec con_flag_shapesT={CON_OFF_shapesT,CON_LEN_shapesT};
constant flagSpec con_flag_geomT  ={CON_OFF_geomT,CON_LEN_geomT};
constant flagSpec con_flag_physT  ={CON_OFF_physT,CON_LEN_physT};
#define CONTACT_FLAG_GET_SET(what) \
	inline int con_##what##_get(global const struct Contact *c){ return flags_get(c->flags,con_flag_##what); } \
	inline void con_##what##_set(global struct Contact *c, int val){ flags_set(&c->flags,con_flag_##what,val); } \
	inline int con_##what##_get_local(const struct Contact *c){ return flags_get(c->flags,con_flag_##what); } \
	inline void con_##what##_set_local(struct Contact *c, int val){ flags_set_local(&c->flags,con_flag_##what,val); } 
CONTACT_FLAG_GET_SET(shapesT);
CONTACT_FLAG_GET_SET(geomT);
CONTACT_FLAG_GET_SET(physT);

static struct Contact Contact_new(){
	struct Contact c;
	c.ids.s0=-1; c.ids.s1=-1;
	c.pos=Vec3_set(0,0,0);
	c.ori=Mat3_identity();
	c.force=c.torque=Vec3_set(0,0,0);
	con_shapesT_set_local(&c,0);
	con_geomT_set_local(&c,0);
	con_physT_set_local(&c,0);
	return c;
}

// must be macros so that they can be used as switch cases
#define SHAPET2_COMBINE(s1,s2) ((s1) | (s2)<<(PAR_LEN_shapeT))
#define GEOMT_PHYST_COMBINE(g,p) ((g) | (p)<<(CON_LEN_geomT))


CLDEM_NAMESPACE_END();

#ifdef __cplusplus
	// needed for py::indexing_suite
	namespace clDem{
		inline bool operator==(const Contact& a, const Contact& b){ return memcmp(&a,&b,sizeof(Contact))==0; }	
	};

	static
	py::object Contact_geom_get(Contact* c){
		int geomT=con_geomT_get(c);
		switch(geomT){
			case 0: return py::object();
			case Geom_L1Geom: return py::object(c->geom.l1g);
			case Geom_L6Geom: return py::object(c->geom.l6g);
			default: throw std::runtime_error("Contact has geom with unknown index "+lexical_cast<string>(geomT));
		}
	}

	static
	void Contact_geom_set(Contact *c,py::object g){
		if(g==py::object()){ con_geomT_set(c,0); return; }
		py::extract<L1Geom> l1g(g);
		if(l1g.check()){ c->geom.l1g=l1g(); con_geomT_set(c,Geom_L1Geom); }
		else throw std::runtime_error("Unknown geom object.");
	}

	static
	py::object Contact_phys_get(Contact* c){
		int physT=con_physT_get(c);
		switch(physT){
			case 0: return py::object();
			case Phys_NormPhys: return py::object(c->phys.normPhys);
			default: throw std::runtime_error("Contact has phys with unknown index "+lexical_cast<string>(physT));
		}
	}

	static
	void Contact_phys_set(Contact *c,py::object p){
		if(p==py::object()){ con_physT_set(c,0); return; }
		py::extract<NormPhys> np(p);
		if(np.check()){ c->phys.normPhys=np(); con_geomT_set(c,Phys_NormPhys); }
		else throw std::runtime_error("Unknown geom object.");
	}
	template<int N>
	par_id_t Contact_id_get(Contact *c){ return (N==0?c->ids.s0:c->ids.s1); }
	
	static
	void Contact_cl_h_expose(){
		py::class_<L1Geom>("L1Geom").def("__init__",L1Geom_new).PY_RW(L1Geom,uN);
		py::class_<L6Geom>("L6Geom").def("__init__",L6Geom_new).PY_RWV(L6Geom,uN).PY_RWV(L6Geom,vel).PY_RWV(L6Geom,angVel);
		py::class_<NormPhys>("NormPhys").def("__init__",NormPhys_new).PY_RW(NormPhys,kN);

		py::class_<Contact>("Contact").def("__init__",Contact_new)
			.PY_RWV(Contact,ids).PY_RWV(Contact,pos).PY_RWV(Contact,ori).PY_RWV(Contact,force).PY_RWV(Contact,torque)
			.def_readonly("flags",&Contact::flags)
			.add_property("geom",Contact_geom_get,Contact_geom_set)
			.add_property("phys",Contact_phys_get,Contact_phys_set)
			.add_property("shapesT",con_shapesT_get)
			.add_property("geomT",con_geomT_get)
			.add_property("physT",con_physT_get)
			.add_property("id1",Contact_id_get<0>)
			.add_property("id2",Contact_id_get<1>)
		;
		// to-python from [...]
		custom_vector_from_seq<Contact>();
		// from-python as ContactList proxy
		py::class_<std::vector<Contact>>("ContactList").def(py::vector_indexing_suite<std::vector<Contact>>());
	};
#endif




#endif
