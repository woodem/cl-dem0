# import the c++ wrapper
import sys
sys.path.append('.')
import _miniDem

import functools
# create new instance of baseClass and set (existing) attributes from keywords
def BaseClassSetter(baseClass,**kw):
	s=baseClass()
	for key,val in kw.items():
		if hasattr(s,key): setattr(s,key,val)
		else: raise AttributeError("No such attribute: %s"%key)
	return s
# functions which look like constructors, but set all attributes from keyword passed to them
# the actual classes are in the _miniDem module, which is not imported into the namespace
# of this module directly
for clss in ('Scene','Simulation','ElastMat','Particle','Sphere','Contact','L1Geom','NormPhys'):
	# global() is reference to the global (module, in this case) dictionary
	globals()[clss]=functools.partial(BaseClassSetter,getattr(_miniDem,clss))
