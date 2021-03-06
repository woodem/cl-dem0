# encoding: utf-8
# import the c++ wrapper
import sys
sys.path.append('.')
import _clDem
import math
import minieigen

import functools
# create new instance of baseClass and set (existing) attributes from keywords
def BaseClassSetter(baseClass,**kw):
	s=baseClass()
	for key,val in kw.items():
		if hasattr(s,key): setattr(s,key,val)
		else: raise AttributeError("No such attribute: %s"%key)
	return s
# functions which look like constructors, but set all attributes from keyword passed to them
# the actual classes are in the _clDem module, which is not imported into the namespace
# of this module directly
for clss in ('Scene','ElastMat','FrictMat','Particle','Sphere','Wall','Contact','L1Geom','L6Geom','NormPhys','FrictPhys'):
	# global() is reference to the global (module, in this case) dictionary
	globals()[clss]=functools.partial(BaseClassSetter,getattr(_clDem,clss))
Simulation=_clDem.Simulation

# utility functions here
def mkSphere(pos,radius,sim,matId=0,groups=None,fixed=False):
	rho=sim.scene.materials[matId].density
	mass=rho*(4/3.)*math.pi*radius**3
	inert=mass*(2/5.)*radius**2
	p=Particle(pos=pos,mass=mass,inertia=(inert,inert,inert),dofs=0 if fixed else 63,shape=Sphere(radius=radius),matId=matId)
	if groups!=None: p.groups=groups
	return p

def mkWall(pos,axis,sim,matId=0,sense=0,groups=None,fixed=True):
	p=Particle(pos=pos,mass=-1,inertia=(0,0,0),dofs=0 if fixed else 63,shape=Wall(axis=axis,sense=sense),matId=matId)
	if groups!=None: p.groups=groups
	return p

def briefOutput():
	# output similar to the c++ code
	minieigen.Vector3.__str__=lambda s: '('+','.join('%g'%e for e in s)+')'
	minieigen.Matrix3.__str__=lambda s: '('+', '.join(','.join('%g'%e for e in s.row(i)) for i in (0,1,2))+')'

def showSim(sim):
	print 'At step %d (t=%g), Δt=%g'%(sim.scene.step,sim.scene.t,sim.scene.dt)
	for i,p in enumerate(sim.par):
		print '#%d x=%s; v=%s, F=%s, T=%s'%(i,p.pos,p.vel,p.force,p.torque)
	for c in sim.con:
		print '##%d+%d: p=%s, F=%s, T=%s, rot=%s'%(c.ids[0],c.ids[1],c.pos,c.force,c.torque,c.ori)
		if c.geom: print '\tgeomT=%d, uN=%g%s'%(c.geomT,c.geom.uN,', v=%s, ω=%s'%(c.geom.vel,c.geom.angVel) if isinstance(c.geom,_clDem.L6Geom) else '')
		if c.phys: print '\tphysT=%d, kN=%g'%(c.physT,c.phys.kN)
Simulation.show=showSim


