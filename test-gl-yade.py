# encoding: utf-8
import sys

from yade import *
import yade.log
import yade.gl
from yade import utils, timing
from yade.dem import *
from yade.core import *

from miniEigen import *
from math import *
import pylab, itertools, random

dim=40,40,10
margin=10
r=.005

ktDivKn=.2


m=FrictMat(young=1e6,density=1e3,tanPhi=0)

import yade.pack
sp=yade.pack.SpherePack()
sp.makeCloud((2*r*margin,2*r*margin,2*r),((dim[0]-margin)*2*r,(dim[1]-margin)*2*r,dim[2]*2*r),r,rRelFuzz=.5)
sp.toSimulation(material=m)
for p in O.dem.par: p.vel=(0,0,-.05)

O.dem.par.append([
	utils.wall((0,0,2*r),sense=1,axis=2,material=m),
	utils.wall((0,0,0),sense=1,axis=0,material=m),
	utils.wall((0,0,0),sense=1,axis=1,material=m),
])

O.scene.dt=.5*utils.pWaveDt()

O.scene.engines=[
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()],verletDist=.05*r),
	ContactLoop([Cg2_Sphere_Sphere_L6Geom(),Cg2_Wall_Sphere_L6Geom()],
		[Cp2_FrictMat_FrictPhys(ktDivKn=ktDivKn)],
		[Law2_L6Geom_FrictPhys_IdealElPl(noSlip=True)],
		applyForces=True
	),
	Gravity(gravity=(-4,-5,-10)),
	Leapfrog(damping=.4,reset=True),
]
O.timingEnabled=True
O.scene.trackEnergy=True
O.saveTmp()

if 1:
	#pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
	#dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

	O.scene.clDev=1,0


	#sys.path.append('.')
	#import clDem
	#sim=clDem.Simulation(1,0,breakTension=True,ktDivKn=ktDivKn)

	if 0:
		sim.scene.materials=[clDem.ElastMat(young=1e6,density=1e3)]
		sim.scene.gravity=(-4,-5,-10)
		sim.scene.damping=.4
		sim.scene.verletDist=.3*r # collision detection in this case
		sim.maxScheduledSteps=10

		sim.par.append(clDem.mkWall(pos=(0,0,2*r),axis=2,sim=sim,matId=0,groups=0b011))
		sim.par.append(clDem.mkWall(pos=(0,0,0),axis=0,sim=sim,matId=0,groups=0b011))
		sim.par.append(clDem.mkWall(pos=(0,0,0),axis=1,sim=sim,matId=0,groups=0b011))
		# walls are in loneGroups, but collide with spheres below (0b001) as well
		sim.scene.loneGroups=0b010

		# ground
		#for x0,x1 in itertools.product(range(0,dim[0]),range(0,dim[1])):
		#	sim.par.append(clDem.mkSphere((x0*2*r,x1*2*r,0),r,sim,matId=0,groups=0b011,fixed=True))

		import yade.pack
		sp=yade.pack.SpherePack()
		sp.makeCloud((2*r*margin,2*r*margin,2*r),((dim[0]-margin)*2*r,(dim[1]-margin)*2*r,dim[2]*2*r),r,rRelFuzz=.5)
		for center,radius in sp:
			sim.par.append(clDem.mkSphere(center,radius,sim,matId=0,groups=0b001,fixed=False))
			sim.par[-1].vel=(0,0,-.05)



