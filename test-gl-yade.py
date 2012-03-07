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
		[Law2_L6Geom_FrictPhys_LinEl6(noSlip=True)]
	),
	IntraForce([In2_Sphere_ElastMat()]),
	Gravity(gravity=(-4,-5,-10)),
	Leapfrog(damping=.4,reset=True),
]
O.timingEnabled=True
O.scene.trackEnergy=True
O.saveTmp()



