# encoding: utf-8
import sys
sys.path.append('.')
import clDem

from miniEigen import *
from math import *
import pylab, itertools, random

dim=30,30,4
margin=12
r=.01
ktDivKn=.2

sim=clDem.Simulation(breakTension=True,ktDivKn=ktDivKn)

sim.scene.materials=[clDem.FrictMat(young=1e6,density=1e3,ktDivKn=.2,tanPhi=.5)]
sim.scene.gravity=(-4,-5,-10)
sim.scene.damping=.4
sim.scene.verletDist=.3*r
sim.maxScheduledSteps=10

#sim.par.append(clDem.mkWall(pos=(0,0,2*r),axis=2,sim=sim,matId=0,groups=0b011))
#sim.par.append(clDem.mkWall(pos=(0,0,0),axis=0,sim=sim,matId=0,groups=0b011))
#sim.par.append(clDem.mkWall(pos=(0,0,0),axis=1,sim=sim,matId=0,groups=0b011))

# ground
for x0,x1 in itertools.product(range(0,dim[0]),range(0,dim[1])):
	sim.par.append(clDem.mkSphere((x0*2*r,x1*2*r,0),r,sim,matId=0,groups=0b011,fixed=True))

# walls/ground spheres are in loneGroups, but collide with spheres below (0b001) as well
sim.scene.loneGroups=0b010

import yade.pack
sp=yade.pack.SpherePack()
sp.makeCloud((2*r*margin,2*r*margin,2*r),((dim[0]-margin)*2*r,(dim[1]-margin)*2*r,dim[2]*2*r),r,rRelFuzz=.5)
for center,radius in sp:
	sim.par.append(clDem.mkSphere(center,radius,sim,matId=0,groups=0b001,fixed=False))
	sim.par[-1].vel=(0,0,-.05)

O.scene.dt=sim.scene.dt=.5*sim.pWaveDt()


from yade import *
import yade.cld
import yade.log
import yade.gl
from yade import timing
if 1: # run on both
	O.scene=yade.cld.CLDemRun.clDemToYade(sim,stepPeriod=20,relTol=1e-5)
	#O.timingEnabled=True
	# remove last engine and the clDem field
	#O.scene.engines=O.scene.engines[0:-1]
	#O.scene.fields=[O.scene.fields[0]]
else: # run via OpenCL only
	O.scene.fields=[yade.cld.CLDemField(sim)]
	O.scene.engines=[yade.cld.CLDemRun(stepPeriod=1),]
#O.scene.ranges=[yade.gl.Gl1_CLDemField.parRange]
yade.gl.Gl1_CLDemField.bboxes=False

O.saveTmp()
#O.timingEnabled=True
#O.save('/tmp/a.xml')
