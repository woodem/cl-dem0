# encoding: utf-8
import sys
sys.path.append('.')
import clDem

from miniEigen import *
from math import *
import pylab, itertools, random

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

dim=40,40,10
margin=10
r=.005

sim=clDem.Simulation(pNum,dNum,"-DL6GEOM_BREAK_TENSION") # -DTRACK_ENERGY")

sim.scene.materials=[clDem.ElastMat(young=1e6,density=1e3)]
sim.scene.gravity=(-4,-5,-10)
sim.scene.damping=.4
sim.scene.verletDist=.3*r # collision detection in this case
sim.maxScheduledSteps=10

sim.par.append(clDem.mkWall(pos=(0,0,2*r),axis=2,sim=sim,matId=0))
sim.par.append(clDem.mkWall(pos=(0,0,0),axis=0,sim=sim,matId=0))
sim.par.append(clDem.mkWall(pos=(0,0,0),axis=1,sim=sim,matId=0))

# ground
#for x0,x1 in itertools.product(range(0,dim[0]),range(0,dim[1])):
#	sim.par.append(clDem.mkSphere((x0*2*r,x1*2*r,0),r,sim,matId=0,groups=0b011,fixed=True))

sim.scene.loneGroups=0b010

import yade.pack
sp=yade.pack.SpherePack()
sp.makeCloud((2*r*margin,2*r*margin,2*r),((dim[0]-margin)*2*r,(dim[1]-margin)*2*r,dim[2]*2*r),r,rRelFuzz=.5)
for center,radius in sp:
	sim.par.append(clDem.mkSphere(center,radius,sim,matId=0,groups=0b001,fixed=False))
	sim.par[-1].vel=(0,0,-.05)

sim.scene.dt=.5*sim.pWaveDt()
O.scene.dt=sim.scene.dt


from yade import *
import yade.cld
import yade.log
import yade.gl
#yade.log.setLevel('GLViewer',yade.log.TRACE)
clf=yade.cld.CLDemField()
clf.sim=sim
O.scene.fields=[clf]
nan=float('nan')
O.scene.loHint=O.scene.hiHint=(nan,nan,nan)
O.scene.engines=[yade.cld.CLDemRun(stepPeriod=10),]
O.scene.ranges=[yade.gl.Gl1_CLDemField.parRange]
yade.gl.Gl1_CLDemField.bboxes=False


#O.step()

if 0:
	for i in range(0,):
		sim.run(1)
		zCoord.append(sim.par[1].pos[2])
		f1z.append(sim.par[1].force[2])
