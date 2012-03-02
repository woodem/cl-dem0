# encoding: utf-8
import sys
sys.path.append('.')
import clDem

from miniEigen import *
from math import *
import pylab, itertools, random

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=50
dim=40,40,5
r=.0051

sim=clDem.Simulation(pNum,dNum,"-DL6GEOM_BREAK_TENSION -DTRACK_ENERGY")

sim.scene.materials=[clDem.ElastMat(young=1e6,density=1e3)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4
sim.scene.verletDist=.2*r # collision detection in this case
sim.maxScheduledSteps=5

# ground
for x0,x1 in itertools.product(range(0,dim[0]),range(0,dim[1])):
	sim.par.append(clDem.mkSphere((x0*2*r,x1*2*r,0),r,sim,matId=0,fixed=True))

import yade.pack
sp=yade.pack.SpherePack()
sp.makeCloud((0,0,2*r),(dim[0]*2*r,dim[1]*2*r,dim[2]*2*r),r,rRelFuzz=.5)
for center,radius in sp:
	sim.par.append(clDem.mkSphere(center,radius,sim,matId=0,fixed=False))
	sim.par[-1].vel=(0,0,random.random()*-.05)

sim.scene.dt=.2*sim.pWaveDt()
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
O.scene.engines=[yade.cld.CLDemRun(stepPeriod=200),]
O.scene.ranges=[yade.gl.Gl1_CLDemField.parRange]
#yade.gl.Gl1_CLDemField.par=False


#O.step()

if 0:
	for i in range(0,):
		sim.run(1)
		zCoord.append(sim.par[1].pos[2])
		f1z.append(sim.par[1].force[2])
