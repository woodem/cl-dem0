# encoding: utf-8
import sys
sys.path.append('.')
import clDem

from miniEigen import *
from math import *
import pylab

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=50
supports=[0,49,];
r=.005
E=1e4
rho=1e3

sim=clDem.Simulation(pNum,dNum,"-DL6GEOM_BREAK_TENSION -DTRACK_ENERGY")

sim.scene.materials=[clDem.ElastMat(young=1e4,density=1e4)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.1
sim.scene.verletDist=.05*r # collision detection in this case
sim.maxScheduledSteps=5

sim.par.append(clDem.mkSphere((0,0,0),.005,sim,matId=0,fixed=True))
sim.par.append(clDem.mkSphere((0,0,.015),.005,sim,matId=0,fixed=False))
# sim.par.append(clDem.mkSphere((1,1,1),0,sim,matId=0,fixed=True)) # tests thin AxBound
sim.par[-1].vel=(0,0,-.02)
sim.scene.dt=.2*sim.pWaveDt()

zCoord,f1z=[],[]

from yade import *
import yade.cld
import yade.log
#yade.log.setLevel('GLViewer',yade.log.TRACE)
clf=yade.cld.CLDemField()
clf.sim=sim
O.scene.fields=[clf]
nan=float('nan')
O.scene.loHint=O.scene.hiHint=(nan,nan,nan)
O.scene.engines=[yade.cld.CLDemRun(stepPeriod=1),]
O.step()

if 0:
	for i in range(0,):
		sim.run(1)
		zCoord.append(sim.par[1].pos[2])
		f1z.append(sim.par[1].force[2])
