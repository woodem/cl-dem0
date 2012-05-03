# encoding: utf-8
import sys
sys.path.append('.')
import clDem

from miniEigen import *
from math import *
import pylab, itertools, random

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

sim=clDem.Simulation(pNum,dNum,breakTension=True,trackEnergy=True)

r=0.1
sim.scene.materials=[clDem.ElastMat(young=1e6,density=1e3)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4
sim.scene.verletDist=.2*r # collision detection in this case
sim.maxScheduledSteps=1

sim.par.append(clDem.mkWall(pos=(-10,-10,5),axis=2,sim=sim,matId=0))
sim.par.append(clDem.mkSphere(pos=(2,2,5+1.1*r),radius=r,sim=sim,matId=0))
sim.par[-1].vel=(0,0,-.0005)

sim.scene.dt=.2*sim.pWaveDt()

sim.run(1)
#for i in range(0,len(sim.par)):
#	print i,sim.getBbox(i)

if 1:
	from yade import *
	import yade.cld
	import yade.log
	import yade.gl
	#yade.log.setLevel('GLViewer',yade.log.TRACE)
	O.scene.fields=[yade.cld.CLDemField(sim)]
	nan=float('nan')
	O.scene.loHint=O.scene.hiHint=(nan,nan,nan)
	O.scene.engines=[yade.cld.CLDemRun(stepPeriod=10),]
	O.scene.ranges=[yade.gl.Gl1_CLDemField.parRange]
	yade.gl.Gl1_CLDemField.bboxes=False



