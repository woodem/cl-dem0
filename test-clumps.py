
import sys
sys.path.append('.')
import clDem

from miniEigen import *
from math import *
import pylab, itertools, random

import yade.core

sim=clDem.Simulation(breakTension=True)

sim.scene.materials=[clDem.ElastMat(young=1e6,density=1e3)] #,ktDivKn=.2,tanPhi=.5)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.1
sim.maxScheduledSteps=10
sim.scene.dt=-.3

# stand-alone sphere
# sim.par.append(clDem.mkSphere([0,0,0],.5,sim=sim))
# clumps
relPos=[(0,-.5,-.5),(0,.5,0),(.5,0,0),(0,0,.5)]
coords=[(-2,0,0),] #,(2,0,0),(0,2,0),(0,-2,0)]
for i,cc in enumerate(coords):
	# This shorthand command does something like this:
	# O.bodies.appendClumped([utils.sphere(...),utils.sphere(...),utils.sphere(...)])
	# and returns tuple of clumpId,[bodyId1,bodyId2,bodyId3]
	clump,spheres=sim.addClump([clDem.mkSphere([relPos[j][0]+coords[i][0],relPos[j][1]+coords[i][1],relPos[j][2]+coords[i][2]],.5,sim=sim) for j in range(0,i+3)])
	print clump,spheres
sim.par.append(clDem.mkWall(pos=(0,0,-1.3),axis=2,sim=sim,matId=0,groups=0b011))

# visualization
#O.scene.fields=[yade.cld.CLDemField(sim)]
#O.scene.engines=[yade.cld.CLDemRun(stepPeriod=1),]

O.scene=yade.cld.CLDemField.clDemToYade(sim,stepPeriod=1,relTol=-1e-5)
#print O.scene.engines[-1].raiseLimit
#O.scene.engines[-1].raiseLimit=1e11
#O.scene.engines=O.scene.engines+[yade.core.PyRunner('import yade; sim=yade.O.scene.fields[-1].sim; print yade.O.dem.par[1].pos,sim.par[1].pos')]
#globals()['sim']=sim

from yade import qt
qt.View()
O.run(60,True)
O.saveTmp()
