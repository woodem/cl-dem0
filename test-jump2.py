# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from minieigen import *
from math import *
import pylab
import itertools, random

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=100
r=.005

sim=clDem.Simulation(pNum,dNum,breakTension=True,trackEnergy=True)

sim.scene.materials=[clDem.ElastMat(young=1e6,density=1e4)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.1
sim.scene.verletDist=r # collision detection in this case
sim.maxScheduledSteps=10

dim=40,40

for x0,x1 in itertools.product(range(0,dim[0]),range(0,dim[1])):
	sim.par.append(clDem.mkSphere((x0*2*r,x1*2*r,0),r,sim,matId=0,fixed=True))
for i in range(0,N):
	sim.par.append(clDem.mkSphere((random.random()*dim[0]*2*r,random.random()*dim[1]*2*r,2*r+.1*i*2*r),r,sim,matId=0,fixed=False))
	sim.par[-1].vel=(0,0,random.random()*-.05)
sim.scene.dt=.2*sim.pWaveDt()

#sim.run(1)
#for i in range(0,500):
#	sim.run(10) # do not reset arrays anymore
#	print 'E',sim.scene.energyTotal(),sim.scene.energyError();
#	print 'saved',sim.saveVtk('/tmp/jump')
import woo, woo.cld
S=woo.master.scene=woo.cld.CLDemField.clDemToWoo(sim,100,relTol=0)

