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
sim.scene.verletDist=2*r # collision detection in this case
sim.maxScheduledSteps=5

sim.par.append(clDem.mkSphere((0,0,0),.005,sim,matId=0,fixed=True))
sim.par.append(clDem.mkSphere((0,0,.015),.005,sim,matId=0,fixed=False))
# sim.par.append(clDem.mkSphere((1,1,1),0,sim,matId=0,fixed=True)) # tests thin AxBound
sim.par[-1].vel=(0,0,-.02)
if 0: # this should be handled by CpuCollider::initialStep now
	# create huge bboxes
	sim.bboxes=[-1,-1,-1,1,1,1,-1,-1,-1,1,1,1]
	sim.con.append(clDem.Contact(ids=(0,1)))
sim.scene.dt=.2*sim.pWaveDt()

sim.run(1)
zCoord,f1z=[],[]
for i in range(0,100):
	sim.run(10,False) # do not reset arrays anymore
	zCoord.append(sim.par[1].pos[2])
	f1z.append(sim.par[1].force[2])
	print '\t',[i.ids for i in sim.con]
	print '\tconFree',sim.conFree, 'pot',sim.pot, 'potFree',sim.potFree,'cLog',sim.cJournal
	#print sim.saveVtk('/tmp/jump')
	print 'E',sim.scene.energyTotal(),sim.scene.energyError();
pylab.plot(zCoord,label='z-coord'); pylab.legend(loc='lower left'); pylab.twinx(); pylab.plot(f1z,c='red',label='F1z'); pylab.legend(); pylab.show()
