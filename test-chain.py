# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=50
supports=[0,49,];
r=.005
E=1e4
rho=1e4

sim=clDem.Simulation(pNum,dNum) # pass from command-line

sim.scene.materials=[clDem.ElastMat(young=E,density=rho)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4
sim.maxScheduledSteps=-1 #unlimited

for i in range(0,N):
	sim.par.append(clDem.mkSphere((2*r*i,0,0),r,sim,matId=0,fixed=(i in supports)))
	print '#%d, flags=%d'%(i,sim.par[-1].flags)
	if i>0:
		sim.con.append(clDem.Contact(ids=(i-1,i)))
sim.scene.dt=.2*sim.pWaveDt()

for i in range(0,5):
	sim.run(10)
	print '10 steps done'
#	print 'Saved',sim.saveVtk('/tmp/chain',compress=False,ascii=True)

#clDem.briefOutput()
#clDem.showSim(sim)
