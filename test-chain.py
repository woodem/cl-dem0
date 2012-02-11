# encoding: utf-8
import sys
sys.path.append('.')
import miniDem
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=50
supports=[0,49,];
r=.005
E=1e4
rho=1e4
mass=rho*(4/3.)*pi*r**3
inertia=mass*(2/5.)*r**2*Vector3.Ones

sim=miniDem.Simulation(pNum,dNum) # pass from command-line

sim.scene.materials=[miniDem.ElastMat(young=E)]
sim.scene.dt=.2*r/sqrt(E/rho)
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4

for i in range(0,N):
	sim.par.append(miniDem.mkSphere(pos=(2*r*i,0,0),radius=r,rho=rho,matId=0,fixed=(i in supports)))
	print '#%d, flags=%d'%(i,sim.par[-1].flags)
	if i>0:
		sim.con.append(miniDem.Contact(ids=(i-1,i)))

for i in range(0,100):
	sim.run(10)
	print 'Saved',sim.saveVtk('/tmp/chain',compress=False,ascii=True)

miniDem.briefOutput()
miniDem.showSim(sim)
