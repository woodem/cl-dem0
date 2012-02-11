# encoding: utf-8
import sys
sys.path.append('.')
import miniDem
from miniEigen import *
from math import *
import itertools

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

M,N=50,50
r=.005
E=1e4
rho=1e4

sim=miniDem.Simulation(pNum,dNum) # pass from command-line

sim.scene.materials=[miniDem.ElastMat(young=1e4)]
sim.scene.dt=.2*r/sqrt(E/rho)
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.2

for n,m in itertools.product(range(0,N),range(0,M)): # m advances the fastest
	isSupp=(m==0 or n==0 or (m==M-1 and n<N/2))
	sim.par.append(miniDem.mkSphere(pos=(2*r*m,2*r*n,0),radius=r,rho=rho,matId=0,fixed=isSupp))
	pid=len(sim.par)-1; assert(pid==n*M+m)
	if n>0: sim.con.append(miniDem.Contact(ids=(pid,pid-M))) # contact below
	if m>0: sim.con.append(miniDem.Contact(ids=(pid,pid-1))) # contact left

for i in range(0,300):
	sim.run(5)
	print 'Saved',sim.saveVtk('/tmp/cloth',)

