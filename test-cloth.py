# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from minieigen import *
from math import *
import itertools

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

M,N=100,100
r=.005
E=1e4
rho=1e4

sim=clDem.Simulation(pNum,dNum) # pass from command-line

sim.scene.materials=[clDem.ElastMat(young=E,density=rho)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.2
sim.scene.verletDist=float('nan') # no collision detection

for n,m in itertools.product(range(0,N),range(0,M)): # m advances the fastest
	isSupp=(m==0 or n==0 or (m==M-1 and n<N/2))
	sim.par.append(clDem.mkSphere((2*r*m*.999999,2*r*n*.9999999,0),r,sim,matId=0,fixed=isSupp))
	pid=len(sim.par)-1; assert(pid==n*M+m)
	if n>0: sim.con.append(clDem.Contact(ids=(pid,pid-M))) # contact below
	if m>0: sim.con.append(clDem.Contact(ids=(pid,pid-1))) # contact left

import woo, woo.cld, woo.gl
S=woo.master.scene=woo.cld.CLDemField.clDemToWoo(sim,stepPeriod=100,relTol=-1e-3)
#S.engines=S.engines[:-1]
#S.fields=S.fields[:-1]

S.ranges=[woo.gl.Gl1_CLDemField.parRange,woo.gl.Gl1_CLDemField.conRange]
S.dt=sim.scene.dt=.2*sim.pWaveDt();

woo.gl.Gl1_CLDemField.parRange.label='|v|'
woo.gl.Gl1_CLDemField.conRange.label='Fn'


woo.master.timingEnabled=True

S.saveTmp()


#for i in range(0,300):
#	sim.run(5)
#	print 'Saved',sim.saveVtk('/tmp/cloth',ascii=False,compress=True)

