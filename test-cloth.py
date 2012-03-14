# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from miniEigen import *
from math import *
import itertools

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

M,N=50,50
r=.005
E=1e4
rho=1e4

sim=clDem.Simulation(pNum,dNum) # pass from command-line

sim.scene.materials=[clDem.ElastMat(young=E,density=rho)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.2
sim.scene.verletDist=-.1*r # no collision detection

for n,m in itertools.product(range(0,N),range(0,M)): # m advances the fastest
	isSupp=(m==0 or n==0 or (m==M-1 and n<N/2))
	sim.par.append(clDem.mkSphere((2*r*m*.999999,2*r*n*.999999,0),r,sim,matId=0,fixed=isSupp))
	pid=len(sim.par)-1; assert(pid==n*M+m)
	if n>0: sim.con.append(clDem.Contact(ids=(pid,pid-M))) # contact below
	if m>0: sim.con.append(clDem.Contact(ids=(pid,pid-1))) # contact left

yade.O.scene.dt=sim.scene.dt=.2*sim.pWaveDt();

import yade.cld
yade.gl.Gl1_CLDemField.parRange.label='|v|'
yade.gl.Gl1_CLDemField.conRange.label='Fn'

O.timingEnabled=True
import yade.cld
yade.O.scene=yade.cld.CLDemRun.clDemToYade(sim,stepPeriod=1,relTol=-1e-3)
O.scene.ranges=[yade.gl.Gl1_CLDemField.parRange,yade.gl.Gl1_CLDemField.conRange]


#for i in range(0,300):
#	sim.run(5)
#	print 'Saved',sim.saveVtk('/tmp/cloth',ascii=False,compress=True)

