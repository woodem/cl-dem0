# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from miniEigen import *
from math import *
import pylab

N=50
supports=[0,49,];
r=.005
E=1e4
rho=1e3

sim=clDem.Simulation(breakTension=True,trackEnergy=True,opts='')
sim.scene.materials=[clDem.ElastMat(young=1e4,density=1e4)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=0.
sim.scene.dt=-.01

sim.par.append(clDem.mkSphere((0,0,0),.005,sim,matId=0,fixed=True))
sim.par.append(clDem.mkSphere((0,0,.015),.005,sim,matId=0,fixed=False))

import woo.cld
woo.O.scene=woo.cld.CLDemRun.clDemToYade(sim,stepPeriod=1,relTol=1e-8)
#pylab.plot(zCoord,label='z-coord'); pylab.legend(loc='lower left'); pylab.twinx(); pylab.plot(f1z,c='red',label='F1z'); pylab.legend(); pylab.show()
