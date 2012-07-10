from woo.dem import *
from woo.core import *
from miniEigen import *
from math import *
from woo import utils,pack
S=woo.master.scene=Scene(fields=[DemField()])

m=FrictMat(young=1e6,density=1e3,ktDivKn=2.,tanPhi=.5)
for axis in 0,1,2: S.dem.par.append(utils.wall((0,0,0),axis=axis,mat=m,mask=0b011))
sp=pack.SpherePack()
sp.load(sys.argv[1])
sp.cellSize=(0,0,0) # make the packing aperiodic
sp.toSimulation(S,mat=m,mask=0b001)
S.dem.collectNodes()
S.dt=.3*utils.pWaveDt()
S.engines=utils.defaultEngines(gravity=(-4,-5,-10),damping=.4,verletDist=-.3)
S.loneGroups=0b010  # no contacts of walls with themselves
S.run(2000)
# S.wait()
