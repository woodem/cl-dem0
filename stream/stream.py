import sys
sys.path.append('..')
import clDem

from miniEigen import *
from math import *

# set platform and device numbers here
# if left at -1, first platform and device will be used
# you will be notified in the terminal as well
sim=clDem.Simulation(platformNum=1,deviceNum=0,breakTension=True,ktDivKn=.2,opts="-cl-strict-aliasing -I..")
sim.scene.materials=[clDem.FrictMat(young=1e6,density=1e3,ktDivKn=.2,tanPhi=.5)]
# add boundaries
for axis in 0,1,2: sim.par.append(clDem.mkWall(pos=(0,0,0),axis=axis,sim=sim,matId=0,groups=0b011))
# load spheres from file
for l in file(sys.argv[1]):
	if l.startswith('#'): continue
	x,y,z,r=[float(n) for n in l.split()]
	sim.par.append(clDem.mkSphere((x,y,z),r,sim,matId=0,groups=0b001,fixed=False))

sim.scene.gravity=(-4,-5,-10)
sim.scene.damping=.4
sim.scene.verletDist=-.3 # negative means relative to minimum sphere size
sim.scene.loneGroups=0b010  # no contacts of walls with themselves

# maximum number of enqueued steps (times 7 kernels/step), before returning back to CPU
#sim.maxScheduledSteps=20

if 1:
	sim.run(10000)
	print 'Done at step %d, bye.'%sim.scene.step
# show inside yade, for visualization
else:
	from yade import *
	import yade.cld
	import yade.gl
	O.scene.fields=[yade.cld.CLDemField(sim)]
	O.scene.engines=[yade.cld.CLDemRun(stepPeriod=10),]
	O.scene.ranges=[yade.gl.Gl1_CLDemField.parRange]
	yade.gl.Gl1_CLDemField.bboxes=False

