import sys
sys.path.append('..')
import clDem

from miniEigen import *
from math import *

# set platform and device numbers here
# if left at -1, first platform and device will be used
# you will be notified in the terminal as well
sim=clDem.Simulation(platformNum=2,deviceNum=0,breakTension=True,ktDivKn=.2,opts="-cl-strict-aliasing -I..")
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

sim.showND=0 # show NDRange for all kernels every 20 steps; set to 0 to disable

# maximum number of enqueued steps (times 7 kernels/step), before returning back to CPU
sim.maxScheduledSteps=100

if 1:
	sim.run(2000)
	print 'Done at step %d, bye.'%sim.scene.step
# show inside woo, for visualization
else:
	from woo import *
	import woo.cld
	import woo.gl
	O.scene.fields=[woo.cld.CLDemField(sim)]
	O.scene.engines=[woo.cld.CLDemRun(stepPeriod=40),]
	O.scene.ranges=[woo.gl.Gl1_CLDemField.parRange]
	woo.gl.Gl1_CLDemField.bboxes=False

