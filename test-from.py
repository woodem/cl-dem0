import sys
sys.path.append('.')
import clDem

from miniEigen import *

from woo import utils
m=utils.defaultMaterial()
O.dem.par.append([
	utils.wall((0,0,0),axis=2,material=m,fixed=True),
	utils.sphere((0,0,1.3),radius=1,material=m)
])
O.dem.par[-1].vel=(-1.,0,0)
O.dem.par[-1].angVel=(0,1,0)
# O.dem.par[-1].material.tanPhi=0. # no friction
O.scene.dt=.003*utils.pWaveDt()
O.scene.trackEnergy=True
from woo.dem import *
from woo.core import*
O.scene.engines=[
	#PyRunner('if len(O.dem.con)>0 and O.dem.con[0].real: O.pause()'),
	Gravity(gravity=(0,0,-10)),
	Leapfrog(damping=.05,reset=True),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()]),
	ContactLoop([Cg2_Sphere_Sphere_L6Geom(),Cg2_Wall_Sphere_L6Geom()],[Cp2_FrictMat_FrictPhys()],[Law2_L6Geom_FrictPhys_IdealElPl()],applyForces=True),
]
import woo.qt
woo.qt.View()

O.scene.clDev=(1,0) # intel
sim=woo.cld.CLDemField.wooToClDem(O.scene,stepPeriod=1,relTol=-1)

#O.scene.engines[-1].raiseLimit=1e4

#from woo import plot
#plot.plots=
woo.gl.Gl1_CLDemField.parWire=True

O.saveTmp()

