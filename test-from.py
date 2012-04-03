import sys
sys.path.append('.')
import clDem

from miniEigen import *

from yade import utils
m=utils.defaultMaterial()
O.dem.par.append([
	utils.wall((0,0,0),axis=2,material=m,fixed=True),
	utils.sphere((0,0,.03),radius=0.01,material=m)
])
O.scene.dt=.3*utils.pWaveDt()
O.scene.trackEnergy=True
from yade.dem import *
from yade.core import*
O.scene.engines=[
	Gravity(gravity=(0,0,-10)),
	Leapfrog(damping=.4,reset=True),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()]),
	ContactLoop([Cg2_Sphere_Sphere_L6Geom(),Cg2_Wall_Sphere_L6Geom()],[Cp2_FrictMat_FrictPhys()],[Law2_L6Geom_FrictPhys_IdealElPl()],applyForces=True),
]
import yade.qt
yade.qt.View()



sim=yade.cld.CLDemField.yadeToClDem(O.scene,stepPeriod=1)
O.saveTmp()


