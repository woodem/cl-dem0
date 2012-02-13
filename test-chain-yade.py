# encoding: utf-8
import sys
sys.path.append('.')
import miniDem
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=10
supports=[0,49,];
r=.005
E=1e4
rho=1e4

sim=miniDem.Simulation(pNum,dNum,useL1Geom=True) # pass from command-line

sim.scene.materials=[miniDem.ElastMat(density=rho,young=E)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4

for i in range(0,N):
	sim.par.append(miniDem.mkSphere((2*r*i,0,0),r,sim,matId=0,fixed=(i in supports)))
	print '#%d, flags=%d'%(i,sim.par[-1].flags)
	if i>0: sim.con.append(miniDem.Contact(ids=(i-1,i)))

sim.scene.dt=.2*sim.pWaveDt()

miniDem.briefOutput()
sim.show()

def yadeCopy():
	import _miniDem
	import yade.utils
	# copy materials
	yadeMat=[]
	for m in sim.scene.materials:
		if isinstance(m,_miniDem.ElastMat):
			yadeMat.append(yade.dem.FrictMat(young=m.young,density=m.density,tanPhi=0))
		else: yadeMat.append(None)
	# copy particles
	for p in sim.par:
		if isinstance(p.shape,_miniDem.Sphere):
			yade.O.dem.par.append(yade.utils.sphere(p.pos,p.shape.radius,material=yadeMat[p.matId],fixed=(p.dofs==0),color=.3,wire=True))
	ids1, ids2=[c.ids[0] for c in sim.con],[c.ids[1] for c in sim.con]
	yade.utils.createContacts(ids1,ids2,[yade.dem.Cg2_Sphere_Sphere_L6Geom()],[yade.dem.Cp2_FrictMat_FrictPhys(ktDivKn=0)])
	yade.dem.BoundDispatcher([yade.dem.Bo1_Sphere_Aabb()])()
	yade.O.dem.par.append([yade.utils.sphere(p.pos,p.shape.radius,material=yadeMat[p.matId],fixed=True,color=.8) for p in sim.par])
	yade.O.scene.dt=sim.scene.dt
	yade.O.scene.engines=[
		yade.core.PyRunner('if O.scene.step==0:\n\tfor i in range(0,len(sim.par)): sim.par[i].pos=O.dem.par[i+len(sim.par)].pos'),
		yade.core.PyRunner('sim.scene.dt=O.scene.dt'),
		yade.dem.Gravity(gravity=(sim.scene.gravity)),
		yade.dem.Leapfrog(damping=sim.scene.damping),
		yade.dem.ForceResetter(),
		yade.dem.ContactLoop([yade.dem.Cg2_Sphere_Sphere_L6Geom()],[yade.dem.Cp2_FrictMat_FrictPhys()],[yade.dem.Law2_L6Geom_FrictPhys_LinEl6(charLen=float('inf'))]),
		yade.dem.IntraForce([yade.dem.In2_Sphere_ElastMat()]),
		yade.core.PyRunner('sim.run(1)',1), # run the OpenCL code
		yade.core.PyRunner('for i in range(0,len(sim.par)): yade.O.dem.par[i+len(sim.par)].pos=sim.par[i].pos',1)
	]
	O.saveTmp()
	yade.qt.Inspector()
	yade.qt.View()
	yade.qt.Controller()
	
yadeCopy()

