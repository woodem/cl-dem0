# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=20 # length of chain
supports=[0,9,19];
r=.005
E=1e4
rho=1e4
nan=float('nan')

useL1Geom,ktDivKn,charLen,dtFrac=(
	(False,.2,nan,.05),
	(False,.8,1e5*r,.01),
)[1] # pick configuration set here

sim=clDem.Simulation(pNum,dNum,charLen=charLen,ktDivKn=ktDivKn,trackEnergy=True)
sim.scene.materials=[clDem.ElastMat(density=rho,young=E)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4
sim.scene.verletDist=.3*r

for i in range(0,N):
	sim.par.append(clDem.mkSphere((2*r*i*.99999,0,0),r,sim,matId=0,fixed=(i in supports)))

sim.scene.dt=dtFrac*sim.pWaveDt()

def wooCopy():
	import _clDem
	import woo.utils
	demField=woo.dem.DemField()
	clField=woo.cld.CLDemField(sim)
	woo.O.scene=woo.core.Scene() # reinitialize everything
	woo.O.scene.fields=[demField,clField]
	# copy materials
	wooMat=[]
	for m in sim.scene.materials:
		if isinstance(m,_clDem.ElastMat):
			wooMat.append(woo.dem.FrictMat(young=m.young,density=m.density,tanPhi=float('inf'),poisson=float('nan')))
		else: wooMat.append(None)
	# copy particles
	for p in sim.par:
		if isinstance(p.shape,_clDem.Sphere):
			woo.O.dem.par.append(woo.utils.sphere(p.pos,p.shape.radius,material=wooMat[p.matId],fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
		elif isinstance(p.shape,_clDem.Wall):
			woo.O.dem.par.append(woo.utils.wall(p.pos,p.shape.axis,material=wooMat[p.matId],fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
	woo.O.dem.collectNodes()
	woo.O.scene.dt=sim.scene.dt
	woo.O.scene.trackEnergy=True
	woo.O.scene.engines=[
		woo.dem.Gravity(gravity=(sim.scene.gravity)),
		woo.dem.Leapfrog(damping=sim.scene.damping,reset=True,kinSplit=True),
		woo.dem.InsertionSortCollider([woo.dem.Bo1_Sphere_Aabb(),woo.dem.Bo1_Wall_Aabb(),woo.dem.Bo1_Facet_Aabb()]),
		woo.dem.ContactLoop([woo.dem.Cg2_Sphere_Sphere_L6Geom(),woo.dem.Cg2_Wall_Sphere_L6Geom(),woo.dem.Cg2_Facet_Sphere_L6Geom()],[woo.dem.Cp2_FrictMat_FrictPhys(ktDivKn=ktDivKn)],[woo.dem.Law2_L6Geom_FrictPhys_LinEl6(charLen=float('inf') if (useL1Geom or not charLen) else charLen)],applyForces=True),
		woo.cld.CLDemRun(stepPeriod=1,compare=True,relTol=1e-5,raiseLimit=10.),
	]

def showEnergies():
	print '$/* Σ = %g / %g, ε = %g / %g'%(sim.scene.energyTotal(),O.scene.energy.total(),sim.scene.energyError(),O.scene.energy.relErr())
	eg,ec=sim.scene.energy,dict(O.scene.energy)
	for g,c in (('Ekt','kinTrans'),('Ekr','kinRot'),('elast','elast'),('damp','nonviscDamp'),('grav','grav')):
		gg='%g'%eg[g] if g in eg.keys() else '-';
		cc='%g'%ec[c] if c in ec.keys() else '-';
		if gg=='-' and cc=='-': continue
		print '\t%10s %15s / %15s'%(g,gg,cc)

# inspect XML files to see that they are identical
woo.O.scene=woo.cld.CLDemRun.clDemToYade(sim,1,relTol=1e-5)
O.save('/tmp/a1.xml')
wooCopy() # this sets O.scene inside
O.save('/tmp/a2.xml')


woo.qt.View()
woo.qt.Renderer().bound=True
woo.qt.Controller()
