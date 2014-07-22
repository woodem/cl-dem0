# encoding: utf-8

#
# TODO: node-particle back-references
#
import sys
sys.path.append('.')
import clDem
import woo, woo.dem, woo.cld, woo.qt, woo.gl
from minieigen import *
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
)[0] # pick configuration set here

sim=clDem.Simulation(pNum,dNum,charLen=charLen,ktDivKn=ktDivKn,trackEnergy=True)
sim.scene.materials=[clDem.ElastMat(density=rho,young=E)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4
sim.scene.verletDist=.3*r

for i in range(0,N):
	sim.par.append(clDem.mkSphere((2*r*i*.9,0,0),r,sim,matId=0,fixed=(i in supports)))

sim.scene.dt=dtFrac*sim.pWaveDt()

def wooCopy(sim):
	'Python-based copy of the cldem scene *sim*; returns new Scene instance'
	import _clDem
	import woo.utils
	demField=woo.dem.DemField(gravity=sim.scene.gravity,loneMask=sim.scene.loneGroups)
	clField=woo.cld.CLDemField(sim)
	S=woo.core.Scene() # reinitialize everything
	S.fields=[demField,clField]
	# copy materials
	wooMat=[]
	for m in sim.scene.materials:
		if isinstance(m,_clDem.ElastMat):
			wooMat.append(woo.dem.FrictMat(young=m.young,density=m.density,tanPhi=float('inf'),ktDivKn=ktDivKn))
		else: wooMat.append(None)
	# copy particles
	for p in sim.par:
		if isinstance(p.shape,_clDem.Sphere):
			S.dem.par.append(woo.utils.sphere(p.pos,p.shape.radius,mat=wooMat[p.matId],mask=p.groups,fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
		elif isinstance(p.shape,_clDem.Wall):
			S.dem.par.append(woo.utils.wall(p.pos,p.shape.axis,mat=wooMat[p.matId],mask=p.groups,fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
		S.dem.par[-1].shape.nodes[0].clDem=woo.cld.CLDemData(clIx=len(S.dem.par)-1)
		S.dem.par[-1].shape.nodes[0].dem.angMom=(0,0,0) # better comparison with clDemToWoo
	S.dem.collectNodes()
	S.dt=sim.scene.dt
	S.trackEnergy=True
	S.engines=[
		woo.dem.Leapfrog(damping=sim.scene.damping,reset=True,kinSplit=True),
		woo.dem.InsertionSortCollider([woo.dem.Bo1_Sphere_Aabb(),woo.dem.Bo1_Wall_Aabb(),woo.dem.Bo1_Facet_Aabb()]),
		woo.dem.ContactLoop([woo.dem.Cg2_Sphere_Sphere_L6Geom(),woo.dem.Cg2_Wall_Sphere_L6Geom(),woo.dem.Cg2_Facet_Sphere_L6Geom()],[woo.dem.Cp2_FrictMat_FrictPhys()],[woo.dem.Law2_L6Geom_FrictPhys_LinEl6(charLen=float('inf') if (useL1Geom or not charLen) else charLen)],applyForces=True),
		woo.cld.CLDemRun(stepPeriod=1,relTol=1e-5,raiseLimit=100),
	]
	return S

def showEnergies():
	print '$/* Σ = %g / %g, ε = %g / %g'%(sim.scene.energyTotal(),S.energy.total(),sim.scene.energyError(),S.energy.relErr())
	eg,ec=sim.scene.energy,dict(S.energy)
	for g,c in (('Ekt','kinTrans'),('Ekr','kinRot'),('elast','elast'),('damp','nonviscDamp'),('grav','grav')):
		gg='%g'%eg[g] if g in eg.keys() else '-';
		cc='%g'%ec[c] if c in ec.keys() else '-';
		if gg=='-' and cc=='-': continue
		print '\t%10s %15s / %15s'%(g,gg,cc)

# inspect XML files to see that they are identical
S=woo.master.scene=woo.cld.CLDemField.clDemToWoo(sim,100,relTol=0)
S.save('/tmp/a1.xml')
S2=wooCopy(sim) # this sets O.scene inside
S2.save('/tmp/a2.xml')


woo.qt.View()
woo.gl.Renderer.bound=True
woo.qt.Controller()
S.one()
