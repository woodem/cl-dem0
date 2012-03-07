# encoding: utf-8
import sys
sys.path.append('.')
import clDem
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=5 # length of chain
supports=[0,9,];
r=.005
E=1e4
rho=1e4

checkEqual=True # write message to the console if position/orientation/uN differs

useL1Geom,ktDivKn,charLen,dtFrac=(
	(True,0,None,.2),
	(False,.2,None,.05),
	(False,.8,1e5*r,.01),
)[2] # pick configuration set here

sim=clDem.Simulation(pNum,dNum,' '.join([
	'-DGEOM_L1GEOM' if useL1Geom else '',
	'-DBEND_CHARLEN=%g'%charLen if charLen else '',
	'-DSHEAR_KT_DIV_KN=%g'%ktDivKn if ktDivKn>0 else '',
	'-DTRACK_ENERGY',
])) # pass from command-line

sim.scene.materials=[clDem.ElastMat(density=rho,young=E)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4
sim.scene.verletDist=.3*r
cons=[]

for i in range(0,N):
	sim.par.append(clDem.mkSphere((2*r*i*.99999,0,0),r,sim,matId=0,fixed=(i in supports)))
	print '#%d, flags=%d'%(i,sim.par[-1].flags)
	if i>0: cons.append((i-1,i))

sim.scene.dt=dtFrac*sim.pWaveDt()
yade.O.scene.dt=sim.scene.dt

demField=yade.dem.DemField()
clField=yade.cld.CLDemField(sim)
yade.O.scene.fields=[demField,clField]

clDem.briefOutput()
sim.show()

def yadeCopy():
	import _clDem
	import yade.utils
	# copy materials
	yadeMat=[]
	for m in sim.scene.materials:
		if isinstance(m,_clDem.ElastMat):
			yadeMat.append(yade.dem.FrictMat(young=m.young,density=m.density,tanPhi=float('inf')))
		else: yadeMat.append(None)
	# copy particles
	for p in sim.par:
		if isinstance(p.shape,_clDem.Sphere):
			yade.O.dem.par.append(yade.utils.sphere(p.pos,p.shape.radius,material=yadeMat[p.matId],fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
		elif isinstance(p.shape,_clDem.Wall):
			yade.O.dem.par.append(yade.utils.wall(p.pos,p.shape.axis,material=yadeMat[p.matId],fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
	# copy contacts
	#ids1, ids2=[c[0] for c in cons],[c[1] for c in cons]
	#yade.utils.createContacts(ids1,ids2,[yade.dem.Cg2_Sphere_Sphere_L6Geom()],[yade.dem.Cp2_FrictMat_FrictPhys()])
	yade.O.dem.collectNodes()
	yade.O.scene.dt=sim.scene.dt
	yade.O.scene.trackEnergy=True
	yade.O.scene.engines=[
		yade.dem.Gravity(gravity=(sim.scene.gravity),field=demField),
		yade.dem.Leapfrog(damping=sim.scene.damping,reset=True,kinSplit=True,field=demField),
		yade.dem.InsertionSortCollider([yade.dem.Bo1_Sphere_Aabb()],field=demField),
		yade.dem.ContactLoop([yade.dem.Cg2_Sphere_Sphere_L6Geom()],[yade.dem.Cp2_FrictMat_FrictPhys(ktDivKn=ktDivKn)],[yade.dem.Law2_L6Geom_FrictPhys_LinEl6(charLen=float('inf') if (useL1Geom or not charLen) else charLen)],field=demField),
		yade.dem.IntraForce([yade.dem.In2_Sphere_ElastMat()],field=demField),
		yade.cld.CLDemRun(stepPeriod=1,compare=True,relTol=1e-4,field=clField),
		#yade.core.PyRunner('showEnergies()',100),
	]
	yade.qt.View()
	#yade.qt.Renderer().bound=True
	yade.qt.Controller()

def showEnergies():
	print '$/* Σ = %g / %g, ε = %g / %g'%(sim.scene.energyTotal(),O.scene.energy.total(),sim.scene.energyError(),O.scene.energy.relErr())
	eg,ec=sim.scene.energy,dict(O.scene.energy)
	for g,c in (('Ekt','kinTrans'),('Ekr','kinRot'),('elast','elast'),('damp','nonviscDamp'),('grav','grav')):
		gg='%g'%eg[g] if g in eg.keys() else '-';
		cc='%g'%ec[c] if c in ec.keys() else '-';
		if gg=='-' and cc=='-': continue
		print '\t%10s %15s / %15s'%(g,gg,cc)

	
yadeCopy()

