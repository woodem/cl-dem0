# encoding: utf-8
import sys
sys.path.append('.')
import miniDem
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=10 # length of chain
supports=[0,9,];
r=.005
E=1e4
rho=1e4

checkEqual=False # write message to the console if position/orientation/uN differs

useL1Geom,ktDivKn,charLen,dtFrac=(
	(True,0,None,.2),
	(False,.2,None,.05),
	(False,.8,1e5*r,.01),
)[1] # pick configuration set here

sim=miniDem.Simulation(pNum,dNum,' '.join([
	'-DGEOM_L1GEOM' if useL1Geom else '',
	'-DBEND_CHARLEN=%g'%charLen if charLen else '',
	'-DSHEAR_KT_DIV_KN=%g'%ktDivKn if ktDivKn>0 else '',
	'-DTRACK_ENERGY',
])) # pass from command-line

sim.scene.materials=[miniDem.ElastMat(density=rho,young=E)]
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4

for i in range(0,N):
	sim.par.append(miniDem.mkSphere((2*r*i,0,0),r,sim,matId=0,fixed=(i in supports)))
	print '#%d, flags=%d'%(i,sim.par[-1].flags)
	if i>0: sim.con.append(miniDem.Contact(ids=(i-1,i)))

sim.scene.dt=dtFrac*sim.pWaveDt()

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
			yade.O.dem.par.append(yade.utils.sphere(p.pos,p.shape.radius,material=yadeMat[p.matId],fixed=(p.dofs==0),color=.5 if p.dofs==0 else .3,wire=True))
	ids1, ids2=[c.ids[0] for c in sim.con],[c.ids[1] for c in sim.con]
	yade.utils.createContacts(ids1,ids2,[yade.dem.Cg2_Sphere_Sphere_L6Geom()],[yade.dem.Cp2_FrictMat_FrictPhys(ktDivKn=0 if useL1Geom else ktDivKn)]) # .2 is constant in the OpenCL code for this case
	yade.dem.BoundDispatcher([yade.dem.Bo1_Sphere_Aabb()])() #create bboxes to give hint to OpenGL on displaying particles
	yade.O.dem.par.append([yade.utils.sphere(p.pos,p.shape.radius,material=yadeMat[p.matId],fixed=True,color=1 if p.dofs==0 else .8) for p in sim.par]) # add particles from OpenCL; will not be moved by the integrator
	for i in range(0,len(sim.par)): yade.O.dem.par[len(sim.par)+i].shape.nodes[0].dem.energySkip=True
	yade.O.scene.dt=sim.scene.dt
	yade.O.scene.trackEnergy=True
	yade.O.dem.collectNodes()
	yade.O.scene.engines=[
		# restore OpenCL sim after reload
		yade.core.PyRunner('if O.scene.step==0:\n\tfor i in range(0,len(sim.par)): sim.par[i].pos,sim.par[i].ori,sim.par[i].vel,sim.par[i].angVel=O.dem.par[i+len(sim.par)].pos,O.dem.par[i+len(sim.par)].ori,O.dem.par[i+len(sim.par)].vel,O.dem.par[i+len(sim.par)].angVel\n\tfor c in sim.con: c.geom,c.phys=None,None\n\tsim.scene.energyReset()'),
		yade.core.PyRunner('sim.scene.dt=O.scene.dt'), # adjust in case it changes meanwhile
		yade.dem.Gravity(gravity=(sim.scene.gravity)),
		yade.dem.Leapfrog(damping=sim.scene.damping,reset=True,kinSplit=True),
		yade.dem.ContactLoop([yade.dem.Cg2_Sphere_Sphere_L6Geom()],[yade.dem.Cp2_FrictMat_FrictPhys()],[yade.dem.Law2_L6Geom_FrictPhys_LinEl6(charLen=float('inf') if (useL1Geom or not charLen) else charLen)]),
		yade.dem.IntraForce([yade.dem.In2_Sphere_ElastMat()]),
		yade.core.PyRunner('sim.run(1)',1), # run the OpenCL code
		# copy the vbisible stuff into DEM
		yade.core.PyRunner('for i in range(0,len(sim.par)): yade.O.dem.par[i+len(sim.par)].pos,yade.O.dem.par[i+len(sim.par)].ori=sim.par[i].pos,sim.par[i].ori',1),
	]+([
		#check normal displacement difference
		yade.core.PyRunner('for c in sim.con:\n\tcDem=O.dem.con[c.ids]\n\tduN=abs(c.geom.uN-cDem.geom.uN)\n\tif duN/r>1e-6: print "##%d+%d: ΔuN=%g"%(c.ids[0],c.ids[1],duN)\n\t#print cDem.geom.node.ori*Vector3.UnitX, c.ori.row(0)'),
		#check position difference
		yade.core.PyRunner('for i in range(0,len(sim.par)):\n\tdx=(yade.O.dem.par[i].pos-sim.par[i].pos).norm()\n\tif dx/r>1e-5: print "#%d Δx=%g"%(i,dx)'),
		#check orientation difference
		yade.core.PyRunner('for i in range(0,len(sim.par)):\n\tdo=(yade.O.dem.par[i].ori-sim.par[i].ori).norm()\n\tif do>1e-5: print "#%d Δori=%g"%(i,do)'),
	] if checkEqual else [])+[
		#yade.core.PyRunner('for c in yade.O.dem.con: print "* ##%d+%d: uN=%g, fN=%g"%(c.id1,c.id2,c.geom.uN,c.phys.force[0])'),
		#yade.core.PyRunner('for c in sim.con: print "$ ##%d+%d: uN=%g, fN=%g"%(c.ids[0],c.ids[1],c.geom.uN,c.force[0])'),
		yade.core.PyRunner('showEnergies()',100), #print "$E:",sim.scene.energy,sim.scene.energyTotal(),sim.scene.energyError()\nprint "*E:",dict(O.scene.energy),O.scene.energy.total(),O.scene.energy.relErr()'),
	]
	O.saveTmp()
	#yade.qt.Inspector()
	yade.qt.View()
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

