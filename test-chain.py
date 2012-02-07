# encoding: utf-8
import miniDem
import sys
sys.path.append('/usr/local/lib/yade-tr2/py')
from miniEigen import *
from math import *

pNum=int(sys.argv[1]) if len(sys.argv)>1 else -1
dNum=int(sys.argv[2]) if len(sys.argv)>2 else -1

N=10
supports=[0,9,];
r=.005
E=1e9
rho=1e4
mass=rho*(4/3.)*pi*r**3
inertia=mass*(2/5.)*r**2*Vector3.Ones

sim=miniDem.Simulation(pNum) # pass from command-line

sim.scene.materials=[miniDem.ElastMat(young=E)]
sim.scene.dt=.2*r/sqrt(E/rho)
sim.scene.gravity=(0,0,-10)
sim.scene.damping=.4

for i in range(0,N):
	isSupport=(i in supports)
	p=miniDem.Particle(pos=(2*r*i,0,0),mass=mass,inertia=inertia,dofs=0 if isSupport else 63,shape=miniDem.Sphere(radius=r),matId=0)
	print '#%d, flags=%d'%(i,p.flags)
	sim.par.append(p)
	if i>0:
		sim.con.append(miniDem.Contact(ids=(i-1,i)))

sim.run(10000)

# output similar to the c++ code
Vector3.__str__=lambda s: '('+','.join('%g'%e for e in s)+')'
Matrix3.__str__=lambda s: '('+', '.join(','.join('%g'%e for e in s.row(i)) for i in (0,1,2))+')'


def showSim(sim):
	print 'At step %d (t=%g), Δt=%g'%(sim.scene.step,sim.scene.t,sim.scene.dt)
	for i,p in enumerate(sim.par):
		print '#%d x=%s; v=%s'%(i,p.pos,p.vel)
	for c in sim.con:
		print '##%d+%d: p=%s'%(c.ids[0],c.ids[1],c.pos)
		print '\tgeomT=%d, uN=%g, rot=%s'%(c.geomT,c.geom.uN,c.ori)
		print '\tphysT=%d, kN=%g, F=%s'%(c.physT,c.phys.kN,c.force)

showSim(sim)
