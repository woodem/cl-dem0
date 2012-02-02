
import miniDem
import sys
sys.path.append('/usr/local/lib/yade-tr2/py')
from miniEigen import *
from math import *

N=10
supports=[0,9,];
r=.005
E=1e9
rho=1e4
mass=rho*(4/3.)*pi*r**3
inertia=mass*(2/5.)*r**2*Vector3.Ones
par,con=[],[]

for i in range(0,N):
	isSupport=(i in supports)
	p=miniDem.Particle(pos=(2*r*i,0,0),mass=mass,inertia=inertia,dofs=0 if isSupport else 63,shape=miniDem.Sphere(radius=r))
	print '#%d, flags=%d'%(i,p.flags)
	par.append(p)
	if i>0:
		con.append(miniDem.Contact(ids=(i-1,i)))


sim=miniDem.Simulation()
sim.par=par
sim.con=con




