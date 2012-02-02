import sys
sys.path.append('.')
import miniDem
sys.path.append('/usr/local/lib/yade-tr2/py')
import miniEigen
sim=miniDem.Simulation() #platformNum,deviceNum
def show(obj):
	import pprint
	d={}
	for a in dir(obj):
		if a.startswith('_'): continue
		if callable(getattr(obj,a)): continue
		d[a]=getattr(obj,a)
	pprint.pprint(d)

show(sim.scene)
print sim.con
print sim.par
par=miniDem.Particle()
show(par)
par.shape=miniDem.Sphere(radius=4)
show(par)
sim.par=[par,]
print sim.par
print sim.scene.materials
sim.scene.materials=[]
sim.scene.materials=[miniDem.ElastMat()]
print sim.par[0].shape.radius
c=miniDem.Contact()
show(c)
