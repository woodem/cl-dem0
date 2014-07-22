import sys
sys.path.append('.')
import clDem
import minieigen
sim=clDem.Simulation() #platformNum,deviceNum
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
par=clDem.Particle()
show(par)
par.shape=clDem.Sphere(radius=4)
show(par)
sim.par=[par,]
print sim.par
print sim.scene.materials
sim.scene.materials=[]
sim.scene.materials=[clDem.ElastMat()]
print sim.par[0].shape.radius
c=clDem.Contact()
show(c)
