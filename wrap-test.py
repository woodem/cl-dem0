import sys
sys.path.append('.')
import miniDem
sim=miniDem.Simulation()
print sim.scene
print sim.con
print sim.par
par=miniDem.Particle()
print par.shape
par.shape=miniDem.Sphere(radius=4)
print par.shape
sim.par=[par,]
print sim.par
print sim.scene.materials
sim.scene.materials=[]
sim.scene.materials=[miniDem.ElastMat()]
print sim.par[0].shape.radius
c=miniDem.Contact()
