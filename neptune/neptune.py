import numpy as np
import matplotlib.pylab as plt
#import scipy as sp

gfac = 4.0*np.pi**2

def potential(x):
  return -gfac/np.sqrt(sum(x*x))  

def accel(x):
  return -gfac*x/sum(x*x)**1.5

class particle:
  # mass, position vector, velocity vector    
  def __init__(self,m,x,v):
    self.m = m
    self.x = x
    self.v = v

  def what(self):
    print( 'mass     = ', self.m)
    print( 'position = ', self.x)
    print( 'velocity = ', self.v)

  def where(self):
    return self.x

  def step(self,accel,dt):
  # one leapfrog jump
    d = self.v + accel(self.x)*0.5*dt
    self.x = self.x + d*dt
    self.v = d + accel(self.x)*0.5*dt

  def angularmomentum(self):
    return np.cross(self.x,self.m*self.v)

  def energy(self,potential):
    return 0.5*self.m*sum(self.v**2) + self.m*potential(self.x)
    
  def rungelenz(self):
    return np.cross(self.m*self.v,self.angularmomentum())-gfac*(self.m**2)*self.x/np.sqrt(sum(self.x**2))

  def trajectory(self,accel,dt,n):
    a = np.empty([3,n+1])
    a[:,0]=self.x
    for i in range(1,n+1):
      self.step(accel,dt)
      a[:,i]=self.x
    return a

def plot_trajectory(a):
  plt.ion()
  plt.axis('equal')
  plt.plot(a[0,:],a[1,:],'r')
  plt.plot([0.0],[0.0],'bo')

a     = particle(0.1,np.array([1.0,0.0,0.0]),np.array([0.0,2.0*np.pi*1.0,0.0]))
print()
print('Start:')
a.what()
print( 'Angular momentum = ', a.angularmomentum())
print( 'Total energy     = ', a.energy(potential))
print( 'Runge Lenz       = ', a.rungelenz())
dt    = 0.0002
#nstep = 2800
nstep = (np.ceil(1.0/dt)).astype(np.int64) * 1
orbit = a.trajectory(accel,dt,nstep)
plot_trajectory(orbit)
print
print( 'End:')
a.what()
print( 'Angular momentum = ', a.angularmomentum())
print( 'Total energy     = ', a.energy(potential))
print( 'Runge Lenz       = ', a.rungelenz())
orbit = a.trajectory(accel,-dt,nstep)
plot_trajectory(orbit)
print
print( 'Reversal End:')
a.what()
print( 'Angular momentum = ', a.angularmomentum())
print( 'Total energy     = ', a.energy(potential))
print( 'Runge Lenz       = ', a.rungelenz())

# save
plt.savefig('./trajectory.png', bbox_inches="tight")