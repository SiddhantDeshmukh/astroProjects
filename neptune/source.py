import numpy
from matplotlib.pyplot import *
import math

i_x = 0
i_y = 1
i_z = 2
i_u = 3
i_v = 4
i_w = 5

rk1      = 0
rk2      = 1
leapfrog = 2

class particles:

   def __init__(self,N,time_integrator,eps,cfl):


         self.m     = numpy.zeros((N))

         self.coord = numpy.zeros((N,3))
         self.vel   = numpy.zeros((N,3))

         self.coord_n   = numpy.zeros((N,3))
         self.vel_n     = numpy.zeros((N,3))
         
         self.R     = numpy.zeros((N,6))
         self.a     = numpy.zeros((N,3))
         self.F     = numpy.zeros((N,3))

         self.dumb  = numpy.zeros((N))

         self.eps   = eps

         self.N     = N
         self.cfl   = cfl
         self.dt    = 0.

         self.V     = 0.
         self.K     = 0.

         self.time_marching_scheme = time_integrator


   def do_timestep(self):

         if self.time_marching_scheme == "rk1":

              self.rk1()

         elif self.time_marching_scheme == "rk2":

              self.rk2()

         elif self.time_marching_scheme == "leapfrog":

              self.leapfrog()

           
   def get_acceleration(self):

            for i in range(self.N):

                      dx = self.coord[i,i_x]-self.coord[:,i_x]
                      dy = self.coord[i,i_y]-self.coord[:,i_y]
                      dz = self.coord[i,i_z]-self.coord[:,i_z]

                      r2 = (dx**2+dy**2+dz**2)
                      r     = r2**0.5
                      
                      self.F[:,:] = 0
                      self.F[r2>0,i_x] = -self.m[i]*self.m[r2>0]/(r[r2>0]**3+self.eps**3)*dx[r2>0]
                      self.F[r2>0,i_y] = -self.m[i]*self.m[r2>0]/(r[r2>0]**3+self.eps**3)*dy[r2>0]
                      self.F[r2>0,i_z] = -self.m[i]*self.m[r2>0]/(r[r2>0]**3+self.eps**3)*dz[r2>0]
                      self.a[i,:] = numpy.sum(self.F,axis=0)/self.m[i]  
  

   def get_residuals(self):


          self.R[:,i_x:i_z+1] = self.vel[:,:]*self.dt
          self.R[:,i_u:i_w+1] = self.a[:,:]*self.dt


   def rk1(self):

          #first stage
          self.get_acceleration()
          self.get_residuals()

          self.coord[:,:] +=  self.R[:,i_x:i_z+1] 
          self.vel[:,:]   +=  self.R[:,i_u:i_w+1]


   def rk2(self):


          self.coord_n[:,:]   = self.coord[:,:]
          self.vel_n[:,:]     = self.vel[:,:]


          #first stage
          self.get_acceleration()
          self.get_residuals()

          self.coord[:,:] = self.coord_n[:,:] + self.R[:,i_x:i_z+1]/2.
          self.vel[:,:]   = self.vel_n[:,:]   + self.R[:,i_u:i_w+1]/2.



          #second stage
          self.get_acceleration()
          self.get_residuals()

          self.coord[:,:] = self.coord_n[:,:] + self.R[:,i_x:i_z+1]
          self.vel[:,:]   = self.vel_n[:,:]   + self.R[:,i_u:i_w+1]


   def leapfrog(self):

          self.coord_n[:,:]   = self.coord[:,:]
          self.vel_n[:,:]     = self.vel[:,:]
          #first stage
          self.vel[:,:]   = self.vel_n[:,:] + self.a[:,:]*self.dt/2.
          #second stage
          self.coord[:,:] = self.coord_n[:,:] + self.vel[:,:]*self.dt
          #third stage
          self.get_acceleration()
          self.vel[:,:] = self.vel[:,:] + self.a[:,:]*self.dt/2.



   def get_dt(self):

            rmin = 1e+33

            for i in range(self.N):

                      dx = self.coord[i,i_x]-self.coord[:,i_x]
                      dy = self.coord[i,i_y]-self.coord[:,i_y]
                      dz = self.coord[i,i_z]-self.coord[:,i_z]

                      r2 = (dx**2+dy**2+dz**2)
                      r     = r2**0.5
                      
                    
                      rmin = min(numpy.min(r[r>0]),rmin)

           
            vmax = numpy.max((self.vel[:,i_x]**2+self.vel[:,i_y]**2+self.vel[:,i_z]**2)**0.5)
            self.dt = self.cfl*rmin/vmax



   def get_V(self):

           self.V = 0.
           for i in range(self.N):

                      dx = self.coord[i,i_x]-self.coord[:,i_x]
                      dy = self.coord[i,i_y]-self.coord[:,i_y]
                      dz = self.coord[i,i_z]-self.coord[:,i_z]

                      r2 = (dx**2+dy**2+dz**2)
                      r     = r2**0.5
                      
                      self.dumb[:] = 0.
                      self.dumb[r2>0]  = -self.m[i]*self.m[r2>0]/r[r2>0]

                      self.V += numpy.sum(self.dumb)

           self.V = self.V/2.


   def get_K(self):

           self.K = numpy.sum(0.5*self.m*(self.vel[:,i_x]**2+self.vel[:,i_y]**2+self.vel[:,i_z]**2))

   def get_Etot(self):
  
           self.get_V()
           self.get_K()
           self.E = self.V + self.K
           return self.E

           

