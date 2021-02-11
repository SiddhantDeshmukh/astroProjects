from source import *


#input parameters
N               = 2
eps             = 1.e-4
cfl             = 0.01
# time_integrator = leapfrog #rk1, rk2
time_integrators = ["leapfrog", "rk1", "rk2"]
tmax            = 12.*math.pi

#create particles class
for integrator in time_integrators:
      p = particles(N,integrator,eps,cfl)

      #setup
      p.m[0] = 1.
      p.m[1] = 1.e-6

      p.coord[0,i_x] = 0.
      p.coord[0,i_y] = 0.

      p.coord[1,i_x] = 1.
      p.coord[1,i_y] = 0.

      p.vel[0,i_x] = 0.
      p.vel[0,i_y] = 0.

      p.vel[1,i_x] = 0.
      p.vel[1,i_y] = 1.

      if p.time_marching_scheme == "leapfrog":
            p.get_acceleration()

      figure()
      #main loop
      t    = 0.
      E    = []
      time = []
      num_iter = 0

      while(t<tmax):
            time.append(t)
            E.append(p.get_Etot())
            p.get_dt()
            p.do_timestep()
            scatter(p.coord[:,i_x],p.coord[:,i_y], color='blue', s=1)
            t+=p.dt
            num_iter += 1

            if num_iter % 100 == 0:
                  print(num_iter, t, tmax)


      xlabel('x [A.U.]')
      ylabel('y [A.U.]')
      title('2-body')
      savefig(f'./figs-2body/2body_long_{integrator}.png')


      figure()
      plot(time,numpy.log10(numpy.abs((E-E[0])/E[0])))
      xlabel('t')
      ylabel(r'log|$\frac{E-E_0}{E_0}$|')
      savefig(f'./figs-2body/energy_2body_long_{integrator}.png')