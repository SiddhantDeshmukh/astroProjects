# Credit: https://github.com/pmocz/sph-python/blob/master/sph.py
# =========================================================================
# Imports
# =========================================================================#
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

# =========================================================================
# Classes
# =========================================================================
class Parameters():
  # Defines parameters for simulation, such as potential and viscosity
  # terms, number of particles and smoothing length
  def __init__(self, mass: float, smoothing_length: float, gamma: float,
               potential: float, viscosity: float):
    self.m = mass
    self.h = smoothing_length
    self.gamma = gamma
    self.pot = potential
    self.nu = viscosity


class Simulation():
  # Parent class API for simulation setups
  def __init__(self, N: int, t_start: float, t_end: float, dt: float) -> None:
    self.N = N
    self.t = t_start
    self.t_end = t_end
    self.dt = dt

# =========================================================================
# Kernel functions and derivatives
# =========================================================================
def gaussian_kernel (x: np.ndarray, y: np.ndarray, z: np.ndarray, h: float):
  # Gaussian smoothing kernel
  r = np.sqrt(x**2+ y**2 + z**2)
  w = (1.0 / (h * np.sqrt(np.pi)))**3 * np.exp(-(r / h)**2)

  return w

def gaussian_derivative(x: np.ndarray, y: np.ndarray, z: np.ndarray, h: float):
  # Gradient of gaussian_kernel()
  r = np.sqrt(x**2 + y**2 + z**2)
  n = -2 * np.exp(-(r / h)**2) / h**5 / np.pi**(3/2)

  # Derivative components
  wx = n * x
  wy = n * y
  wz = n * z

  return wx, wy, wz


def cubic_spline_kernel(x: np.ndarray, y: np.ndarray, z: np.ndarray, h: float):
  # Cubic spline smoothing kernel (Schoenberg B-spline M_4)
  pass

# =========================================================================
# Density and pressure calculation
# =========================================================================
def get_pairwise_separations(ri: np.ndarray, rj: np.ndarray):
  # Calculate pairwise separations 'r_i - r_j'
  M = ri.shape[0]
  N = rj.shape[0]

  # Get x, y, z positions from ri and rj
  rix = ri[:,0].reshape((M,1))
  riy = ri[:,1].reshape((M,1))
  riz = ri[:,2].reshape((M,1))
  
  rjx = rj[:,0].reshape((N,1))
  rjy = rj[:,1].reshape((N,1))
  rjz = rj[:,2].reshape((N,1))
  
  # matrices that store all pairwise particle separations: r_i - r_j
  dx = rix - rjx.T
  dy = riy - rjy.T
  dz = riz - rjz.T
  
  return dx, dy, dz

def get_density(r: np.ndarray, pos: np.ndarray, m: float, h: float):
  # 'r' is sampling locations, 'pos' is SPH particle positions
  M = r.shape[0]
  dx, dy, dz = get_pairwise_separations(r, pos)
  rho = np.sum(m * gaussian_kernel(dx, dy, dz, h), 1).reshape((M, 1))

  return rho

# Should be a method of the 'Simulation()'
def get_pressure(rho: np.ndarray, vel: np.ndarray, gamma=5/3,
                 k=None, n=None):
  if k is not None and n is not None: # Polytropic EOS
    P  = k * rho**(1 + 1 / n)

  else:  # ideal gas EO
    P = (gamma - 1) * rho * vel

  return P

# =========================================================================
# Acceleration calculations
# =========================================================================
def get_acceleration(pos: np.ndarray, vel: np.ndarray, m: float, h: float,
                     k: float, n: float, potential: float, nu: float):
  N = pos.shape[0]

  # Calculate densities at position of particles
  rho = get_density(pos, pos, m, h)

  # Calculate pressures
  P = get_pressure(rho, vel, k=k, n=n)

  # Get pairwise distances and gradients
  dx, dy, dz = get_pairwise_separations(pos, pos)
  dWx, dWy, dWz = gaussian_derivative(dx, dy, dz, h)

  # Add pressure contribution to accelerations
  ax = -np.sum(m * (P / rho**2 + P.T / rho.T**2) * dWx, 1).reshape((N, 1))
  ay = -np.sum(m * (P / rho**2 + P.T / rho.T**2) * dWy, 1).reshape((N, 1))
  az = -np.sum(m * (P / rho**2 + P.T / rho.T**2) * dWz, 1).reshape((N, 1))

  a = np.hstack((ax, ay, az))

  # Add external potential force and viscosity
  a += -potential * pos - nu * vel

  return a

# =========================================================================
# Energy calculations
# =========================================================================
def get_internal_energy(pressure: np.ndarray, rho: np.ndarray, mass: float, n=1):
  gamma_idx = 1 + 1 / n
  energy = np.mean(pressure / (rho * mass)) / (gamma_idx - 1)

  return energy

# =========================================================================
# Time integration schemes
# =========================================================================
def leapfrog_step(vel: np.ndarray, pos: np.ndarray, acc: np.ndarray,
                  dt: float, pars: Parameters):
  # half-timestep kick
  vel += acc * dt / 2

  # drift
  pos += vel * dt

  # update accelerations
  acc = get_acceleration(pos, vel, pars.m, pars.h, pars.gamma, pars.pot, pars.nu)

  # half-timestep kick
  vel += acc * dt / 2

  return pos, vel, acc

# =========================================================================
# Main loop
# =========================================================================
if __name__ == "__main__":
  # Simulation parameters for star formation
  N = 400  # number of particles
  t = 0  # current sim time
  t_end = 12 # final sim time
  dt = 0.04  # timestep
  M = 2  # total sim mass (all stars)
  R = 0.75  # sim radius start (cloud radius)
  h = 0.1  # smoothing length
  k = 0.1  # EOS constant
  n = 1  # polytropic index
  nu = 1  # damping

  # Initial conditions
  np.random.seed(42)

  potential = 2 * k * (1 + n) * np.pi**(-3/(2*n)) * (M * gamma(5/2 + n) / R**3 / gamma (1 + n))** (1 / n)  / R**2  # ~2.01
  m = M / N  # single star mass
  pos = np.random.randn(N, 3)  # random positions and velocities
  vel = np.zeros(pos.shape)
  
  # Initial accelerations, densities, pressures
  acc = get_acceleration(pos, vel, m, h, k, n, potential, nu)
  rho = get_density(pos, pos, m, h)
  pressure = get_pressure(rho, vel, k=k, n=n)
  initial_energy = get_internal_energy(pressure, rho, m, n=n)

  # Quantity arrays for evolution plots
  times = []
  energies = []

  # Set up plot
  fig = plt.figure(figsize=(4, 5), dpi=80)
  grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
  ax1 = plt.subplot(grid[0:2, 0])
  ax2 = plt.subplot(grid[2, 0])

  rr = np.zeros((100, 3))
  rlin = np.linspace(0, 1, 100)
  rr[:, 0] = rlin
  rho_analytic = potential / (4*k) * (R**2 - rlin**2)

  # Evolve through time
  num_timesteps = int(np.ceil(t_end / dt))
  for i in range(num_timesteps):
    # (1/2) kick
    vel += acc * dt/2
    
    # drift
    pos += vel * dt
    
    # update accelerations
    acc = get_acceleration(pos, vel, m, h, k, n, potential, nu)
    
    # (1/2) kick
    vel += acc * dt/2
    
    # update time
    t += dt
    
    # get density for plotting
    rho = get_density(pos, pos, m, h)

    # pressure and energy
    pressure = get_pressure(rho, vel, k=k, n=n)
    energy = get_internal_energy(pressure, rho, m, n=n)

    times.append(t)
    energies.append(energy)

    if i % 1 == 0:
      if i % 100 == 0:
        print(f"Done {i+1} / {num_timesteps}")

      # Plot and save
      # Positions in real-time with colour as density
      plt.sca(ax1)
      plt.cla()
      cval = np.minimum((rho - 3) / 3, 1).flatten()
      ax1.scatter(pos[:, 0], pos[:, 1], c=cval, cmap=plt.cm.autumn, s=10, alpha=0.5)
      ax1.set(xlim=(-1.4, 1.4), ylim=(-1.2, 1.2))
      ax1.set_aspect('equal', 'box')
      ax1.set_xticks([-1, 0, 1])
      ax1.set_yticks([-1, 0, 1])
      ax1.set_facecolor('black')
      ax1.set_facecolor((0.1, 0.1, 0.1))

      # Density plot vs analytic
      plt.sca(ax2)
      plt.cla()
      ax2.set(xlim=(0, 1), ylim=(0, 3))
      # ax2.set_aspect(0.1)
      ax2.plot(rlin, rho_analytic, color='grey', linewidth=2)
      rho_radial = get_density(rr, pos, m, h)
      ax2.plot(rlin, rho_radial, color='blue')

      # Labels / legend / title
      ax2.set_xlabel("Radius")
      ax2.set_ylabel("Density")

      plt.suptitle(f"N = {N}, t = {t:0.2f}")

      plt.savefig(f'./figs/frames-star/sph-{str(i).zfill(3)}.png', dpi=240, bbox_inches="tight")

  # Check energy conservation (plot)
  fig, axes = plt.subplots(figsize=(4, 4))
  axes.plot(times, energies)
  axes.axhline(initial_energy, color='k')
  axes.set_xlabel("Time")
  axes.set_ylabel("Internal Energy")
  plt.savefig('./figs/energy.png', bbox_inches="tight")
