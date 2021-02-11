# An N-body system integrator using the leapfrog method
# =========================================================================
# Imports
# =========================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# =========================================================================
# Classes
# =========================================================================
G = 6.67E-11  # SI


class NbodySystem:
  def __init__(
      self,
      input_file=None,
      N=1,
      names=None,
      masses="random",
      positions="random",
      velocities="random",
      softening_length=1e-4
  ) -> None:
    self.N = None
    self.names = []
    self.masses = None
    self.positions = None
    self.velocities = None

    self.eps = softening_length
    self.dt_max = 10000
    self.dt_min = 1

    # Initialise positions, velocities and accelerations of system
    if input_file:  # read nbody details from a text file
      df = pd.read_csv(input_file, skipinitialspace=True, comment='#')

      self.N = len(df)
      self.masses = np.zeros((self.N))
      self.positions = np.zeros((self.N, 3))
      self.velocities = np.zeros((self.N, 3))

      for i, row in enumerate(df.itertuples(index=False)):
        self.names.append(row.body)
        self.masses[i] = row.mass
        self.positions[i] = np.array([row.x, row.y, row.z])
        self.velocities[i] = np.array([row.vx, row.vy, row.vz])

    else:
      self.N = N
      for i in range(N):
        if names is None:
          self.names.append(f"body_{i+1}")

        if masses == "random":
          self.masses[i] = np.random.random()

        if positions == "random":
          x = np.random.normal()
          y = np.random.normal()
          z = np.random.normal()

          self.positions[i] = np.array([x, y, z])

        if velocities == "random":
          vx = np.random.normal()
          vy = np.random.normal()
          vz = np.random.normal()

          self.velocities[i] = np.array([[vx, vy, vz]])

    self.accelerations = np.zeros((self.N, 3))

    self.update_accelerations()
    self.update_timestep()

    # Histories
    self.position_history = []
    self.velocity_history = []

  def update_timestep(self):
    min_distance = 1e+33
    for i in range(self.N):
      dx = self.positions[i, 0] - self.positions[:, 0]
      dy = self.positions[i, 1] - self.positions[:, 1]
      dz = self.positions[i, 2] - self.positions[:, 2]

      r2 = (dx**2 + dy**2 + dz**2)
      r = r2**0.5

      min_distance = min(np.min(r[r > 0]), min_distance)

    max_velocity = np.max(
        (self.velocities[:, 0]**2 + self.velocities[:, 1]**2 + self.velocities[:, 2])**0.5)

    dt = min_distance / max_velocity * 0.01  # CFL = 0.01
    if dt > self.dt_max:
      dt = self.dt_max

    # if dt < self.dt_min:
    #   dt = self.dt_min

    self.dt = dt

  def update_accelerations(self):
    self.accelerations = np.zeros((self.N, 3))
    # iterate over bodies
    for i in range(self.N):
      dx = self.positions[i, 0] - self.positions[:, 0]
      dy = self.positions[i, 1] - self.positions[:, 1]
      dz = self.positions[i, 2] - self.positions[:, 2]

      r2 = (dx**2 + dy**2 + dz**2)
      r = r2**0.5

      forces = np.zeros((self.N, 3))

      forces[r2 > 0, 0] = -G * self.masses[i] * self.masses[r2 > 0] / \
          (r[r2 > 0]**3 + self.eps**3) * dx[r2 > 0]
      forces[r2 > 0, 1] = -G * self.masses[i] * self.masses[r2 > 0] / \
          (r[r2 > 0]**3 + self.eps**3) * dy[r2 > 0]
      forces[r2 > 0, 2] = -G * self.masses[i] * self.masses[r2 > 0] / \
          (r[r2 > 0]**3 + self.eps**3) * dz[r2 > 0]

      self.accelerations[i, :] = np.sum(forces, axis=0) / self.masses[i]

  def leapfrog_step(self):
    # Kick half-timestep
    self.velocities += self.accelerations * self.dt / 2

    # Drift full timestep
    self.positions += self.velocities * self.dt

    # Kick half-timestep
    self.update_accelerations()
    self.velocities += self.accelerations * self.dt / 2

  def evolve(self, initial_time: float, final_time: float, ax):
    current_time = initial_time
    num_steps = 0
    max_steps = 10000
    while current_time <= final_time and num_steps < max_steps:
      if num_steps % 100 == 0:
        print(num_steps, current_time)

      self.position_history.append(self.positions.copy())
      self.velocity_history.append(self.velocities.copy())

      self.leapfrog_step()
      current_time += self.dt
      num_steps += 1
      self.update_timestep()

    self.plot_trajectory(ax)

  def plot_trajectory(self, ax: plt.axes):  # x-y plane projection
    for i in range(self.N):
      ax.scatter(np.array(self.position_history)[:, i, 0],
                 np.array(self.position_history)[:, i, 1], s=1,
                 label=self.names[i])

    ax.legend()


# =========================================================================
# ESM system
# =========================================================================
if __name__ == "__main__":
  esm_path = "./res/solar-system-input.csv"
  ESM_System = NbodySystem(input_file=esm_path)

  fig, axes = plt.subplots(1, 1, figsize=(4, 4))
  ESM_System.evolve(0., 1e10, axes)

  plt.savefig('./solar-system.png', bbox_inches="tight")
