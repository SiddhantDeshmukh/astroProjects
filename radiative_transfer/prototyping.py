# %%
# # general workspace for radiative transfer solution
# let's try a monte-carlo method
# We want to see how an ensemble of N photons behaves in terms of
# statistical properties
# 1. emit N photon packets (hereafter 'photons')
# 2. Track progress of photons one-by-one through medium. Locations of
#	 interactions found by sampling optical depth from distribution
#	 P(tau) = 1 - e^(-tau)
#	 Scattering and absorption determined from sampling albedo and phase
# 	 functions
# 3. As photons exit medium, capture on a pixelated image plane (like how
#  	 real photons are captured on a CCD)
import numpy as np
import matplotlib.pyplot as plt
from itertools import product


k_b = 1.3806503e-23  # [m^2 kg s^-2 K^-1]
sigma = 5.670367e-8  # [W m^-2 K^-4]


def sample_tau():
	# exact inversion of probability P(tau) = 1 - e^(-tau)
	return -np.log(1 - np.random.random())


def luminosity(radius: float, effective_temperature: float):
	# point source, stefan-boltzmann law
	return 4 * np.pi * radius**2 * sigma * effective_temperature**4


def compute_distance_in_cell(position: np.ndarray, direction: np.ndarray,
                             x_face: np.ndarray, y_face: np.ndarray,
                             z_face: np.ndarray):
	# for a given photon in a grid cell, find the distance to the closest
	# wall along the photon's path
	# 'direction' is a 3-component unit vector (nx, ny, nz)
	# 'x,y,z_face' are the positions of the nearest cell faces
	faces = np.array([x_face[0], y_face[1], z_face[2]])
	return np.min((faces - position) / direction)


def optical_depth_in_cell(density: float, opacity: float, distance: float):
	# for distance travelled along photon's trajectory
	return density * opacity * distance


def create_grid(x_limits=(-5., 5.),
                y_limits=(-5., 5.),
                z_limits=(0., 20.),
                n_x_points=10,
                n_y_points=10,
                n_z_points=20):
	x = np.linspace(*x_limits, num=n_x_points)
	y = np.linspace(*y_limits, num=n_y_points)
	z = np.linspace(*z_limits, num=n_z_points)
	return np.meshgrid(x, y, z, indexing='ij')

def determine_grid_spacing(grid: np.ndarray):
	# assumes uniform grid spacing
	return np.diff(grid)[:, 0]

def cell_from_position(grid: np.ndarray, position: np.ndarray):
	# get coords of cell 'position' is in on 'grid'
	n = grid.shape[:3]
	extents = np.array([grid[-1, -1, -1, i] for i in range(0, 3)])
	cell_idx = tuple(map(int, position * n / extents))

	return cell_idx

def main():
	# Dummy grid to test distance calculations
	grid = np.array(create_grid()).transpose(1, 2, 3, 0)

	# isotropic density grid
	density = np.zeros(shape=(1, *grid.shape[1:]))
	photon_position = np.array([-3., -2., 3.])
	photon_direction = np.array([2., 1., 4.])
	photon_direction /= np.linalg.norm(photon_direction)

	# plot grid and photon
	fig = plt.figure(figsize=(12, 12))
	ax = plt.axes(projection ='3d')
	x, y, z = grid[:, 0, 0, 0], grid[0, :, 0, 1], grid[0, 0, :, 2]
	points = product(*[x, y, z])
	ax.plot3D(*zip(*points), c='gray', marker='.', ls='none')
	# origin
	ax.plot3D([0.], [0.], [5.], 'kx')
	pos = [[p] for p in photon_position]
	ax.plot(*pos, 'ro')
	ax.quiver(*photon_position, *photon_direction, color='r')
	# TODO:
	# Fix index-finding, also figure out how to slice better with the index
	idx = cell_from_position(grid, photon_position)
	photon_grid_pos = grid[idx]
	print(photon_position)
	print(idx)
	print(photon_grid_pos)
	true_idx = (2, 3, 3)
	print(grid[true_idx])
	ax.plot3D([photon_grid_pos[0]], [photon_grid_pos[1]], [photon_grid_pos[2]], 'gs')
	

	plt.show()


if __name__ == "__main__":
	main()

# %%
