import numpy as np
import matplotlib.pyplot as plt


def main():
	# Create some fake data
	x = np.linspace(-10., 10., num=100)
	y1 = np.sin(x)
	y2 = np.cos(x)

	# Single plot
	fig_single, ax_single = plt.subplots(nrows=1, ncols=1)  # my preferred way of plotting! Also see 'plt.figure(); plt.plot();'
	ax_single.plot(x, y1, ls='-', c='r', label="sin")
	ax_single.plot(x, y2, ls='--', c='#587427', label="cos")

	ax_single.axhline(ls=':', c='k')  # horizontal line, without args, defaults at 'y=0'
	for point in np.arange(-3*np.pi, 3*np.pi, step=np.pi):
		# see https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
		# for more info on linestyles
		ax_single.axvline(point, ls=(0, (1, 10)), c='k')

	ax_single.legend()  # populate legend in 'ax_single' with 'label's from above
	ax_single.set_xlabel("x")
	# Use 'r' (raw) before strings to also do LaTeX magic
	ax_single.set_ylabel(r"y$_x$")

	# fig.show()  # doesn't work for me on this PC, but uncomment to check
	# plt.show() # Alternative display that shows all open figures
	fig_single.savefig('./single.png', bbox_inches="tight")
	# fig.close()  # close the figure

	# 2x2 subplots with random data
	nrows, ncols = 2, 2
	# subplots() returns fig and:
	# - a numpy array when dims > 1x2 or 2x1 
	# - a list of subplots when one of the dims is size 1
	# - AxesSubplot when nrows=ncols=1 (as above)
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
	colours = ['m', 'g']
	for i in range(nrows):
		for j in range(ncols):
			if j == 0:
				axes[i, j].plot(np.random.random(size=1000), marker='o', ls='none',
																										c=colours[i])
			else:
				# histogram
				axes[i, j].hist(np.random.normal(size=100*(i+1)))

	fig.savefig('./2x2.png', bbox_inches="tight")

	"""
	So much for the basics, you can also check out cool things like
	pcolormesh & imshow (heatmaps), contourf (contours), 3D plots and
	quiver (plots arrows on top for vector fields).
	"""

# This little snippet means if you call the file from the command line
# as the 'main' file (not an import), it runs the following
if __name__ == "__main__":
	# main() is a separate function just for good practice sake
	main()