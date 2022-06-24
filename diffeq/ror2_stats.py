from typing import Union
import numpy as np
import matplotlib.pyplot as plt


def tougher_times(num_stacks: Union[int, np.ndarray]):
	# hyperbolic stacking
	alpha = 0.15  # block chance
	a = alpha * num_stacks
	xn = a / (a + 1)

	return xn

def safer_spaces(num_stacks: Union[int, np.ndarray]):
	# geometric?
	x0 = 15.  # initial cooldown
	alpha = 0.1  # per-stack reduction
	# x1 = x0 - (x0 * a) = x0 * (1 - a)
	# x2 = x1 - (x1 * a) = x0 * (1 - a)^2
	# ...
	# xn = x0 * (1 - a)^n; note that here 'n' is num stacks!
	xn = x0 * (1 - alpha)**num_stacks

	return xn

# Simple diffeq for ror2
def main():
	num_stacks = np.array(range(0, 26))
	safer_cooldowns = safer_spaces(num_stacks)
	tougher_chance = tougher_times(num_stacks)

	fig, axes = plt.subplots(2, 1)

	axes[0].plot(num_stacks, tougher_chance, label="Tougher Times",
					marker='o', mfc='none', c="orange")
	axes[0].plot(num_stacks, 1 / safer_cooldowns, label="Safer Spaces",
					marker='o', mfc='none', c="purple")
	axes[1].plot(num_stacks, safer_cooldowns, marker='o', mfc='none',
								c="purple")
	axes[0].legend()
	
	axes[1].set_xlabel("Num Stacks")
	axes[0].set_ylabel("Hits avoided per second")
	axes[1].set_ylabel("Cooldown")
	plt.savefig("bears.png", bbox_inches="tight")

if __name__ == "__main__":
	main()