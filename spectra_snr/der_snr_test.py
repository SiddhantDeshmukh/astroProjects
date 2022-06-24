# Test the "DER_SNR" algorithm on a dummy spectrum
import numpy as np
import matplotlib.pyplot as plt


def compute_der_snr(flux: np.ndarray, return_components=True):
		n    = len(flux)      

		# For spectra shorter than this, no value can be returned
		if (n>4):
			signal = np.nanmedian(flux)
			components = [2. * flux[2:n-2], flux[0:n-4], flux[4:n]]
			# noise  = 0.6052697 * np.nanmedian(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
			noise  = 1.482602 * np.nanmedian(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
			snr = float(signal / noise)  
		else:
			snr =  0.0
		
		if return_components:
			return snr, components
		else:
			return snr


def load_flux(file_: str):
	arr = np.loadtxt(file_, delimiter=" ").T
	wavelength, flux, error = arr
	return wavelength, flux, error

def main(seed=42):
	np.random.seed(42)
	spectrum_file = "./spec_save_ADP_procyon_snr493_UVES_21.007g.txt"
	target_snr = 493.
	wavelength, flux, error = load_flux(spectrum_file)
	# Values that are exactly zero (padded) are skipped
	idxs = (flux != 0.0)
	wavelength = wavelength[idxs]
	second_arm_start = 5800.  # Angstroms
	mask = (wavelength > second_arm_start)  # only use blue arm
	wavelength = wavelength[mask]
	flux = flux[idxs][mask]
	error = error[idxs][mask]
	num_slices = 50  # number of slices to get an average SNR
	plot_every = 10
	num_plot_slices = int(num_slices / plot_every)  # how many to plot
	mean_snr = 0.
	ax_idx = 0
	fig, axes = plt.subplots(num_plot_slices, 1, figsize=(8, 12))
	for i in range(num_slices):
		# Choose a random slice
		slice_length = np.random.randint(1000, len(flux))
		slice_start = np.random.randint(0, len(flux) - slice_length)
		plot_idxs = slice(slice_start, slice_start + slice_length)
		snr, components = compute_der_snr(flux[plot_idxs],
																			return_components=True)
		mean_snr += snr

		if i % plot_every == 0:
			axes[ax_idx].plot(wavelength[plot_idxs], flux[plot_idxs], marker="o", mfc="none",
								label=f"DER_SNR = {snr:.2f} from {slice_length} points"
								f"\ntarget SNR = {target_snr}")
			axes[ax_idx].legend()
			ax_idx += 1

	mean_snr /= num_slices
	full_snr = compute_der_snr(flux, return_components=False)
	print(f"Average DER_SNR from {num_slices} slices = {mean_snr:.2f}")
	print(f"DER_SNR from entire spectrum = {full_snr:.2f}")
	plt.savefig("./output.png", bbox_inches="tight")
	

if __name__ == "__main__":
	main()
