# %%
from scipy.io import readsav
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, ExpSineSquared, WhiteKernel
import numpy as np

# Interested in
# mid: underlying 3D model ID
# Teff, logg, moh: as perhaps clear moh=[M/H]
# h1, h2: I guess clear
# sh1, sh2: 1sigma uncertainty on h1, and h2
# ch1h2: linear correlation coefficient between h1 and h2
# Use data.dtype to see the possible records
data = readsav('./shh.idlsave', verbose=True)['shh'][0]
print(data.dtype)

# We only care about models in the range 12-44 (zero-indexed)
h1, h2, sh1, sh2 =\
    data.h1[12:44], data.h2[12:44], data.sh1[12:44], data.sh2[12:44]
teff, logg, moh = data.teff[12:44], data.logg[12:44], data.moh[12:44]

# Creating training data
# Fit 'h1' and 'h2' separately using 'teff', 'logg', 'moh'
X = np.vstack((teff, logg, moh)).T  # shape (n_samples, n_features)
h1_kernel = RBF(length_scale=1) + WhiteKernel(1e-1)
gp_h1 = GaussianProcessRegressor(kernel=h1_kernel)
gp_h1.fit(X, h1)
h1_pred, h1_std = gp_h1.predict(X, return_std=True)

h2_kernel = RBF(length_scale=1) + WhiteKernel(1e-1)
gp_h2 = GaussianProcessRegressor(kernel=h2_kernel)
gp_h2.fit(X, h2)
h2_pred, h2_std = gp_h2.predict(X, return_std=True)

print(h1_std, h2_std)

fig, axes = plt.subplots(1, 1)
axes.errorbar(h1, h2, xerr=sh1, yerr=sh2, fmt='k.', label="Data")
axes.plot(h1_pred, h2_pred, 'r-', label="GP fits")
axes.fill_between(h1, h2_pred - h2_std, h2_pred +
                  h2_std, color='darkorange')

# %%
# Let's just try a random example
# Sample data with noise
rng = np.random.RandomState(42)
X = 15 * rng.rand(100, 1)
y1 = np.sin(X).ravel()
y2 = np.cos(X).ravel()

y1 += 3 * (0.5 - rng.rand(X.shape[0]))
y2 += 2 * (0.1 - rng.rand(X.shape[0]))

targets = np.vstack((y1, y2)).T
print(targets.shape)

# GP Regressor
gp_kernel = ExpSineSquared(1., 5.0, periodicity_bounds=(1e-2, 1e1)) +  \
    WhiteKernel(1e-1)
gpr = GaussianProcessRegressor(kernel=gp_kernel)
gpr.fit(X, targets)
X_plot = np.linspace(0, 15, 1000).reshape(-1, 1)
X = np.sort(X, axis=0)
h1_pred, y_std = gpr.predict(X, return_std=True)

fig, axes = plt.subplots(1, 2)
axes[0].plot(X, y1, 'k.', label="Data")
axes[0].plot(X_plot, np.sin(X_plot), 'b-', label="True")
axes[0].plot(X, h1_pred[:, 0], 'r-', label="GP Fit")
axes[0].fill_between(X[:, 0], h1_pred[:, 0] - y_std,
                     h1_pred[:, 0] + y_std, color='darkorange', alpha=0.2)

axes[0].set_title("sin")
axes[0].legend()

axes[1].set_title("cos")
axes[1].plot(X, y1, 'k.', label="Data")
axes[1].plot(X_plot, np.cos(X_plot), 'b-', label="True")
axes[1].plot(X, h1_pred[:, 1], 'r-', label="GP Fit")
axes[1].fill_between(X[:, 0], h1_pred[:, 1] - y_std,
                     h1_pred[:, 1] + y_std, color='darkorange', alpha=0.2)

axes[1].legend()

fig = plt.figure()
ax2 = fig.add_subplot(111, projection='3d')
ax2.scatter(X, y1, y2, 'k.')
ax2.plot(X_plot[:, 0], np.sin(X_plot[:, 0]), np.cos(X_plot[:, 0]), 'b-')
ax2.plot(X[:, 0], h1_pred[:, 0], h1_pred[:, 1], 'r-')
