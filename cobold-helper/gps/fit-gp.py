# %%
from scipy.io import readsav
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
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
X_train = np.vstack((teff, logg, moh)).T  # shape (n_samples, n_features)

X_train, X_pred = X_train[:25], X_train[25:]
h1_train, h1_pred = h1[:25], h1[25:]

print(X_train.shape)
# %%
# Plot the data
fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharey=False)

# Effective temperature
axes[0][0].plot(teff, h1, '.')
axes[0][0].set_xlabel(r"T$_\mathrm{eff}$")
axes[0][0].set_ylabel(r"h1")

# Surface gravity
axes[0][1].plot(logg, h1, '.')
axes[0][1].set_xlabel(r"$\log{g}$")
axes[0][0].set_ylabel(r"h1")

# Metallicity
axes[1][0].plot(moh, h1, '.')
axes[1][0].set_xlabel("[M/H]")
axes[0][0].set_ylabel(r"h1")

# h1 vs h2
axes[1][1].errorbar(data.h1[12:44], data.h2[12:44],
                    yerr=data.sh2[12:44], xerr=data.sh1[12:44], fmt='x')

axes[1][1].set_xlabel("h1")
axes[1][1].set_ylabel("h2")
# %%
# Gaussian process regression
kernel = RBF(1.0, (1e-2, 1e-2))
gp_regressor = GaussianProcessRegressor(
    kernel=kernel, n_restarts_optimizer=9)
gp_regressor.fit(X_train, h1_train)

y_pred, sigma = gp_regressor.predict(X_pred, return_std=True)

# Plot fit
fig, axes = plt.subplots(2, 2, sharey=True)

axes[0][0].plot(h1_train, '.')
axes[0][0].plot(gp_regressor.predict(X_train))
axes[0][0].set_xlabel("Training instance")
axes[0][0].set_ylabel("h1")

# vs temperature
x = teff[25:]
axes[0][1].plot(x, h1_pred, '.')
axes[0][1].plot(x, y_pred)
axes[0][1].set_xlabel("Temperature")
# axes[0][1].fill(np.concatenate([x, x[::-1]]),
#                 np.concatenate([y_pred - 1.9600 * sigma,
#                                 (y_pred + 1.9600 * sigma)[::-1]]),
#                 alpha=.5, fc='b', ec='None', label='95% confidence interval')

# vs logg
x = logg[25:]
axes[1][0].plot(x, h1_pred, '.')
axes[1][0].plot(x, y_pred)
axes[1][0].set_xlabel("Surface gravity")
# axes[1][0].fill(np.concatenate([x, x[::-1]]),
#                 np.concatenate([y_pred - 1.9600 * sigma,
#                                 (y_pred + 1.9600 * sigma)[::-1]]),
#                 alpha=.5, fc='b', ec='None', label='95% confidence interval')

# vs metallicity
x = moh[25:]
axes[1][1].plot(x, h1_pred, '.')
axes[1][1].plot(x, y_pred)
axes[1][1].set_xlabel("[M/H]")
# axes[1][1].fill(np.concatenate([x, x[::-1]]),
#                 np.concatenate([y_pred - 1.9600 * sigma,
#                                 (y_pred + 1.9600 * sigma)[::-1]]),
#                 alpha=.5, fc='b', ec='None', label='95% confidence interval')
# %%
