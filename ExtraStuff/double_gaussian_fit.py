import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.io import fits

# ----------------------------
# 1. Define a double Gaussian
# ----------------------------
def two_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3):
    g1 = A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2))
    g2 = A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2))
    g3 = A3 * np.exp(-(x - mu3)**2 / (2 * sigma3**2))
    return g1 + g2 + g3

# ---------------------------------------------------
# 2. Insert your data here (replace x_data, y_data)
# ---------------------------------------------------

hdul = fits.open(f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/1561/repro/01561_sgra_2-8keV_lc300_pileup.fits')
data = hdul[1].data

#flare_evts_mask = ((50814 + data['time']/86400) > (58715.156961405)) & ((50814 + data['time']/86400) < (58715.2262554391))
flare_evts_mask = ((50814 + data['time']/86400) > (51844.1581720461)) & ((50814 + data['time']/86400) < (51844.2795123523))

x_data = data['TIME'][flare_evts_mask]   # fill in
y_data = data['RATE_PILEUP'][flare_evts_mask]   # fill in

quiesc_idx = np.argwhere(y_data < 3*0.005)
print(len(y_data))
print(quiesc_idx)

# ---------------------------------------------------
# 3. Initial guesses — IMPORTANT for convergence
# ---------------------------------------------------
A1_guess = np.max(y_data)
A2_guess = np.max(y_data) * 0.5
A3_guess = np.max(y_data) * 0.1
mu1_guess = x_data[np.argmax(y_data)]
mu2_guess = mu1_guess + (x_data.max() - x_data.min()) * 0.25
mu3_guess = mu1_guess - (x_data.max() - x_data.min()) * 0.25
sigma1_guess = (x_data.max() - x_data.min()) / 10
sigma2_guess = (x_data.max() - x_data.min()) / 10
sigma3_guess = (x_data.max() - x_data.min()) / 10

p0 = [A1_guess, mu1_guess, sigma1_guess,
      A2_guess, mu2_guess, sigma2_guess,
      A3_guess, mu3_guess, sigma3_guess]

# ---------------------------------------------------
# 4. Fit the model to the data
# ---------------------------------------------------
popt, pcov = curve_fit(two_gaussian, x_data, y_data, p0=p0)

A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = popt
perr = np.sqrt(np.diag(pcov))  # parameter errors

print("Fit results:")
print(f"A1 = {A1:.3f} ± {perr[0]:.3f}")
print(f"mu1 = {mu1:.3f} ± {perr[1]:.3f}")
print(f"sigma1 = {sigma1:.3f} ± {perr[2]:.3f}")
print(f"A2 = {A2:.3f} ± {perr[3]:.3f}")
print(f"mu2 = {mu2:.3f} ± {perr[4]:.3f}")
print(f"sigma2 = {sigma2:.3f} ± {perr[5]:.3f}")
print(f"A3 = {A3:.3f} ± {perr[6]:.3f}")
print(f"mu3 = {mu3:.3f} ± {perr[7]:.3f}")
print(f"sigma3 = {sigma3:.3f} ± {perr[8]:.3f}")

# ---------------------------------------------------
# 5. Plot the result
# ---------------------------------------------------
plt.figure(figsize=(7,5))
plt.scatter(x_data, y_data, color='black', label="Data")

# fitted curve
x_fit = np.linspace(x_data.min(), x_data.max(), 500)
y_fit = two_gaussian(x_fit, *popt)
plt.plot(x_fit, y_fit, 'r-', label="Double Gaussian Fit")

# individual components
plt.plot(x_fit, A1 * np.exp(-(x_fit - mu1)**2 / (2 * sigma1**2)),
         'b--', label="Gaussian 1")
plt.plot(x_fit, A2 * np.exp(-(x_fit - mu2)**2 / (2 * sigma2**2)),
         'g--', label="Gaussian 2")
plt.plot(x_fit, A3 * np.exp(-(x_fit - mu3)**2 / (2 * sigma3**2)),
         'g--', label="Gaussian 2")

plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.show()