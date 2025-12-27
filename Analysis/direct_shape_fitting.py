import os
import numpy as np
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit

df = pd.read_csv('flare_properties.csv')
flare_numbers = df['flare_number']
flare_ids = df['obs_id']
max_rate = df['rate_max']
mjd_starts = df['start_mjd']
mjd_ends = df['end_mjd']

def single_peak_model(x, A, width, position):
    y = A*np.exp(-(x-position)**2/(2*width**2))
    return y

def double_peak_model(x, A1, width1, position1, A2, width2, position2):
    if np.abs(position1 - position2) < 500:
        return np.full_like(x, 1e6)
    if np.abs(A1) < 0.005 or np.abs(A2) < 0.005:
        return np.full_like(x, 1e6)
    y = A1*np.exp(-(x-position1)**2/(2*width1**2)) + A2*np.exp(-(x-position2)**2/(2*width2**2))
    return y

for i, flare in enumerate(flare_ids):
    flare = str(flare)
    if flare == '1561':
        lc = f"../1561/1561_000/repro/01561_sgra_2-8keV_lc300_pileup.fits"
    elif len(flare) < 5:
        flare_5digit = str(flare).zfill(5)
        lc = f"../{flare}/repro/{flare_5digit}_sgra_2-8keV_lc300_pileup.fits"
    else:
        lc = f"../{flare}/repro/{flare}_sgra_2-8keV_lc300_pileup.fits"
    f = fits.open(lc)
    data = f[1].data

    chandra_time = data['TIME']
    time = 50814.0 + (chandra_time / 86400.0)
    
    rate = data['RATE_PILEUP']
    rate_err = data['PILEUP_ERR']
    max_rate_reported = max_rate[i]
    
    mjd_start_flare = mjd_starts[i]
    mjd_end_flare = mjd_ends[i]

    # Identify points in the flare window
    flare_mask = (time >= mjd_start_flare) & (time <= mjd_end_flare)

    if flare == str(21454):
        flare_mask[4] = True

    flare_times = time[flare_mask]
    flare_rates = rate[flare_mask]

    
    if len(flare_times) == 0:
        print(f"No lightcurve data inside flare {i+1}")
        continue

    # Find brightest datapoint
    peak_idx = np.argmax(flare_rates)

    

    peak_time = flare_times[peak_idx]
    peak_rate = flare_rates[peak_idx]

    # Extract ±4200s around the peak
    window_mask = (time >= peak_time - 4200/86400) & (time <= peak_time + 4200/86400)
    
    flare_window = rate[window_mask]
    flare_err = rate_err[window_mask]

    #Find the peak idx in the flare time (not of the total observation)
    peak_idx_flare_window = np.where(flare_window == peak_rate)

    #Need to grab 29 total datapoints - 14 on each side. Problem if flare is near the edge of the observation.
    if len(flare_window) != 29:
        if len(peak_idx_flare_window) == 1:
            flare_window = np.concatenate([np.zeros(14 - peak_idx_flare_window[0][0]), flare_window])
            flare_window = np.concatenate([flare_window, np.zeros(29 - len(flare_window))])
            flare_err = np.concatenate([np.zeros(14 - peak_idx_flare_window[0][0]), flare_err])
            flare_err = np.concatenate([flare_err, np.zeros(29 - len(flare_err))])
        else:
            peak_rate_single_flare_idx = np.min(np.abs(14 - peak_idx_flare_window))[0]
            flare_window = np.concatenate([np.zeros(14 - peak_rate_single_flare_idx), flare_window])
            flare_window = np.concatenate([flare_window, np.zeros(29 - len(flare_window))])
            flare_err = np.concatenate([np.zeros(14 - peak_rate_single_flare_idx), flare_err])
            flare_err = np.concatenate([flare_err, np.zeros(29 - len(flare_err))])

    
    #Normalize by peak rate
    flare_window = flare_window/(peak_rate)


    popt, pcov = curve_fit(single_peak_model, np.arange(-4200, 4201, 300), flare_window, sigma=flare_err, absolute_sigma=True, p0=[np.max(flare_window), 600, 0], maxfev=10000)
    model_rate = single_peak_model(np.arange(-4200, 4201, 300), *popt)

    valid = flare_err > 0
    flare_window_valid = flare_window[valid]
    model_rate_valid = model_rate[valid]
    flare_err_valid = flare_err[valid]

    chi2 = np.sum(((flare_window_valid - model_rate_valid) / flare_err_valid) ** 2)
    dof = len(flare_window_valid) - len(popt)
    reduced_chi2 = chi2 / dof




    # Time grid (your model x-values)
    x = np.arange(-4200, 4201, 300)

    # Estimate peak indices
    sorted_peaks = np.argsort(flare_window)[::-1]

    # Take two largest peaks (most likely candidates for double flare)
    peak1_idx = sorted_peaks[0]
    peak2_idx = sorted_peaks[1] if len(sorted_peaks) > 1 else peak1_idx + 3  # fallback if only one clear peak

    # Ensure peak2 is not too close to peak1
    if np.abs(x[peak1_idx] - x[peak2_idx]) < 600:
        peak2_idx = (peak2_idx + 3) % len(x)

    # Extract initial guesses
    A1_guess = flare_window[peak1_idx]
    A2_guess = flare_window[peak2_idx]
    position1_guess = x[peak1_idx]
    position2_guess = x[peak2_idx]

    width1_guess = 600
    width2_guess = 600

    p02 = [A1_guess, width1_guess, position1_guess,
        A2_guess, width2_guess, position2_guess]


    popt2, pcov2 = curve_fit(double_peak_model, np.arange(-4200, 4201, 300), flare_window, sigma=flare_err, absolute_sigma=True, p0=p02, maxfev=10000)
    model_rate2 = double_peak_model(np.arange(-4200, 4201, 300), *popt2)
    print(*popt2)

    valid2 = flare_err > 0
    flare_window_valid2 = flare_window[valid2]
    model_rate_valid2 = model_rate[valid2]
    flare_err_valid2 = flare_err[valid2]

    chi22 = np.sum(((flare_window_valid2 - model_rate_valid2) / flare_err_valid2) ** 2)
    dof2 = len(flare_window_valid2) - len(popt2)
    reduced_chi2_2 = chi22 / dof2

    

    print(f'---- FLARE {flare} ---- {i}')
    print(f"SINGLE PEAK Reduced Chi-squared: {reduced_chi2:.2f}")
    print(f"DOUBLE PEAK Reduced Chi-squared: {reduced_chi2_2:.2f}")
    if reduced_chi2 < reduced_chi2_2:
        print('Single Peak Optimal')
    elif reduced_chi2_2 < reduced_chi2:
        print('Double Peak Optimal')


    plt.plot(np.arange(-4200, 4201, 300), model_rate, label='Single Peak Model')
    plt.plot(np.arange(-4200, 4201, 300), model_rate2, label='Double Peak Model')
    plt.plot(np.arange(-4200, 4201, 300), flare_window, label='Data', c='black')
    plt.legend()
    plt.show()
