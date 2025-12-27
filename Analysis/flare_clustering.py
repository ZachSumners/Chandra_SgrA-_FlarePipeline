import os
import numpy as np
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit

df = pd.read_csv('flare_properties.csv')
flare_ids = df['obs_id']
max_rate = df['rate_max']
mjd_starts = df['start_mjd']
mjd_ends = df['end_mjd']
j = 0

quiescent_df = pd.read_csv('quiescent_rates.csv')
all_plots = []

for i, flare in enumerate(flare_ids):
    flare = str(flare)
    if flare == '1561':
        lc = f"./1561/1561_000/repro/01561_sgra_2-8keV_lc300_pileup.fits"
    elif len(flare) < 5:
        flare_5digit = str(flare).zfill(5)
        lc = f"./{flare}/repro/{flare_5digit}_sgra_2-8keV_lc300_pileup.fits"
    else:
        lc = f"./{flare}/repro/{flare}_sgra_2-8keV_lc300_pileup.fits"
    f = fits.open(lc)
    data = f[1].data 

    chandra_time = data['TIME']
    time = 50814.0 + (chandra_time / 86400.0)
    rate = data['RATE_PILEUP']
    mjd_start_flare = mjd_starts[i]
    mjd_end_flare = mjd_ends[i]
    max_rate_reported = max_rate[i]


    mjd_start = 51844.0
    sec_per_day = 86400

    # Identify points in the flare window
    flare_mask = (time >= mjd_start_flare) & (time <= mjd_end_flare)
    flare_times = time[flare_mask]
    flare_rates = rate[flare_mask]


    max_idx = np.where(np.around(flare_rates, decimals = 4) == max_rate_reported)
    
    if len(flare_times) == 0:
        print(f"No lightcurve data inside flare {i+1}")
        continue

    # Find brightest datapoint
    peak_idx = np.argmax(flare_rates)
    peak_time = flare_times[peak_idx]
    peak_rate = flare_rates[peak_idx]

    # Extract ±4200s around the peak
    window_mask = (time >= peak_time - 4200/86400) & (time <= peak_time + 4200/86400)
    
    flare_window = rate[window_mask]#/peak_rate

    peak_idx_flare_window = np.where(flare_window == peak_rate)

    if len(flare_window) != 29:
        if len(peak_idx_flare_window) == 1:
            flare_window = np.concatenate([np.zeros(14 - peak_idx_flare_window[0][0]), flare_window])
            flare_window = np.concatenate([flare_window, np.zeros(29 - len(flare_window))])
        else:
            peak_rate_single_flare_idx = np.min(np.abs(14 - peak_idx_flare_window))[0]
            flare_window = np.concatenate([np.zeros(14 - peak_rate_single_flare_idx), flare_window])
            flare_window = np.concatenate([flare_window, np.zeros(29 - len(flare_window))])

    quiescent_rate = quiescent_df.loc[quiescent_df['obs_id'] == int(flare), 'quiescent_rate'].values
    flare_window -= quiescent_rate[0]/1000

    flare_window = flare_window/(peak_rate-quiescent_rate[0]/1000)
    #flare_window = (flare_window - np.min(flare_window)) / (np.max(flare_window) - np.min(flare_window))
    
    all_plots.append(flare_window)

    flare_window = np.where(flare_window > 1, 0, flare_window)
    #plt.plot(np.arange(-4200, 4201, 300), flare_window, c='black', linewidth=1, alpha=0.3)


X = np.vstack(all_plots)
X_centered = X - np.mean(X, axis=0)

pca2 = PCA(2)
Xp2 = pca2.fit_transform(X_centered)

plt.scatter(Xp2[:, 0], Xp2[:, 1], s=3)
for i, obsid in enumerate(flare_ids):
    plt.text(Xp2[i, 0] + 0.02, Xp2[i, 1], str(obsid), fontsize=6, alpha=0.7)
plt.show()
