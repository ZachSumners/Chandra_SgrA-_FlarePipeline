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
j = 0

quiescent_df = pd.read_csv('quiescent_rates.csv')
all_plots = []

singlepeak_ids = [11843, 10556, 13838, 13839, 13842, 13854, 13852, 13857, 14463, 14466, 15041, 15042, 16963, 16508, 18731, 20446, 21454, 22230, 23666, 23665, 23739, 28231, 3393, 4684, 6363]

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
    #flare_window -= quiescent_rate[0]/1000

    #flare_window = flare_window/(peak_rate-quiescent_rate[0]/1000)
    flare_window = flare_window/(peak_rate)
    #flare_window = (flare_window - np.min(flare_window)) / (np.max(flare_window) - np.min(flare_window))
    
    if (flare == str(13854)) or (flare == str(15042) and flare_numbers[i] == 1):
        flare_window[:9] = 0
        flare_window[19:] = 0

    all_plots.append(flare_window)

    flare_window = np.where(flare_window > 1, 0, flare_window)
    #plt.plot(np.arange(-4200, 4201, 300), flare_window, linewidth=1, c='black', alpha=0.3)#, label=f'{flare} - {flare_number}')


X = np.vstack(all_plots)
X_centered = X - np.mean(X, axis=0)
pca = PCA(n_components=2)
X_PCA = pca.fit_transform(X_centered)






def classify_flare(obs_id, flare_number):
    key = (obs_id, flare_number)

    single_peak = {
        (10556, 1), (10556, 3), (10556, 4), (13838, 1), (13839, 1),
        (13842, 2), (13854, 1), (13854, 2), (13854, 3), (13854, 4), (13852, 1), (13852, 2),
        (13857, 1), (14463, 1), (14466, 1), (15041, 2), (15042, 1), (15045, 2),
        (16963, 1), (16508, 1), (18731, 1), (20446, 1), (21454, 1), (22230, 1), (23666, 1),
        (23665, 1), (23739, 1), (28231, 1), (3393, 3), (4684, 1), (6363, 1), (1561, 1)
    }

    double_peak = {
        (10556, 2), (13843, 1), (13847, 1), (13851, 1), (13849, 1), (14462, 2), (14468, 1),
        (14944, 1), (15043, 1), (15042, 2), (15570, 1), (16966, 1), (18732, 1), (20041, 1),
        (20346, 1), (21454, 2), (21456, 1), (25297, 1), (26760, 1), (3393, 1), (3393, 2),
        (5953, 1), (22595, 2), (22595, 1), (14392, 1), (13839, 2)
    }

    higher_substructure = {
        (13842, 1), (14432, 2), (14392, 2), (16218, 1), (20751, 1), (22594, 1), (22595, 3), (22592, 1), (1561, 2), (11843, 1)
    }

    ambiguous = {
        (13016, 1), (13840, 1), (13842, 3), (13845, 1), (14427, 1), (14427, 2),
        (14462, 1), (14439, 1), (14465, 1), (14465, 2), (14468, 2), (15041, 1), (15045, 1),
        (16217, 1), (18055, 1), (19703, 1), (19703, 2),
        (22937, 1), (26759, 1), (3392, 1), (3392, 2), (3392, 3),
        (3663, 1), (5952, 1), (9173, 1), (18055, 2), (14432, 1)
    }
    
    if key in single_peak:
        return "single_peak"
    elif key in double_peak:
        return "double_peak"
    elif key in higher_substructure:
        return "higher_substructure"
    elif key in ambiguous:
        return "ambiguous"
    else:
        return "unclassified"

plotted = {
    'single_peak': False,
    'double_peak': False,
    'higher_substructure': False,
    'ambiguous': False,
    'unclassified': False
}

for i, flare in enumerate(flare_ids):
    flaretype = classify_flare(flare, flare_numbers[i])

    label = None
    if not plotted[flaretype]:
        label = flaretype.replace('_', ' ').title()
        plotted[flaretype] = True
    

    if flaretype=='single_peak':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='r', ecolor="black", alpha=0.8, label=label)
        #single.append(flare)
        plt.scatter(X_PCA[i, 0], X_PCA[i, 1], s=12, c='r', alpha=0.8, label=label)
    elif flaretype=='double_peak':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='b', ecolor="black", alpha=0.8, label=label)
        #double.append(flare)
        plt.scatter(X_PCA[i, 0], X_PCA[i, 1], s=12, c='b', alpha=0.8, label=label)
    elif flaretype=='higher_substructure':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='g', ecolor="black", alpha=0.8, label=label)
        #higher.append(flare)
        plt.scatter(X_PCA[i, 0], X_PCA[i, 1], s=12, c='g', alpha=0.8, label=label)
    elif flaretype=='ambiguous':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='orange', ecolor="black", alpha=0.8, label=label)
        #ambiguous.append(flare)
        plt.scatter(X_PCA[i, 0], X_PCA[i, 1], s=12, c='orange', alpha=0.8, label=label)
    else:
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='gray', ecolor="black", alpha=0.8, label=str(flare))
        #unclassified.append(flare)
        plt.scatter(X_PCA[i, 0], X_PCA[i, 1], s=12, c='gray', alpha=0.8, label=label)
   

plt.show()

