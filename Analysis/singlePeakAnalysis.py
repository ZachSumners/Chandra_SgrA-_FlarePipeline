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

all_plots = []

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

for i, flare in enumerate(flare_ids):
    flaretype = classify_flare(flare, flare_numbers[i])

    if flaretype=='single_peak':
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

        #Find the peak idx in the flare time (not of the total observation)
        peak_idx_flare_window = np.where(flare_window == peak_rate)

        #Need to grab 29 total datapoints - 14 on each side. Problem if flare is near the edge of the observation.
        if len(flare_window) != 29:
            if len(peak_idx_flare_window) == 1:
                flare_window = np.concatenate([np.zeros(14 - peak_idx_flare_window[0][0]), flare_window])
                flare_window = np.concatenate([flare_window, np.zeros(29 - len(flare_window))])
            else:
                peak_rate_single_flare_idx = np.min(np.abs(14 - peak_idx_flare_window))[0]
                flare_window = np.concatenate([np.zeros(14 - peak_rate_single_flare_idx), flare_window])
                flare_window = np.concatenate([flare_window, np.zeros(29 - len(flare_window))])

        
        #Normalize by peak rate
        flare_window = flare_window/(peak_rate)
        
        #Crop the flare for a few observations that have flares show up in the 70 minutes on each side
        if flare == str(13854) or (flare == str(15042) and flare_numbers[i] == 1):
            flare_window[:9] = 0
            flare_window[19:] = 0

        all_plots.append(flare_window)

        #flare_window = np.where(flare_window > 1, 0, flare_window)
        plt.plot(np.arange(-4200, 4201, 300), flare_window, c='black', alpha = 0.3, linewidth=1)#, label=f'{flare} - {flare_numbers[i]}')

all_flares_stacked = np.vstack(all_plots)

def exp_decay(t, tau, C):
    return np.exp(-t / tau) + C
def exp_rise(t, tau, C):
    return np.exp(t / tau) + C

def exp_fit(average_flare):
    # Fit the model
    popt, _ = curve_fit(exp_decay, np.arange(0, 4201, 300), average_flare[14:], p0=[1000, 0])
    tau_fit, C = popt

    time_rise = np.arange(-4200, 1, 300)
    popt2, _ = curve_fit(exp_rise, time_rise, average_flare[:15], p0=[500, 0])
    tau_fit2, C2 = popt2

    return tau_fit, C, tau_fit2, C2


def uncert_bootstrap(flare_array):
    sample_indices = np.random.choice(flare_array.shape[0], size=21, replace=True)
    flare_array_sampled = flare_array[sample_indices]   

    average_flare = np.median(flare_array_sampled, axis=0)

    tau_fit, C, tau_fit2, C2 = exp_fit(average_flare)
    return tau_fit, C, tau_fit2, C2

tau_fits = np.zeros(100)
tau_fits2 = np.zeros(100)
for i in range(100):
    tau_fit, C, tau_fit2, C2 = uncert_bootstrap(all_flares_stacked)
    tau_fits[i] = tau_fit
    tau_fits2[i] = tau_fit2

average_flare = np.median(all_flares_stacked, axis=0)
tau_fit, C, tau_fit2, C2 = exp_fit(average_flare)

time_rise = np.arange(-4200, 1, 300)
plt.plot(np.arange(-4200, 4201, 300), average_flare, c='red', linewidth=3, label='Median')
#plt.plot(np.arange(-4200, 4201, 300), np.median(X, axis=0), c='blue', linewidth=3, label='median')
plt.plot(np.arange(0, 4201, 300), exp_decay(np.arange(0, 4201, 300), tau_fit, C), c='orange', linewidth=3, label=f'Fit to median, decay factor = {tau_fit/60:.2f} ± {np.std(tau_fits)/60:.2f} min')
plt.plot(time_rise, exp_rise(time_rise, tau_fit2, C2), c='orange', linewidth=3, label=f'Fit to median, rise factor = {tau_fit2/60:.2f} ± {np.std(tau_fits2)/60:.2f} min')#plt.ylim(-0.05, 1.05)
plt.legend()
plt.xlabel('Time (seconds)')
plt.ylabel('Normalized Flux')
plt.title('Characteristic Flare Shape')
plt.show()