from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/flare_properties.csv')
flare_ids = df['obs_id']
max_rate = df['rate_max']
mjd_starts = df['start_mjd']
mjd_ends = df['end_mjd']

for i, flare in enumerate(flare_ids):
    flare = str(flare)
    if len(flare) < 5:
        flare_5digit = str(flare).zfill(5)
        lc = f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare}/repro/{flare_5digit}_sgra_2-8keV_lc300_pileup.fits"
    else:
        lc = f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare}/repro/{flare}_sgra_2-8keV_lc300_pileup.fits"
    f = fits.open(lc)
    data = f[1].data 

    chandra_time = data['TIME']
    time = 50814.0 + (chandra_time / 86400.0)
    rate = data['RATE_PILEUP']
    err = data['PILEUP_ERR']
    mjd_start_flare = mjd_starts[i]
    mjd_end_flare = mjd_ends[i]
    max_rate_reported = max_rate[i]

    flare_mask = (time >= mjd_start_flare - (300*5)/86400) & (time <= mjd_end_flare + (300*5)/86400)
    flare_times = time[flare_mask]
    flare_rates = rate[flare_mask]
    flare_errors = err[flare_mask]

    plt.errorbar(flare_times, flare_rates, yerr=flare_errors, marker='o', color="red", mfc="black",mec="black", ecolor="grey")
    plt.title(f'Flare {i}, OBS ID {flare}')
    plt.gca().set_ylim(bottom=-0.01)
    plt.axvline(x=mjd_start_flare, color='black', linestyle='--') 
    plt.axvline(x=mjd_end_flare, color='black', linestyle='--') 
    plt.axhspan(0, 0.005*3, facecolor='grey', alpha=0.3)
    plt.show()