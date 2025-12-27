import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from lightcurve_extract import extract_lightcurve, extract_lightcurve_magnetar, extract_lightcurve_grating
import subprocess
import os

df = pd.read_csv('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/flare_properties_nov19.csv')
flare_ids = df['obs_id']
flare_numbers = df['flare_number']
start = df['start_mjd']
end = df['end_mjd']
mean_rate = df['rate_mean']
mean_rate_err = df['rate_mean_err']
max_rate = df['rate_max']
max_rate_err = df['rate_max_err']
duration = df['duration_s']
fluence = df['fluence_ct']

wd = os.getcwd()

medians = []
median_errs = []
snrs = []
for i in range(len(flare_ids)):
    observationID_5digit = str(flare_ids[i]).zfill(5)
    hdul = fits.open(f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/{observationID_5digit}_sgra_2-8keV_lc300_pileup.fits")
    data = hdul[1].data 
    
    mask = ((50814 + data['time']/86400) > float(start[i])) & ((50814 + data['time']/86400) < float(end[i]))
    flare = data[mask]['RATE_PILEUP']
    quiescent = data[~mask]['RATE_PILEUP']

    median_cr = np.median(flare)
    medians.append(median_cr)

    median_err = 1.253*np.std(flare)/np.sqrt(len(flare))
    median_errs.append(median_err)

    snrs.append(np.mean(flare)/(np.std(quiescent)))
    print(np.std(quiescent))

df['median_rate'] = medians
df['median_rate_err'] = median_errs
df['snr'] = snrs
df.to_csv("flare_properties_nov19.csv", index=False)
