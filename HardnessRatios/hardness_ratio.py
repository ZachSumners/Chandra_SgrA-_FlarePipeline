import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from lightcurve_extract import extract_lightcurve, extract_lightcurve_magnetar, extract_lightcurve_grating
import subprocess
import os

df = pd.read_csv('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/flare_properties_dec2.csv')
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
fluence_err = mean_rate_err * duration

wd = os.getcwd()

hr = np.zeros(len(flare_ids))
hr_err = np.zeros(len(flare_ids))

for i in range(len(flare_ids)):

    #if flare_ids[i] != 14463:
    #    hr[i] = 0
    #    continue
    observationID_5digit = str(flare_ids[i]).zfill(5)

    try:
        subprocess.call(f'dmcopy infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_bary_evt2.fits[EVENTS][sky=region(/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/sgra.reg)]" outfile="{observationID_5digit}_sgra_evt2.fits" clobber=yes option="all"', shell=True, cwd=wd)
        hdul = fits.open(f'{observationID_5digit}_sgra_evt2.fits')
        data = hdul[1].data

        flare_evts_mask = ((50814 + data['time']/86400) > float(start[i])) & ((50814 + data['time']/86400) < float(end[i]))
        soft_mask = (data['energy'] > 2000) & (data['energy'] < 4000)
        hard_mask = (data['energy'] >= 4000) & (data['energy'] <= 8000)


        soft_evts = data[flare_evts_mask*soft_mask]
        err_soft = np.sqrt(len(soft_evts))
        hard_evts = data[flare_evts_mask*hard_mask]
        err_hard = np.sqrt(len(hard_evts))
        
        print(flare_ids[i], len(hard_evts), len(soft_evts))
        hr[i] = len(hard_evts)/len(soft_evts)
        hr_err[i] = hr[i]*np.sqrt((err_hard/len(hard_evts))**2 + (err_soft/len(soft_evts))**2)
        
        print(flare_ids[i], err_hard, err_soft, len(hard_evts), len(soft_evts), hr[i], hr_err[i])
    except:
        subprocess.call(f'dmcopy infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_bary_evt2.fits[EVENTS][sky=region(/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/order1and0.reg)]" outfile="{observationID_5digit}_sgra_evt2.fits" clobber=yes option="all"', shell=True, cwd=wd)
        hdul = fits.open(f'{observationID_5digit}_sgra_evt2.fits')
        data = hdul[1].data

        flare_evts_mask = ((50814 + data['time']/86400) > float(start[i])) & ((50814 + data['time']/86400) < float(end[i]))
        soft_mask = (data['energy'] > 2000) & (data['energy'] < 4000)
        hard_mask = (data['energy'] >= 4000) & (data['energy'] <= 8000)


        soft_evts = data[flare_evts_mask*soft_mask]
        err_soft = np.sqrt(len(soft_evts))
        hard_evts = data[flare_evts_mask*hard_mask]
        err_hard = np.sqrt(len(hard_evts))
        
        print(flare_ids[i], len(hard_evts), len(soft_evts))
        hr[i] = len(hard_evts)/len(soft_evts)
        hr_err[i] = hr[i]*np.sqrt((err_hard/len(hard_evts))**2 + (err_soft/len(soft_evts))**2)
        
        print(flare_ids[i], err_hard, err_soft, len(hard_evts), len(soft_evts), hr[i], hr_err[i])

print(hr, hr_err)

df['hardness_ratio'] = hr
df['hardness_ratio_err'] = hr_err
df.to_csv("flare_properties_dec2.csv", index=False)

plt.figure(figsize=(8, 8))
plt.errorbar(hr, fluence, xerr=hr_err, yerr=fluence_err , ecolor='gray', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='steelblue', markersize=6,)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Hardness Ratio')
plt.ylabel('Fluence (ct)')
plt.grid()
plt.show()
