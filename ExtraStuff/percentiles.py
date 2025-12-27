import pandas as pd
import numpy as np
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
fluence_err = mean_rate_err * duration

wd = os.getcwd()

hr = np.zeros(len(flare_ids))
hr_err = np.zeros(len(flare_ids))

compacts = []

for i in range(len(flare_ids)):

    #if flare_ids[i] == 1561:
    #    hr[i] = 0
    #    continue
    observationID_5digit = str(flare_ids[i]).zfill(5)

    if flare_ids[i] >= 13838 and flare_ids[i] <= 14468:
        subprocess.call(f'dmcopy infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_bary_evt2.fits[EVENTS][sky=region(/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/order1and0.reg)]" outfile="{observationID_5digit}_sgra_evt2.fits" clobber=yes option="all"', shell=True, cwd=wd)
        hdul = fits.open(f'{observationID_5digit}_sgra_evt2.fits')
        data = hdul[1].data

        flare_evts_mask = ((50814 + data['time']/86400) > float(start[i])) & ((50814 + data['time']/86400) < float(end[i]))
        mask = (data['energy'] >= 2000) & (data['energy'] <= 8000)

        evts = data[flare_evts_mask*mask]
        err_evts = np.sqrt(len(evts))

        t = evts['time']

        t5 = np.percentile(t, 5)
        t95 = np.percentile(t, 95)
        T90 = t95 - t5

        compact = (86400*(end[i] - start[i]))/T90
        compacts.append(compact)
        print(compact, 86400*(end[i] - start[i]), T90, 'here2')
    else:
        subprocess.call(f'dmcopy infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_bary_evt2.fits[EVENTS][sky=region(/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/sgra.reg)]" outfile="{observationID_5digit}_sgra_evt2.fits" clobber=yes option="all"', shell=True, cwd=wd)
        hdul = fits.open(f'{observationID_5digit}_sgra_evt2.fits')
        data = hdul[1].data

        flare_evts_mask = ((50814 + data['time']/86400) > float(start[i])) & ((50814 + data['time']/86400) < float(end[i]))
        mask = (data['energy'] >= 2000) & (data['energy'] <= 8000)

        evts = data[flare_evts_mask*mask]
        err_evts = np.sqrt(len(evts))

        t = evts['time']

        t5 = np.percentile(t, 5)
        t95 = np.percentile(t, 95)
        T90 = t95 - t5

        compact = (86400*(end[i] - start[i]))/T90
        compacts.append(compact)
        print(compact, 86400*(end[i] - start[i]), T90, 'here')


    

#df['hardness_ratio'] = hr
#df['hardness_ratio_err'] = hr_err
#df.to_csv("flare_properties_nov19.csv", index=False)

plt.figure(figsize=(8, 8))
plt.errorbar(df['obs_id'], compacts, ecolor='gray', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='steelblue', markersize=6,)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('ID')
plt.ylabel('Compacts')
plt.grid()
plt.show()
