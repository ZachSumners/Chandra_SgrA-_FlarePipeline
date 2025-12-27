import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from lightcurve_extract import extract_lightcurve, extract_lightcurve_magnetar, extract_lightcurve_grating
import subprocess
import os

df = pd.read_csv('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/flare_properties.csv')
df = df[df['rate_max'] > 0.2].reset_index(drop=True)

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


binning = 300/86400

for i in range(len(flare_ids)):

    #if flare_ids[i] == 1561:
    #    hr[i] = 0
    #    continue
    hr_timeseries = []
    hr_timeseries_err = []
    lc_timeseries = []
    observationID_5digit = str(flare_ids[i]).zfill(5)
    observationID = flare_ids[i]
    fileName = 'bary'
    erange_low = [2, 4]
    erange_high = [4, 8]
    tbin = 300

    try:
        subprocess.call(f'dmextract infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_{fileName}_evt2.fits[energy={int(erange_low[0])*1000}:{int(erange_low[1])*1000},sky=region("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/sgra.reg")][bin time=::{tbin}]" outfile="{observationID}_sgra_{erange_low[0]}-{erange_low[1]}keV_lc{tbin}.fits" bkg="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_{fileName}_evt2.fits[sky=region("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/bkg.reg")]" opt="ltc1" clobber = yes', shell=True, cwd=wd)
        subprocess.call(f'dmextract infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_{fileName}_evt2.fits[energy={int(erange_high[0])*1000}:{int(erange_high[1])*1000},sky=region("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/sgra.reg")][bin time=::{tbin}]" outfile="{observationID}_sgra_{erange_high[0]}-{erange_high[1]}keV_lc{tbin}.fits" bkg="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_{fileName}_evt2.fits[sky=region("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/bkg.reg")]" opt="ltc1" clobber = yes', shell=True, cwd=wd)
        
        hdul_low = fits.open(f'{observationID}_sgra_{erange_low[0]}-{erange_low[1]}keV_lc{tbin}.fits')
        data_low = hdul_low[1].data['COUNT_RATE']
        time = hdul_low[1].data['TIME']
        hdul_high = fits.open(f'{observationID}_sgra_{erange_high[0]}-{erange_high[1]}keV_lc{tbin}.fits')
        data_high = hdul_high[1].data['COUNT_RATE']


        starttime = start[i]
        endtime = end[i]
        counter = 0

        flare_mask = ((50814 + time/86400) > float(starttime)) & ((50814 + time/86400) < float(endtime))

    except:
        subprocess.call(f'dmextract infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_{fileName}_evt2.fits[energy={int(erange_low[0])*1000}:{int(erange_low[1])*1000},sky=region("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/order1and0.reg")][bin time=::{tbin}]" outfile="{observationID}_sgra_{erange_low[0]}-{erange_low[1]}keV_lc{tbin}.fits" opt="ltc1" clobber = yes', shell=True, cwd=wd)
        subprocess.call(f'dmextract infile="/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/acisf{observationID_5digit}_{fileName}_evt2.fits[energy={int(erange_high[0])*1000}:{int(erange_high[1])*1000},sky=region("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/order1and0.reg")][bin time=::{tbin}]" outfile="{observationID}_sgra_{erange_high[0]}-{erange_high[1]}keV_lc{tbin}.fits" opt="ltc1" clobber = yes', shell=True, cwd=wd)
        
        hdul_low = fits.open(f'{observationID}_sgra_{erange_low[0]}-{erange_low[1]}keV_lc{tbin}.fits')
        data_low = hdul_low[1].data['COUNT_RATE']
        time = hdul_low[1].data['TIME']
        hdul_high = fits.open(f'{observationID}_sgra_{erange_high[0]}-{erange_high[1]}keV_lc{tbin}.fits')
        data_high = hdul_high[1].data['COUNT_RATE']


        starttime = start[i]
        endtime = end[i]
        counter = 0

        flare_mask = ((50814 + time/86400) > float(starttime)) & ((50814 + time/86400) < float(endtime))
        

    fig, axs = plt.subplots(2, 1, figsize=(8, 6))

    hr = data_high/(data_low + 1e-5)
    hr = np.where(hr > 5, np.nan, hr)

    lc = data_high + data_low

    axs[0].plot(time[flare_mask], hr[flare_mask], marker='.')
    axs[1].plot(time[flare_mask], lc[flare_mask], marker='.')
    axs[0].set_ylim(0, 5)
    axs[0].grid()
    axs[1].grid()
    axs[1].set_xlabel('Time (s)')
    axs[0].set_ylabel('Hardness Ratio')
    axs[1].set_ylabel('Count Rate (Not Pileup Corrected)')
    axs[0].set_title(f'OBS ID: {flare_ids[i]}')
    plt.show()

    

#df['hardness_ratio'] = hr
#df['hardness_ratio_err'] = hr_err
#df.to_csv("flare_properties_hr.csv", index=False)


