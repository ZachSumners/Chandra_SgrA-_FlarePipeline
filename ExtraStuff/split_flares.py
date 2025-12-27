import numpy as np
import os
import glob
from astropy.io import fits
from astropy.time import Time
from datetime import datetime, timedelta
import matplotlib.dates as mdates

from astropy.io import fits
from astropy.table import Table

def read_block_info(filename):
    """
    Read a Chandra BayesianBlocks info file and return a list of 
    rows, each containing six float values.
    """
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            # skip empty lines or comments
            if not line or line.startswith("#"):
                continue
            # split into fields and convert to float
            parts = line.split()
            if len(parts) < 6:
                continue
            vals = list(map(float, parts[:6]))
            rows.append(vals)
    return rows


def gather_all_obsid_results(base_dir, obsID):
    """
    Look for all *_sgrA_bayesianBlocks_info.txt files 
    inside base_dir/**/ directories.
    Returns dict: { obsid : [ [vals], [vals], ... ] }
    """
    results = {}

    # Pattern assumes each ObsID folder contains the file:
    pattern = os.path.join(base_dir, "Results", f"0{obsID}_sgra_bayesianBlocks_info_pileup.txt")
    files = glob.glob(pattern, recursive=True)

    for f in files:
        # extract ObsID from folder name
        obsid = os.path.basename(os.path.dirname(f))
        rows = read_block_info(f)
        results[obsid] = rows

    return results


def characterizeFlare(bp, evtfile, lcfile, gratingtype):
    f = fits.open(evtfile)

    primary = f[0].header
    data_info = f[1]
    evt = Table.read(f[1])
    gti = Table.read(f[2])
    
    g = fits.open(lcfile)
    dat = Table.read(g, hdu=1)
    count_rate_lc = dat['RATE_PILEUP']
    err = dat['PILEUP_ERR']

    mjdref_lc = primary['MJDREF']
    timezero_lc = primary['TIMEZERO']
    time_lc = mjdref_lc + (timezero_lc + dat["TIME"] ) / 86400

    start = bp[0] #start time of the flare
    end = bp[1]  #end time of the flare 

    '''Picking out the mximum count rate from lightcurve:  '''
    #within this time region, need to find the max count rate:
    upper = np.where(time_lc <= bp[1])
    lower = np.where(time_lc >= bp[0]) #indices that say where 
    #within the light curve the lightcurves start/finish 


    #the maximum count rate can be found for times between indices min(lower) and max(upper)
    #then the maximum count rate is:

    if (np.min(lower) == np.max(upper)+1):
        ct_max = np.max(count_rate_lc[np.min(lower)])
        ct_max_err = err[np.where(count_rate_lc == ct_max)][0]
    else:
        ct_max = np.max(count_rate_lc[np.min(lower):np.max(upper)+1])
        ct_max_err = err[np.where(count_rate_lc == ct_max)][0]

    '''Picking out Luminosity & Energy emitted by each flare: '''
    #According to Eli's paper there is a proportionality between 
    #counts in a block and the energy 
    #for non gratings: 0.013 ct / (10^34 erg) 
    #then the energy is given by: 
    ct_mean = bp[4]/86400.0 #count rate in each block

    ct_mean_err = bp[5]/86400.0 #err in mean ct rate

    l = 0.005275
    ct_quies_sub = ct_mean - 0.005275 #quiescence subtracted countrates
    dur = bp[3]*86400.0 #duration of flare in s 

    if gratingtype == 'False':
        #Determining values
        luminosity = (ct_quies_sub*(10**(34)/0.0118766)) # erg/s  
        energy = luminosity * dur  #erg
        #Error propagation
        lum_err = (np.sqrt((ct_mean_err)**2 + (l/86400.0)**2))*(10**(34)/0.0008598192516)
        energy_err = energy * np.sqrt((lum_err/luminosity)**2)
    elif gratingtype == 'True':
        #Determining values
        luminosity = (ct_quies_sub*(10**(34)/0.008485)) # erg/s  
        energy = luminosity * dur  #erg
        #Error propagation
        lum_err = (np.sqrt((ct_mean_err)**2 + (l/86400.0)**2))*(10**(34)/0.0006965207202)
        energy_err = energy * np.sqrt((lum_err/luminosity)**2)


    #Finding the Flux of each flare: 
    #flux is related to the distance between sgrA and the detector 
    #which is about radius = 2.45*10^22 cm 
    #the surface area over which we are dividing for the 
    #flux is 4*pi*radius^2 
    radius = 2.523467125*10**(22)
    err_rad = 4 *10**(19)
    flux = ((luminosity)/(4*np.pi*radius**2))# 10^-12 erg/s/cm2

    #Error propagation: 
    flux_err = flux*np.sqrt((lum_err/luminosity)**2 + (err_rad/radius)**2)

    #Finding the fluence of each flare: 
    #fluence is given in ct so I am assuming it is the total 
    #count number incident upon the detector 
    fluence = ct_mean*bp[3]*86400.0

    print("-------------------------------------------------------------------------------------", "\n")
    print("Start Time: ", start, "(MJD) \n")
    print("End Time: ", end, "(MJD) \n")
    print("Duration: ", dur, "(s) \n")
    print("Count Rate (mean): ",np.around(ct_mean,4) , "+/-", np.around(ct_mean_err, decimals=4), "(ct/s) \n")
    print("Count Rate (max): ",np.around(ct_max,4) ,"+/-", np.around(ct_max_err,4), "(ct/s) \n")
    print("Energy: ", np.around(energy/10**(37), 4), "+/-", np.around(energy_err/10**(37),4), "10^37 ergs \n")
    print("Luminosity: ", np.around(luminosity/10**(34),4), "+/-", np.around(lum_err/10**(34),4), "10^34 erg/s \n")
    print("Flux: ", np.around(flux/(10**(-12)),4) ,"+/-" ,np.around(flux_err/(10**(-12)),4) , "10^-12 erg/s/cm2 \n")
    print("Fluence: ", np.around(fluence,4), "ct \n")


# ---- RUN ----

obsID = 1561
base = f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{obsID}/repro"   # or the path to your Chandra_SgrA_-_FlarePipelineNoGithub directory
obs_data = gather_all_obsid_results(base, obsID)

observationID_5digit = str(obsID).zfill(5)

gratingtype = 'False'

# Example: print everything
for obsid, rows in obs_data.items():
    print(f"ObsID {obsid}:")
    for r in rows:
        print("    ", r)

newflare = [51844.2409347652, 51844.2795123523, 232.0, 0.038578, 6013.856388, 864.468099]
characterizeFlare(newflare, f'{base}/acisf{observationID_5digit}_bary_evt2.fits', f'{base}/{observationID_5digit}_sgra_2-8keV_lc300_pileup.fits', gratingtype)