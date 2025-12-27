#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 21:22:24 2020

@author: Elisa

This is a modification of Hope Boyce's and Eli Bouffard's algorithm 
In turn, each of these was based on Peter William's work. 

In this version of the file, I contain all the necessary functions to 
get all the information we want from the bayseian block analysis 
for the x-ray flare database. 

Functions in this code are for plotting & creating/writing files.

Needs to be used in conjunction with xbblocks.py file. 

FUNCTIONS: 
    
1. Process ->> creates an 'output file' with bb blocks info. 
2. plot_bb ->> plots the  bayesian blocks on axes
3. plot_lc ->> plots the lightcurve on axes 
3. getInfo ->> creates an output file that contains all the
                information needed to populate the flare database.
 
"""

import matplotlib.pyplot as plt
import sys
import numpy as np
import numpy.ma as ma
import math

from astropy.time import Time
from datetime import datetime, timedelta
import matplotlib.dates as mdates

from astropy.io import fits
from astropy.table import Table

import xbblocks #this is a set of functions written by Eli Bouffard
from marx_pileup_simulation import marx_pileup_interpolation_block

def nongrating_equation(x, alpha=0.5, dividezero_epsilon=1e-8):
	fd = 1-(((np.exp(alpha*x)-1)*np.exp(-alpha*x))/(alpha*x + dividezero_epsilon))
	observed_count_rate = x/(1-(fd + dividezero_epsilon))

	return observed_count_rate

def nongrating_intersection(x):
	return nongrating_equation(x) - count_rate_equation(x)

def grating_equation(x):
	K = 0.94
	alpha = 1
	beta = 0.52

	y = x * (1 + K/(alpha * x)*(np.exp(alpha*x) - 1)*np.exp(-x))*(beta)

	return y

def intersection(x):
	return grating_equation(x) - count_rate_equation(x)

def block_pileup_correction(rate, exptime):
    rate = rate/86400 * exptime

    #Additional factor so no division by 0.
    dividezero_epsilon = 1e-8

    pileup_rate = np.zeros(len(count_rate))
    for i, count_rate_point in enumerate(count_rate):
        if count_rate_point == 0.0:
            pileup_rate[i] = count_rate_point
        else:
            intersection = lambda x: nongrating_equation(x) - count_rate_point

            result = root_scalar(intersection, bracket=[0.00001, 2], method='brentq')
            pileup_rate[i] = result.root

    #Need to derive pileup error rate
    pileup_rate = pileup_rate*1/exptime
    return pileup_rate


def block_pileup_correction_grating(rate, exptime, repro_wd):
    #Additional factor so no division by 0.
    dividezero_epsilon = 1e-8

    pileup_rate = np.zeros(len(count_rate))
    for i, count_rate_point in enumerate(count_rate):
        if count_rate_point < 0.06:
            pileup_rate[i] = count_rate_point
        else:
            intersection = lambda x: grating_equation(x) - count_rate_point

            result = root_scalar(intersection, bracket=[0.00001, 2], method='brentq')
            pileup_rate[i] = result.root

    #Need to derive pileup error rate
    pileup_rate = pileup_rate*1/exptime
    return pileup_rate


def block_pileup_correction_marx(rate, exptime, repro_wd):
    rate = rate/86400

    marx_conversions = np.loadtxt(f"{repro_wd}/marx_pileup_conversion.txt")
    marx_observed_flux = marx_conversions[:, 0]
    marx_true_flux = marx_conversions[:, 1]
    pileup_rate = marx_pileup_interpolation_block(marx_observed_flux, marx_true_flux, rate)

    return pileup_rate * 86400

def process(infile, outfile, pileup_outfile, pileup_correction, repro_wd, p0=0.05):
    
    ''' 
    Generates the bayesian blocks and creates file with results 
    
    Parameters
    ----------
    
        infile   : name of input fits file (event file)
        
        p0       : default is 0.05
        
        outfile  : name of output text file into which results are 
                   written. 
                   
    Returns
    ----------
    Nothing, but creates a file with bayesian block results. 
    These results contain: 
        p0 value 
        timesys 
        start time of observation (MJD)
        end time of observation (MJD)
        time array length 
        
        Block info in order: 
            left edge of block mjd
            right edge of block  mjd
            counts in the block
            width of the block mjd
            count rate (avg) in block ct/day 
            count rate error ct/day 
    
    '''
    
    
    f = fits.open(infile)
    o = open(outfile, mode='w')
    op = open(pileup_outfile, mode='w')

    timesys = f[0].header['timesys']
    mjdref = f[0].header['mjdref']
    
    timeunit = f[0].header['timeunit']
    timezero = f[0].header['timezero']
    tstart = f[0].header['tstart']
    tstop = f[0].header['tstop']

    if timeunit == 's':
        tscale = 1. / 86400 #to convert s to days 
        #24*60*60 = 86400s/day so that timestamps can be added to mjdref
    else:
        die ('can\'t handle time unit "%s" in input "%s"', timeunit)

    eventhdu = None

    for hdu in f[1:]:
        if hdu.name == 'EVENTS':
            if eventhdu is None:
                eventhdu = hdu
            else:
                die('input "%s" has multiple EVENTS sections; don\'t know '
                     'which to use',)

    if eventhdu is None:
        die('input "%s" has no EVENTS sections')

    values, counts = np.unique(eventhdu.data.ccd_id, return_counts=True)
    most_common_ccd = values[np.argmax(counts)]
    
    ccdid = most_common_ccd#eventhdu.data.ccd_id.min ()
    #if eventhdu.data.ccd_id.max () != ccdid:
    #    die ('can\'t handle data from multiple CCDs in input "%s"')

    gtihdu = None

    for hdu in f[1:]:
        if hdu.name == 'GTI' and hdu.header['ccd_id'] == ccdid:
            if gtihdu is None:
                gtihdu = hdu
            else:
                die ('input "%s" has multiple matching GTI sections; don\'t know '
                     'which to use')

    if gtihdu is None:
        print (sys.stderr, 'warning: no GTI info for active CCD %d; trusting ' \
            'TSTART and TSTOP instead' % ccdid)
        tstarts = np.asarray ([tstart])
        tstops = np.asarray ([tstop])
    else:
        tstarts = gtihdu.data.START #defined in terms of good time interval
        tstops = gtihdu.data.STOP

    times = (eventhdu.data.time + timezero) * tscale + mjdref
    tstarts = (tstarts + timezero) * tscale + mjdref
    tstops = (tstops + timezero) * tscale + mjdref


    exptime = float(f[1].header['exptime'])
    info = xbblocks.bsttbblock (times, tstarts, tstops, exptime, pileup_correction, repro_wd, p0=p0, nbootstrap=256)


    ######### Drops blocks with 0 count and if they are too thin ##########
    #MIN_DURATION = 3/86400  # seconds
    keep = [True] * info.nblocks  # a boolean flag per block
    merged_counts = info.counts.copy()
    merged_exposure = [info.widths[i] for i in range(info.nblocks)]
    # note: if info.widths is already exposure in seconds, use that. 
    # Otherwise, convert lengths to actual exposure sums if needed.

    # First round: mark any block that is too thin
    #for i in range(info.nblocks):
    #    if info.widths[i] < MIN_DURATION:
    #        keep[i] = False
    #    if info.counts[i] == 0:
    #        keep[i] = False

    # Second round: for each block i that is “too thin,” merge its counts/widths into the closest neighbor
    for i in range(info.nblocks):
        if not keep[i]:
            # Decide whether to merge with block to the left (i-1) or right (i+1).
            left_exists  = (i-1 >= 0) and keep[i-1]
            right_exists = (i+1 < info.nblocks) and keep[i+1]

            if left_exists and right_exists:
                # Compare rates to decide which neighbor is “closer”
                left_rate  = info.counts[i-1] / info.widths[i-1]
                right_rate = info.counts[i+1] / info.widths[i+1]
                this_rate  = info.counts[i] / info.widths[i]

                # Merge into the neighbor whose rate is closer to this block’s rate
                if abs(left_rate - this_rate) <= abs(right_rate - this_rate):
                    merged_counts[i-1] += info.counts[i]
                    merged_exposure[i-1] += info.widths[i]
                else:
                    merged_counts[i+1] += info.counts[i]
                    merged_exposure[i+1] += info.widths[i]

            elif left_exists:
                merged_counts[i-1] += info.counts[i]
                merged_exposure[i-1] += info.widths[i]
            elif right_exists:
                merged_counts[i+1] += info.counts[i]
                merged_exposure[i+1] += info.widths[i]
            # If neither neighbor “exists” (edge case), just leave it—dropping will happen later.

    print('# p0 = %g' % p0, file=o)
    print('# timesys =', timesys, file=o)
    print('# tstarts =', ' '.join ('%.5f' % t for t in tstarts), file=o)
    print('#tstops  =', ' '.join ('%.5f' % t for t in tstops), file=o)
    print('# n = %d' % times.size, file=o)

    if pileup_outfile != None:
        print('# p0 = %g' % p0, file=op)
        print('# timesys =', timesys, file=op)
        print('# tstarts =', ' '.join ('%.5f' % t for t in tstarts), file=op)
        print('#tstops  =', ' '.join ('%.5f' % t for t in tstops), file=op)
        print('# n = %d' % times.size, file=op)

    # Finally, print only those blocks where keep[i] == True, using the merged stats:
    for i in range(info.nblocks):
        if not keep[i]:
            continue

        # The “merged” block has:
        #    counts = merged_counts[i]
        #    exposure (seconds) = merged_exposure[i]
        # so the rate becomes merged_counts[i] / merged_exposure[i].
        merged_rate = merged_counts[i] / merged_exposure[i]

        s = '%.10f %.10f %10d %10f %10f %10f' % (
            info.ledges[i],
            info.redges[i],
            merged_counts[i],        # total counts (including any merged‐in bits)
            merged_exposure[i],      # total duration in seconds
            merged_rate,             # updated rate
            info.bsrstds[i]          # you can either leave the old block STD or recompute
        )
        print(s, file=o)

        if pileup_outfile != None:
            if pileup_correction == 'analytic':
                merged_rate = block_pileup_correction(merged_rate, exptime)
                merged_counts[i] = merged_rate*merged_exposure[i]
            elif pileup_correction == 'marx':
                merged_rate = block_pileup_correction_marx(np.array([merged_rate]), exptime, repro_wd)
                merged_counts[i] = merged_rate*merged_exposure[i]
                
            sp = '%.10f %.10f %10d %10f %10f %10f' % (
                info.ledges[i],
                info.redges[i],
                merged_counts[i],        # total counts (including any merged‐in bits)
                merged_exposure[i],      # total duration in seconds
                merged_rate,             # updated rate
                info.bsrstds[i]          # you can either leave the old block STD or recompute
            )
            print(sp, file=op)

    o.close()
    if pileup_outfile != None:
        op.close()

    #print('# p0 = %g' % p0, file=o)
    #print('# timesys =', timesys, file=o)
    #print('# tstarts =', ' '.join ('%.5f' % t for t in tstarts), file=o)
    #print('#tstops  =', ' '.join ('%.5f' % t for t in tstops), file=o)
    #print('# n = %d' % times.size, file=o)
    #for i in range (info.nblocks):
    #    #s = '%.5f %.5f %4d %g %g %g' % (info.ledges[i], info.redges[i],
    #                                    #info.counts[i], info.widths[i],
    #                                    #info.rates[i], info.bsrstds[i])
    #    s = '%.10f %.10f %10f %10f %10f %10f' % (info.ledges[i], info.redges[i],
    #                                             info.counts[i], info.widths[i],
    #                                             info.rates[i], info.bsrstds[i])
        
    #    print(s, file=o)
    #o.close()



'''___________________________________________________________________
'''



def plot_bb(file, ax):
    '''
    Plot Bayesian Blocks onto a plot

    Parameters
    ----------   
    file     : string
        location and name of Bayesian Blocks results

    time_start : float
        MJD start of observation in days

    time_end : float
        MJD end of observation in days

    Returns
    -------
        nothing, plots the bayesian blocks. 

    '''
    ### read Bayesian Blocks data in 
    ledges, redges, counts, widths, rates, bsrstds = np.transpose(np.loadtxt(file))

    bstart = ledges ### rename
    bend = redges ### rename
    rates = rates/86400.0 #convert from counts/day to counts/s
    #####
    ### Recasting Bblocks results into arrays that are easier for plotting
    ### x = time array
    ### rates = the levels of blocks
    #####
    if hasattr(rates, "__len__"):
        ### If there is more than one Bblock
        x = [bstart[0]]
        for j in range(0,len(rates)):
            x= np.concatenate([x,[bstart[j]]])
        x= np.concatenate([x,[max(bend)]])
        rates = np.concatenate([[0],rates,[0]])
    else:
        ### If there's only one Bblock
        x = np.array([bstart,bstart,(bend)])

        rates = np.array([rates])
        rates = np.concatenate([[0],rates,[0]])
        

    dt = Time(x, format='mjd')
    datetime_block = dt.datetime

    ##x_utc = Time(x, format='mjd', scale='utc') ## convert time array into UTC format

    ###plotting


    ax.plot(datetime_block,rates,drawstyle='steps-post',color='#ff7b0f', lw=1.5,zorder=5) ### plot with MJD axis - timestart
    # plt.plot(x_utc.datetime,rates,drawstyle='steps-post',color='#ff7b0f', lw=1.5,zorder=5) ### plot with UTC axis
    return ax


'''___________________________________________________________________
'''



def plot_lc(file, rate_header, rate_err_header, ax):
    ''' Plot the lightcurve file 
    
    Parameters
    ----------
    
    file: fits file 
        File that contains the binned light curve data. 
        
    Since the bayesian blocks algorithm bins data by 300s, it is most convinent to use 300s 
    bins when plotting both on overlap. 
    '''
    
    f = fits.open(file)
    
    #obtain the mjdref from the primary header file : 
    #obtain the timezero correction from primary header: 
    primary = f[0].header

    try:
        mjdref = primary['MJDREF']
        timezero = primary['TIMEZERO']
    except KeyError:
        primary = f[1].header
        mjdref = primary['MJDREF']
        timezero = primary['TIMEZERO']
    
    #obtain necessary data for plotting: 
    data = Table.read(f, hdu=1)
    time = (mjdref + (timezero + data["TIME"] )/86400) #convert time data to mjd 
    counts = data[rate_header]
    err = data[rate_err_header]
    
    
    
    #plot : 
    dt = Time(time, format='mjd')
    datetime_ct = dt.datetime
    #ax.spines['right'].set_position(('outward', 80))

    plt.errorbar(datetime_ct, counts, yerr=err, marker=".", color="blue", mfc="black",mec="black", ecolor="navy")
    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)



'''___________________________________________________________________
'''



    
def getInfo(evtfile, lcfile, bbfile, bbfile_p, outfile, rate_header, rate_err_header, gratingtype, countmin=8, amp_crit=3):
    ''' 
    Generate a document that contains all the data required for the flare
    table database. 
    
    Parameters
    ----------
    evtfile: fits file
        File on which the bayesian blocks was run. 
        
    lcfile: fits file 
        File with the lightcurve of the observation. Ideally, bin 300s
        
    bbfile: txt file 
        File that was output of the process function. Contains block info
    
    outfile: str 
        Name of the output file to generate. 
    
    countmin: int 
        Minimum number of counts in a block. If lower, it will be combined with another block
        Default value is 8 *(arbitrary)
    
    amp_crit: int 
        Sigma range above quiesence for a block to be considered a flare
        Default value is 3
        
    
    Returns
    -------
    Nothing, but generates output file with information needed for the 
    database. 
    
    '''
    
    #Open the event file to get the basic information about observation
    
    f = fits.open(evtfile)
    o = open(outfile, mode='w')
    primary = f[0].header
    data_info = f[1]
    evt = Table.read(f[1])
    gti = Table.read(f[2])
    
    g = fits.open(lcfile)
    dat = Table.read(g, hdu=1)
    count_rate_lc = dat[rate_header]
    err = dat[rate_err_header]
    #convert times to mjd:
    
    mjdref_lc = primary['MJDREF']
    timezero_lc = primary['TIMEZERO']
    time_lc = mjdref_lc + (timezero_lc + dat["TIME"] ) / 86400
    
    mjdref_evt = primary['MJDREF']
    timezero_evt = primary['TIMEZERO']
    time_evt = mjdref_evt + (timezero_evt + evt["time"] ) / 86400
    
    ''' Get data from the primary header '''
    
    obsid = primary["OBS_ID"]
    date = primary["DATE-OBS"]
    instrument = primary["INSTRUME"]
    telescope = primary["TELESCOP"]
    exp = (gti["STOP"][0]-gti["START"][0])/1000 #the good time interval 
    
    #Add this data to the file: 
    print("-------------------------------------------------------------------------------------", "\n", file=o)
    print("BASIC INFORMATION", "\n", file=o)
    print("-------------------------------------------------------------------------------------", "\n", file=o)
    print("Obs_ID:", obsid, "\n", file=o)
    print("Obs Date:", date, "\n", file=o)
    print("Telescope:", telescope, "\n", file=o)
    print("Instrument:", instrument, "\n", file=o)
    print("Exposure (ks):", exp, "\n", file=o)
    
    '''Get data about the flares and Quiescent count rate'''
    
    #load bayesian blocks data 
    ledges, redges, counts, widths, rates, bsrstds = np.transpose(np.loadtxt(bbfile))
    #NOTE: the data that is given by 'process' function is in mjd 
    
    #Get the flare and block information:
    splits = None#[56516.3577665360 + 2360/86400, 56516.3577665360 + 4119/86400, 56516.3577665360 + 9109/86400]#
    manual_groups = None#[[1], [2], [4]]
    d,b,l = get_flare_bb_nobsnopcr(ledges, redges, counts, widths, rates, amplitude_criteria=3, minflu=5, manual_groups=manual_groups, split_times=splits)#get_flare_bb_nobsnopcr(ledges, redges, counts, widths, rates, amp_crit, countmin)
    
    if bbfile_p != None:
        ledges_p, redges_p, counts_p, widths_p, rates_p, bsrstds_p = np.transpose(np.loadtxt(bbfile_p))
        dp,bp,lp = get_flare_bb_nobsnopcr(ledges, redges, counts, widths, rates, amplitude_criteria=3, minflu=5, manual_groups=manual_groups, split_times=splits)#get_flare_bb_nobsnopcr(ledges_p, redges_p, counts_p, widths_p, rates_p, amp_crit, countmin)
    else:
        dp = d
        bp = b
        lp = l
    
    print("Quiescent Count Rate (10**-3 ct/s):", np.around((lp[0][0]/86400.0)/(10**(-3)), 3), "+/-", np.around((lp[0][1]/86400.0)/(10**(-3)), 3), "\n", file=o)
   # print("Quiescent Count Rate (10^-3 ct/s):", np.around((269.678/86400.0)/(10**(-3)), 3), "+/-", np.around((131.224/86400.0)/(10**(-3)), 3), "\n", file=o)
    
    quies_rate_cts = lp[0][0] / 86400.0

    #to find the number of flares: 
    if type(b) is type(None):
        print("No Flares in Observation \n", file=o)
      
    else:
        
        b = np.asarray(np.transpose(b))
        bp = np.asarray(np.transpose(bp))
        num_flare = b[0].size
    
        #for each flare, get the information 
        for i in range (0,num_flare):
            start = bp[0][i] #start time of the flare
            end = bp[1][i]  #end time of the flare 
        
            '''Picking out the mximum count rate from lightcurve:  '''
            #within this time region, need to find the max count rate:
            upper = np.where(time_lc <= bp[1][i])
            lower = np.where(time_lc >= bp[0][i]) #indices that say where 
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
            ct_mean = bp[4][i]/86400.0 #count rate in each block

            ct_mean_err = bp[5][i]/86400.0 #err in mean ct rate
        
            ct_quies_sub = ct_mean - (lp[0][0]/86400.0) #quiescence subtracted countrates
            dur = bp[3][i]*86400.0 #duration of flare in s 
        
            if gratingtype == 'False':
                #Determining values
                luminosity = (ct_quies_sub*(10**(34)/0.0118766)) # erg/s  
                energy = luminosity * dur  #erg
                #Error propagation
                lum_err = (lp[0][1]/86400)*(10**(34)/0.0118766)
                energy_err = dur * lum_err
            elif gratingtype == 'True':
                #Determining values
                luminosity = (ct_quies_sub*(10**(34)/0.008485)) # erg/s  
                energy = luminosity * dur  #erg
                #Error propagation
                lum_err = (lp[0][1]/86400)*(10**(34)/0.008485)
                energy_err = dur * lum_err

            
            
        
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
            fluence = ct_mean*bp[3][i]*86400.0
        
        
        
            print("-------------------------------------------------------------------------------------", "\n", file=o)
            print("FLARE NUMBER ", i+1,"\n", file=o)
            print("-------------------------------------------------------------------------------------", "\n", file=o)
            print("Start Time: ", start, "(MJD) \n", file=o)
            print("End Time: ", end, "(MJD) \n", file=o)
            print("Duration: ", dur, "(s) \n", file=o)
            print("Count Rate (mean): ",np.around(ct_mean,4) , "+/-", np.around(ct_mean_err, decimals=4), "(ct/s) \n", file=o)
            print("Count Rate (max): ",np.around(ct_max,4) ,"+/-", np.around(ct_max_err,4), "(ct/s) \n", file=o)
            print("Energy: ", np.around(energy/10**(37), 4), "+/-", np.around(energy_err/10**(37),4), "10**37 ergs \n", file =o)
            print("Luminosity: ", np.around(luminosity/10**(34),4), "+/-", np.around(lum_err/10**(34),4), "10**34 erg/s \n", file=o)
            print("Flux: ", np.around(flux/(10**(-12)),4) ,"+/-" ,np.around(flux_err/(10**(-12)),4) , "10**-12 erg/s/cm2 \n", file=o)
            print("Fluence: ", np.around(fluence,4), "ct \n", file=o)

            # ------------------------------------------------------------------
            # Times within this flare where LC falls below 3 * quiescence
            # ------------------------------------------------------------------
            in_flare = (time_lc >= start) & (time_lc <= end)

            # boolean mask for points inside the flare that are below 3 * quies
            low_mask = in_flare & (count_rate_lc <= 3.0 * quies_rate_cts)

            times_below_3q = time_lc[low_mask]          # MJD times
            rates_below_3q = count_rate_lc[low_mask]    # corresponding ct/s

            # Print to the output file (MJD + rate); adjust formatting as you like
            if times_below_3q.size > 0:
                print("Times in flaring block below 3*quiescence (MJD, ct/s):", file=o)
                for t_b, r_b in zip(times_below_3q, rates_below_3q):
                    print(f"  {t_b:.10f}  {r_b:.5e}", file=o)
                print("", file=o)  # blank line for readability
            else:
                print("No points in flaring block below 3*quiescence within this flare.\n", file=o)
        
    f.close()
    o.close()
    




'''___________________________________________________________________
'''


def get_flare_bb_nobsnopcr(ledges, redges, counts, widths, rates,
                           amplitude_criteria=3, minflu=7,
                           manual_groups=None, split_times=None):
    '''
    Outputs flare parameters for a given set of parameters from a 
    Bayesian Blocks analysis without bootstrap. For subarray mode 
    (not gratings!). 
    
    This function is taken from Eli Bouffard's code and lightly
    modified to allow:
      - manual grouping of blocks into flares (manual_groups)
      - manual splitting of blocks at specified times (split_times)
    
    Parameters
    ----------
    ledges : array-like of float
        Time of the beginning of each block (s)
        
    redges : array-like of float
        Time of the end of each block (s)
        
    counts : array-like of int or float
        Total counts in each block 
         
    widths : array-like of float
        Total length of each block (s)
         
    rates : array-like of float
        Mean count rate of each block (ct/s)

    amplitude_criteria : float, optional
        Sigma range above quiescence for a block to be considered a 
        flare. Default is 4.
        
    minflu : int, optional
        Minimum number of counts in a block. If lower, it is combined 
        with a nearby block. Default is 7.
        
    manual_groups : list of lists of int, optional
        If given, overrides the automatic flare grouping.
        Each inner list is a set of block indices (after low-count
        merging AND any splitting) that form one flare, e.g.:
            [[3, 4, 5], [6, 7, 8]]
        means:
            flare 1 = blocks 3,4,5
            flare 2 = blocks 6,7,8

    split_times : list of float, optional
        Times (in same units as ledges/redges) at which to split blocks.
        For each t_split, the block that contains t_split is divided
        into two blocks at t_split, with counts split in proportion to
        the time width (so the rate stays the same).
        This happens after low-count merging and before flare grouping.

    Returns
    -------
    data : 2D ndarray (Nblocks, 6)
        For each block (after low-count merging and optional splitting):
        [t_start, t_end, counts, width, rate, rate_error]
        
    block : 2D ndarray (Nflares, 6) or None
        Same columns as `data` but for each FLARE instead of each block
        (flares can be made of multiple blocks). If no flare is found,
        this can be None.
        
    LoRate : 2D ndarray (1, 2)
        [quiescent_rate, quiescent_rate_error]
    '''
    ledges = np.asarray(ledges, dtype=float)
    redges = np.asarray(redges, dtype=float)
    counts = np.asarray(counts, dtype=float)
    widths = np.asarray(widths, dtype=float)
    rates = np.asarray(rates, dtype=float)

    num_blocks = np.size(redges)
    l = 0  # counts the number of flares
    rateserr = np.sqrt(counts) / widths
    block = None

    # ------------------------------------------------------------------
    # If there are more than 1 block, we can look for flares
    # ------------------------------------------------------------------
    if num_blocks > 1:
        # Merge low-fluence blocks (counts < minflu) into neighbours
        del_blocks = np.array([], dtype=int)

        if (counts < minflu).any():
            lowfluence = np.where(counts < minflu)[0]

            for i in lowfluence:
                if i == 0:
                    # merge first block into next
                    counts[i + 1] += counts[i]
                    ledges[i + 1] = ledges[i]
                    widths[i + 1] += widths[i]
                    rates[i + 1] = counts[i + 1] / widths[i + 1]
                    rateserr[i + 1] = np.sqrt(counts[i + 1]) / widths[i + 1]
                    del_blocks = np.append(del_blocks, i)
                else:
                    # merge into previous block
                    counts[i - 1] += counts[i]
                    redges[i - 1] = redges[i]
                    widths[i - 1] += widths[i]
                    rates[i - 1] = counts[i - 1] / widths[i - 1]
                    rateserr[i - 1] = np.sqrt(counts[i - 1]) / widths[i - 1]
                    del_blocks = np.append(del_blocks, i)

            # delete merged-away blocks
            keep = np.ones(len(counts), dtype=bool)
            keep[del_blocks.astype(int)] = False
            counts = counts[keep]
            ledges = ledges[keep]
            redges = redges[keep]
            widths = widths[keep]
            rates = rates[keep]
            rateserr = rateserr[keep]

        # ------------------------------------------------------------------
        # OPTIONAL: split blocks at specified times
        # ------------------------------------------------------------------
        if split_times is not None:
            # ensure list-like
            split_times = np.atleast_1d(split_times)
            # sort by time so splits proceed in temporal order
            split_times = np.sort(split_times)

            for t_split in split_times:
                # find block that contains t_split
                idx = np.where((t_split > ledges) & (t_split < redges))[0]
                if idx.size != 1:
                    # either no block or ambiguous; skip this split
                    continue
                idx = idx[0]

                tL = ledges[idx]
                tR = redges[idx]
                W = widths[idx]
                C = counts[idx]

                # sanity check
                if not (tL < t_split < tR):
                    continue

                W1 = t_split - tL
                W2 = tR - t_split
                if W1 <= 0 or W2 <= 0:
                    continue

                # split counts proportional to width so rate stays constant
                C1 = C * (W1 / W)
                C2 = C - C1

                # new blocks
                ledges_new = []
                redges_new = []
                counts_new = []
                widths_new = []
                rates_new = []
                rateserr_new = []

                for j in range(len(ledges)):
                    if j != idx:
                        ledges_new.append(ledges[j])
                        redges_new.append(redges[j])
                        counts_new.append(counts[j])
                        widths_new.append(widths[j])
                        rates_new.append(rates[j])
                        rateserr_new.append(rateserr[j])
                    else:
                        # insert split block A
                        ledges_new.append(tL)
                        redges_new.append(t_split)
                        counts_new.append(C1)
                        widths_new.append(W1)
                        rates_new.append(C1 / W1)
                        rateserr_new.append(np.sqrt(C1) / W1)
                        # insert split block B
                        ledges_new.append(t_split)
                        redges_new.append(tR)
                        counts_new.append(C2)
                        widths_new.append(W2)
                        rates_new.append(C2 / W2)
                        rateserr_new.append(np.sqrt(C2) / W2)

                ledges = np.array(ledges_new)
                redges = np.array(redges_new)
                counts = np.array(counts_new)
                widths = np.array(widths_new)
                rates = np.array(rates_new)
                rateserr = np.array(rateserr_new)

        # Update block count after merging + splitting
        num_blocks = np.size(redges)

        # Quiescent block is the one with the longest duration
        quies_id = np.argmax(widths)

        # Identify flaring blocks: significantly above quiescent rate
        
        flares_id = (rates - amplitude_criteria * rateserr) > \
                    (rates[quies_id] + amplitude_criteria * rateserr[quies_id])
        flares = ma.array(rates, mask=~flares_id)

        # Build the "data" array: one row per block
        data = np.ones(6).reshape(1, 6)
        for h in range(num_blocks):
            data = np.concatenate(
                (
                    data,
                    np.array(
                        [
                            ledges[h],
                            redges[h],
                            counts[h],
                            widths[h],
                            rates[h],
                            rateserr[h],
                        ]
                    ).reshape(1, 6),
                )
            )
        data = np.delete(data, 0, axis=0)
        LoRate = np.array([rates[quies_id], rateserr[quies_id]]).reshape(1, 2)

        # ------------------------------------------------------------------
        # MANUAL GROUPING OPTION
        # ------------------------------------------------------------------
        if manual_groups is not None and len(manual_groups) > 0:
            block_list = []
            for grp in manual_groups:
                grp = np.asarray(grp, dtype=int)
                if grp.size == 0:
                    continue
                # ensure sorted indices
                grp = np.sort(grp)
                # start and end of flare
                t_start = ledges[grp[0]]
                t_end = redges[grp[-1]]
                # sum quantities over the selected blocks
                csum = counts[grp].sum()
                wsum = widths[grp].sum()
                rate = csum / float(wsum)
                rateerr = math.sqrt(csum) / float(wsum)
                block_list.append([t_start, t_end, csum, wsum, rate, rateerr])

            block = np.array(block_list) if block_list else None

        else:
            # ------------------------------------------------------------------
            # ORIGINAL AUTOMATIC GROUPING (Eli's logic)
            # ------------------------------------------------------------------
            # If the obsid has only one flaring block
            if np.size(flares[~flares.mask]) == 1:
                indice = np.where(flares_id)[0][0]
                if l == 0:
                    block = np.array(
                        [
                            ledges[indice],
                            redges[indice],
                            counts[indice],
                            widths[indice],
                            rates[indice],
                            rateserr[indice],
                        ]
                    ).reshape(1, 6)
                    l += 1
                else:
                    block = np.concatenate(
                        (
                            block,
                            np.array(
                                [
                                    ledges[indice],
                                    redges[indice],
                                    counts[indice],
                                    widths[indice],
                                    rates[indice],
                                    rateserr[indice],
                                ]
                            ).reshape(1, 6),
                        ),
                        axis=0,
                    )
            # If multiple flaring blocks, group contiguous ones into flares
            else:
                k = 0  # used to spot the first flare of the obsid
                j = 0  # index over blocks
                while j < num_blocks:
                    if not flares_id[j]:
                        j += 1
                        continue
                    else:
                        if k == 0:
                            k += 1
                            if j < num_blocks - 1:
                                if not flares_id[j + 1]:
                                    # flare consists of a single block
                                    if l == 0:
                                        block = np.array(
                                            [
                                                ledges[j],
                                                redges[j],
                                                counts[j],
                                                widths[j],
                                                rates[j],
                                                rateserr[j],
                                            ]
                                        ).reshape(1, 6)
                                        l += 1
                                    else:
                                        block = np.concatenate(
                                            (
                                                block,
                                                np.array(
                                                    [
                                                        ledges[j],
                                                        redges[j],
                                                        counts[j],
                                                        widths[j],
                                                        rates[j],
                                                        rateserr[j],
                                                    ]
                                                ).reshape(1, 6),
                                            ),
                                            axis=0,
                                        )
                                    j += 1
                                else:
                                    # multiple contiguous flaring blocks
                                    flare_block = np.ones(6).reshape(1, 6)
                                    while flares_id[j] and j < (num_blocks - 1):
                                        flare_block = np.concatenate(
                                            (
                                                flare_block,
                                                np.array(
                                                    [
                                                        ledges[j],
                                                        redges[j],
                                                        counts[j],
                                                        widths[j],
                                                        rates[j],
                                                        rateserr[j],
                                                    ]
                                                ).reshape(1, 6),
                                            ),
                                            axis=0,
                                        )
                                        j += 1

                                    if flares_id[j] and j == (num_blocks - 1):
                                        flare_block = np.concatenate(
                                            (
                                                flare_block,
                                                np.array(
                                                    [
                                                        ledges[j],
                                                        redges[j],
                                                        counts[j],
                                                        widths[j],
                                                        rates[j],
                                                        rateserr[j],
                                                    ]
                                                ).reshape(1, 6),
                                            ),
                                            axis=0,
                                        )

                                    flare_block = np.delete(flare_block, 0, axis=0)

                                    if l == 0:
                                        l += 1
                                        block = np.array(
                                            [
                                                flare_block[0, 0],
                                                flare_block[-1, 1],
                                                np.sum(flare_block[:, 2]),
                                                np.sum(flare_block[:, 3]),
                                                np.sum(flare_block[:, 2])
                                                / float(np.sum(flare_block[:, 3])),
                                                math.sqrt(np.sum(flare_block[:, 2]))
                                                / float(np.sum(flare_block[:, 3])),
                                            ]
                                        ).reshape(1, 6)
                                    else:
                                        block = np.concatenate(
                                            (
                                                block,
                                                np.array(
                                                    [
                                                        flare_block[0, 0],
                                                        flare_block[-1, 1],
                                                        np.sum(flare_block[:, 2]),
                                                        np.sum(flare_block[:, 3]),
                                                        np.sum(flare_block[:, 2])
                                                        / float(
                                                            np.sum(flare_block[:, 3])
                                                        ),
                                                        math.sqrt(
                                                            np.sum(flare_block[:, 2])
                                                        )
                                                        / float(
                                                            np.sum(flare_block[:, 3])
                                                        ),
                                                    ]
                                                ).reshape(1, 6),
                                            ),
                                            axis=0,
                                        )
                            else:
                                # this is the last block → single-block flare
                                if l == 0:
                                    l += 1
                                    block = np.array(
                                        [
                                            ledges[j],
                                            redges[j],
                                            counts[j],
                                            widths[j],
                                            rates[j],
                                            rateserr[j],
                                        ]
                                    ).reshape(1, 6)
                                else:
                                    block = np.concatenate(
                                        (
                                            block,
                                            np.array(
                                                [
                                                    ledges[j],
                                                    redges[j],
                                                    counts[j],
                                                    widths[j],
                                                    rates[j],
                                                    rateserr[j],
                                                ]
                                            ).reshape(1, 6),
                                        ),
                                        axis=0,
                                    )
                                j += 1
                        else:
                            # not the first flare
                            if j < (num_blocks - 1):
                                if not flares_id[j + 1]:
                                    # single-block flare
                                    block = np.concatenate(
                                        (
                                            block,
                                            np.array(
                                                [
                                                    ledges[j],
                                                    redges[j],
                                                    counts[j],
                                                    widths[j],
                                                    rates[j],
                                                    rateserr[j],
                                                ]
                                            ).reshape(1, 6),
                                        ),
                                        axis=0,
                                    )
                                    j += 1
                                else:
                                    # multiple contiguous flaring blocks
                                    flare_block = np.ones(6).reshape(1, 6)
                                    while flares_id[j] and j < (num_blocks - 1):
                                        flare_block = np.concatenate(
                                            (
                                                flare_block,
                                                np.array(
                                                    [
                                                        ledges[j],
                                                        redges[j],
                                                        counts[j],
                                                        widths[j],
                                                        rates[j],
                                                        rateserr[j],
                                                    ]
                                                ).reshape(1, 6),
                                            ),
                                            axis=0,
                                        )
                                        j += 1

                                    if flares_id[j] and j == (num_blocks - 1):
                                        flare_block = np.concatenate(
                                            (
                                                flare_block,
                                                np.array(
                                                    [
                                                        ledges[j],
                                                        redges[j],
                                                        counts[j],
                                                        widths[j],
                                                        rates[j],
                                                        rateserr[j],
                                                    ]
                                                ).reshape(1, 6),
                                            ),
                                            axis=0,
                                        )

                                    flare_block = np.delete(flare_block, 0, axis=0)

                                    block = np.concatenate(
                                        (
                                            block,
                                            np.array(
                                                [
                                                    flare_block[0, 0],
                                                    flare_block[-1, 1],
                                                    np.sum(flare_block[:, 2]),
                                                    np.sum(flare_block[:, 3]),
                                                    np.sum(flare_block[:, 2])
                                                    / float(
                                                        np.sum(flare_block[:, 3])
                                                    ),
                                                    math.sqrt(
                                                        np.sum(flare_block[:, 2])
                                                    )
                                                    / float(
                                                        np.sum(flare_block[:, 3])
                                                    ),
                                                ]
                                            ).reshape(1, 6),
                                        ),
                                        axis=0,
                                    )
                            else:
                                # last block
                                if not flares_id[j - 1]:
                                    block = np.concatenate(
                                        (
                                            block,
                                            np.array(
                                                [
                                                    ledges[j],
                                                    redges[j],
                                                    counts[j],
                                                    widths[j],
                                                    rates[j],
                                                    rateserr[j],
                                                ]
                                            ).reshape(1, 6),
                                        ),
                                        axis=0,
                                    )
                                j += 1

    else:
        # Only one block in the obsid
        data = np.array([ledges, redges, counts, widths, rates, rateserr]).reshape(
            1, 6
        )
        LoRate = np.array([rates, rateserr]).reshape(1, 2)

    # ----------------------------------------------------------------------
    # Sort arrays in time
    # ----------------------------------------------------------------------
    data = data[np.argsort(data[:, 0])]
    if block is not None:
        block = block[np.argsort(block[:, 0])]

    return data, block, LoRate