import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import CompleteBB as bb
from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt
import sys
import numpy as np
import numpy.ma as ma
import math

def order0_correct(bb_info, table_res, evt):
	ledges, redges, counts, widths, rates, bsrstds = np.transpose(np.loadtxt(bb_info))

	f = open(table_res, mode='r')
	lines = f.read().split('\n')
	start = float(lines[24].split()[2])
	end = float(lines[26].split()[2])

	start_block = np.argwhere(ledges == start)[0][0]
	end_block = np.argwhere(redges == end)[0][0]

	hdu_list = fits.open(evt)
	exptime = float(hdu_list[1].header['exptime'])
	table = Table(hdu_list[1].data)

	rates_order0 = []
	rates_order1 = []
	counts_order0 = []
	counts_order1 = []
	#for i in range(end_block - start_block + 1):

	for i in range(len(rates)):
		block_start_idx = table['time'] > ((ledges[i] - 50814)*86400)
		block_end_idx = table['time'] < ((redges[i] - 50814)*86400)
		
		
		
		block_mask = block_start_idx * block_end_idx
		
		
		block_events = table[block_mask]
		count_order0 = sum(block_events['tg_m'] == 0)
		count_order1 = sum(block_events['tg_m'] == -1) + sum(block_events['tg_m'] == 1)
		
		rates_order0.append(count_order0/widths[i])
		rates_order1.append(count_order1/widths[i])
		counts_order0.append(count_order0)
		counts_order1.append(count_order1)

	def unpile(rates_order0):
		rate = (np.array(rates_order0))/86400 * exptime
		fd = 1-(((np.exp(rate)-1)*np.exp(-rate))/rate)
		pileup_rate = (rate/(1-fd)) * 86400 * 1/exptime
		return pileup_rate


	pileup_rate = unpile(rates_order0)
	combined_rate = pileup_rate + rates_order1

	pileup_counts = unpile(counts_order0)
	combined_counts = pileup_counts + counts_order1

	return ledges, redges, combined_counts, widths, combined_rate, bsrstds

def process(infile, outfile, ledges, redges, combined_counts, widths, combined_rate, bsrstds, p0=0.05):
    
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

    ccdid = 7#eventhdu.data.ccd_id.min ()
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
    print(combined_counts)
    print('# p0 = %g' % p0, file=o)
    print('# timesys =', timesys, file=o)
    print('# tstarts =', ' '.join ('%.5f' % t for t in tstarts), file=o)
    print('#tstops  =', ' '.join ('%.5f' % t for t in tstops), file=o)
    print('# n = %d' % times.size, file=o)
    for i in range (len(ledges)):
        s = '%.5f %.5f %.5f %.5f %.5f %.5f' % (ledges[i], redges[i],
                                        combined_counts[i], widths[i],
                                        combined_rate[i], bsrstds[i])
        
        print(s, file=o)
    o.close()


def grating_pileup(observationID):
	obsid = observationID

	bb_info = "./"  + str(obsid) + "/repro/" + "Results/" + str(obsid) + "_sgra_bayesianBlocks_info.txt" #block info
	table_res = "./" + str(obsid) + "/repro/" + "Results/" + str(obsid) + "_SGRA_TABLE_RESULTS.txt"

	evt = "./" + str(obsid) + "/repro/" + (str(obsid) +  "_sgra_2-8keV_evt.fits") 
	bb_info_new = "./"  + str(obsid) + "/repro/" + "Results/combined_" + str(obsid) + "_sgra_bayesianBlocks_info_pileupcorr.txt" #block info 
	lc = "./" +  str(obsid) + "/repro/" + (str(obsid) + "_sgra_2-8keV_lc300_pileup.fits")
	plot = "./" + str(obsid) + "/repro/" + "Results/combined_" + str(obsid) + "_PLOT_sgra_pileupcorr.png" #plot 
	table_res_new = "./" + str(obsid) + "/repro/" + "Results/combined_" + str(obsid) + "_SGRA_TABLE_RESULTS_pileupcorr.txt" #info for flare table

	rate_header = 'RATE_PILEUP'
	rate_err_header = 'PILEUP_ERR'
	
	#Take block info from unpiled BB run, perform pileup correction.
	ledges, redges, combined_counts, widths, combined_rate, bsrstds = order0_correct(bb_info, table_res, evt)
	#Create a new BB info file.
	process(evt, bb_info_new, ledges, redges, combined_counts, widths, combined_rate, bsrstds)

	#Plot new Bayesian blocks alongside the pileup corrected lightcurve from pileup.py
	fig = plt.figure()
	bb.plot_bb(bb_info_new) 
	bb.plot_lc(lc, rate_header, rate_err_header) 
	plt.xlabel("Time (days)")
	plt.ylabel("Count Rate")
	plt.title("ObsID " + str(obsid) + " - SgrA")
	fig.savefig(plot)

	#Create a new BB flare info file. 
	bb.getInfo(evt , lc , bb_info_new, table_res_new, rate_header, rate_err_header)
