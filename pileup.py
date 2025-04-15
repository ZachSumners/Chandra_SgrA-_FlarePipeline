import subprocess
import numpy as np
from astropy.io import fits
from astropy.table import Table	
from astropy.io.fits import BinTableHDU

def pileup_calc(table, f, exptime):
	'''This function does the actual calculation of the pileup rate.'''
	#Convert count rates to counts per frame.
	count_rate = table['NET_RATE']*exptime
	count_error = table['ERR_RATE']*exptime

	#Additional factor so no division by 0.
	dividezero_epsilon = 1e-8

	#Apply the pileup correction to determine fraction of counts lost. Calculate the error too.
	fd = 1-(((np.exp(count_rate)-1)*np.exp(-count_rate))/(count_rate + dividezero_epsilon))
	fd_error = (-np.exp(count_rate)/(count_rate + dividezero_epsilon)**2 + np.exp(2*count_rate)/(count_rate + dividezero_epsilon)**2 + np.exp(count_rate)/(count_rate + dividezero_epsilon) - 2*np.exp(2*count_rate)/(count_rate + dividezero_epsilon))*(count_error + dividezero_epsilon)

	#Reverse the effect of fd and put back in binned rates. Calculate the error as well.
	pileup_rate = count_rate/(1-(fd + dividezero_epsilon))*1/exptime
	pileup_error = np.sqrt((1/(1-(fd+dividezero_epsilon)))**2*count_error**2 + (count_rate/(1-(fd+dividezero_epsilon))**2)**2*fd_error**2)*1/exptime
	
	#Add the new pileup rate and pileup error to the lightcurve fits file.
	table.add_column(pileup_rate, index=21, name='RATE_PILEUP')
	table.add_column(pileup_error, index=22, name='PILEUP_ERR')

	#Save this data into a new fits file.
	# Create new HDU with updated table
	hdu1 = BinTableHDU(data=table, header=f[1].header)
	hdu1.header['MJDREF'] = 5.08140000000000E+04  # Add/update MJDREF

	return hdu1





def pileup_correction(observationID, repro_wd, erange, tbin, fileName):
	'''This file does a pileup correction based on the work of Bouffard (2019) (Equation 2).
 	Note that we do the pileup correction on the *lightcurves*. This will save a new lightcurve with the pileup correction applied.'''

	#Open the lightcurve.
	f = fits.open(f'{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc{tbin}.fits')
	table = Table(f[1].data)

	#Open the event files to get the exposure time of the observation.
	f_evt = fits.open(f'{repro_wd}/acisf{observationID}_{fileName}_evt2.fits')
	exptime = float(f_evt[1].header['EXPTIME'])

	hdu1 = pileup_calc(table, f, exptime)

	# Assemble final file with original PRIMARY and GTI extensions
	hdul = fits.HDUList([f[0], hdu1, f[2]])
	hdul.writeto(f'{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits', overwrite=True)




def pileup_correction_magnetar(observationID, repro_wd, erange, tbin, fileName):
	'''This file does a pileup correction based on the work of Bouffard (2019) (Equation 2) on the MAGNETAR (just a different file name).
 	Note that we do the pileup correction on the *lightcurves*. This will save a new lightcurve with the pileup correction applied.'''

	#Open the lightcurve.
	f = fits.open(f'{repro_wd}/{observationID}_magnetar_{erange[0]}-{erange[1]}keV_lc{tbin}.fits')
	table = Table(f[1].data)

	#Open the event files to get the exposure time of the observation.
	f_evt = fits.open(f'{repro_wd}/acisf{observationID}_magnetar_{erange[0]}-{erange[1]}keV_evt.fits')
	exptime = float(f_evt[1].header['EXPTIME'])

	hdu1 = pileup_calc(table, f, exptime)

	# Assemble final file with original PRIMARY and GTI extensions
	hdul = fits.HDUList([f[0], hdu1, f[2]])
	hdul.writeto(f'{repro_wd}/{observationID}_magnetar_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits', overwrite=True)


def pileup_correction_contam(observationID, repro_wd, erange, tbin, fileName):
	'''This file does a pileup correction based on the work of Bouffard (2019) (Equation 2) on the CONTAMINATION REGION (just a different file name).
 	Note that we do the pileup correction on the *lightcurves*. This will save a new lightcurve with the pileup correction applied.'''

	#Open the lightcurve.
	f = fits.open(f'{repro_wd}/{observationID}_contam_{erange[0]}-{erange[1]}keV_lc{tbin}.fits')
	table = Table(f[1].data)

	#Open the event files to get the exposure time of the observation.
	f_evt = fits.open(f'{repro_wd}/acisf{observationID}_contam_{erange[0]}-{erange[1]}keV_evt.fits')
	exptime = float(f_evt[1].header['EXPTIME'])

	hdu1 = pileup_calc(table, f, exptime)

	# Assemble final file with original PRIMARY and GTI extensions
	hdul = fits.HDUList([f[0], hdu1, f[2]])
	hdul.writeto(f'{repro_wd}/{observationID}_contam_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits', overwrite=True)



def pileup_correction_eff(observationID, repro_wd, erange, tbin, fileName):
	'''This file does a pileup correction based on the work of Bouffard (2019) (Equation 2) on the ORIGINAL SGR A* REGION - NOT MAGNETAR SUBTRACTED (just a different file name).
 	Note that we do the pileup correction on the *lightcurves*. This will save a new lightcurve with the pileup correction applied.'''

	#Open the lightcurve.
	f = fits.open(f'{repro_wd}/{observationID}_eff_{erange[0]}-{erange[1]}keV_lc{tbin}.fits')
	table = Table(f[1].data)

	#Open the event files to get the exposure time of the observation.
	f_evt = fits.open(f'{repro_wd}/acisf{observationID}_eff_{erange[0]}-{erange[1]}keV_evt.fits')
	exptime = float(f_evt[1].header['EXPTIME'])

	hdu1 = pileup_calc(table, f, exptime)

	# Assemble final file with original PRIMARY and GTI extensions
	hdul = fits.HDUList([f[0], hdu1, f[2]])
	hdul.writeto(f'{repro_wd}/{observationID}_eff_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits', overwrite=True)

