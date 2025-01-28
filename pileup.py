import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
from astropy.visualization import simple_norm
#from pycrates import read_file
from astropy.time import Time
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from astropy.table import Table
import pandas as pd

def pileup_correction_grating(observationID, repro_wd, fileName):
	f = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_order0.fits')
	table = Table(f[1].data)

	f_evt = fits.open(f'{repro_wd}/acisf{observationID}_{fileName}_evt2.fits')
	exptime = float(f_evt[1].header['EXPTIME'])

	count_rate = table['COUNT_RATE']*exptime
	count_error = table['COUNT_RATE_ERR']*exptime
	fd = 1-(((np.exp(count_rate)-1)*np.exp(-count_rate))/count_rate)
	fd_error = (-np.exp(count_rate)/count_rate**2 + np.exp(2*count_rate)/count_rate**2 + np.exp(count_rate)/count_rate - 2*np.exp(2*count_rate)/count_rate)*count_error

	pileup_rate_order0 = count_rate/(1-fd)*1/exptime
	pileup_error_order0 = np.sqrt((1/(1-fd))**2*count_error**2 + (count_rate/(1-fd)**2)**2*fd_error**2)*1/exptime
	
	f_order1 = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_order1.fits')
	table_order1 = Table(f_order1[1].data)
	count_rate_order1 = table_order1['COUNT_RATE']
	count_error_order1 = table_order1['COUNT_RATE_ERR']
	
	pileup_rate = np.nan_to_num(pileup_rate_order0) + count_rate_order1
	pileup_error = np.sqrt(np.nan_to_num(pileup_error_order0)**2 + count_error_order1**2)
	
	pileup_rate_num = np.nan_to_num(pileup_rate)
	pileup_error_num = np.nan_to_num(pileup_error)

	table.add_column(pileup_rate_num, index=21, name='RATE_PILEUP')
	table.add_column(pileup_error_num, index=22, name='PILEUP_ERR')

	hdu = table.write(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup_TABLE.fits', overwrite=True, format='fits')
	ft = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup_TABLE.fits')
	ft[1].header['MJDREF'] = 5.08140000000000E+04

	hdul = fits.HDUList([f[0], ft[1], f[2]])
	hdul.writeto(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup.fits', overwrite=True)
	
	

def pileup_correction(observationID, repro_wd, fileName):
	f = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300.fits')
	table = Table(f[1].data)

	f_evt = fits.open(f'{repro_wd}/acisf{observationID}_{fileName}_evt2.fits')
	exptime = float(f_evt[1].header['EXPTIME'])

	count_rate = table['NET_RATE']*exptime
	count_error = table['ERR_RATE']*exptime
	fd = 1-(((np.exp(count_rate)-1)*np.exp(-count_rate))/count_rate)
	fd_error = (-np.exp(count_rate)/count_rate**2 + np.exp(2*count_rate)/count_rate**2 + np.exp(count_rate)/count_rate - 2*np.exp(2*count_rate)/count_rate)*count_error

	pileup_rate = count_rate/(1-fd)*1/exptime
	pileup_error = np.sqrt((1/(1-fd))**2*count_error**2 + (count_rate/(1-fd)**2)**2*fd_error**2)*1/exptime
	

	table.add_column(pileup_rate, index=21, name='RATE_PILEUP')
	table.add_column(pileup_error, index=22, name='PILEUP_ERR')

	hdu = table.write(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup_TABLE.fits', overwrite=True, format='fits')
	ft = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup_TABLE.fits')
	ft[1].header['MJDREF'] = 5.08140000000000E+04

	hdul = fits.HDUList([f[0], ft[1], f[2]])
	hdul.writeto(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup.fits', overwrite=True)
