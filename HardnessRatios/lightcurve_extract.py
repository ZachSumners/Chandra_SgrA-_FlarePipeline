import subprocess
import os
import numpy as np
from astropy.io import fits

def general_lightcurve_extraction(infile, outfile, bkg, repro_wd):
	'''The commands we use to extract a lightcurve are generally the same. This function are those commands, we just need to change file names.'''

	#Extract a general lightcurve given some file names.
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	if bkg != None:
		subprocess.call(f'dmextract infile={infile} outfile={outfile} bkg={bkg} opt="ltc1" clobber = yes', shell=True, cwd=repro_wd)
	else:
		subprocess.call(f'dmextract infile={infile} outfile={outfile} opt="ltc1" clobber = yes', shell=True, cwd=repro_wd)

def extract_lightcurve_magnetar(observationID, repro_wd, erange, tbin, fileName):
	'''This function extracts the light curve from the specific region for observations where the magnetar is present. This requires special treatment
 	because there are 5 regions instead of 2.
  
  	A lightcurve is created for the Sgr A* region, magnetar region, and contamination region individually. We copy the events used in the Sgr A* region
   	to a new fits file for archiving.
	'''

	#Open the events files and prepare for lightcurve extraction. Calculates which events belong to Sgr A* and which to the background.
	#All events need to fall on the same CCD (this is assumed for now).
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	p = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(sgra.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	sgra_ccd_id = 7#result[16]

	#The same for the background.
	pb = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(bkg.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result_bkg = pb.stdout.decode("utf-8")
	bkg_ccd_id = 7#result_bkg[16]

	#Sgr A* lightcurve extraction as given by the Guide to Analyzing Flares.
	#The double quotes in the name ARE necessary because the "" is actually sent to the command line.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(sgra.reg)][bin time=::{tbin}]"', f'"{observationID}_eff_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Magnetar lightcurve extraction as given by the Guide to Analyzing Flares.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(mag.reg)][bin time=::{tbin}]"', f'"{observationID}_magnetar_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Contamination region lightcurve extraction as given by the Guide to Analyzing Flares.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(contam.reg)][bin time=::{tbin}]"', f'"{observationID}_contam_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Copies events used in Sgr A* lightcurve to new file.
	#subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	#subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]"', shell=True, cwd=repro_wd)
	#subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_{erange[0]}-{erange[1]}keV_evt.fits"', shell=True, cwd=repro_wd)
	#subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	#subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	#subprocess.call('dmcopy', shell=True, cwd=repro_wd)	


def extract_lightcurve(observationID, repro_wd, erange, tbin, fileName):
	'''This function extracts a lightcurve for Sgr A* given the regions found in other files. We also copy the events used to create the lightcurve into
 	a new fits files for archiving.'''

	#Open the events files and prepare for lightcurve extraction. Calculates which events belong to Sgr A* and which to the background.
	#All events need to fall on the same CCD (this is assumed for now).
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	p = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(sgra.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	sgra_ccd_id = result[16]

	#Do the same for the background.
	pb = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(bkg.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result_bkg = pb.stdout.decode("utf-8")
	bkg_ccd_id = result_bkg[16]

	#Sgr A* lightcurve extraction as given by the Guide to Analyzing Flares.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(sgra.reg)][bin time=::{tbin}]"', f'"{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)

	#Copies events used in Sgr A* lightcurve to new file.
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="{observationID}_sgra_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes option=all', shell=True, cwd=repro_wd)

def extract_lightcurve_grating(observationID, repro_wd, erange, tbin, fileName):
	'''This function extracts a lightcurve for Sgr A* order 0 and order 1 combined in the HETG observations.'''
		#Open the events files and prepare for lightcurve extraction. Calculates which events belong to Sgr A* and which to the background.
	#All events need to fall on the same CCD (this is assumed for now).
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	p = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(order1and0.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	sgra_ccd_id = result[16]

	#Sgr A* order 0 lightcurve extraction.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(order0.reg)][tg_m=0][bin time=::{tbin}]"', f'"{observationID}_sgra_order0_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"', None, repro_wd)
	#Sgr A* order 1 lightcurve extraction.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(order1.reg)][tg_m=-1,1][bin time=::{tbin}]"', f'"{observationID}_sgra_order1_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"', None, repro_wd)

	# Input lightcurve files
	lc0_file = f"{repro_wd}/{observationID}_sgra_order0_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"
	lc1_file = f"{repro_wd}/{observationID}_sgra_order1_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"
	combined_file = f"{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc{tbin}.fits"

	# Open FITS files
	with fits.open(lc0_file) as lc0, fits.open(lc1_file) as lc1:
		lc0_data = lc0[1].data
		lc1_data = lc1[1].data

		# Check time alignment
		if not np.allclose(lc0_data['TIME'], lc1_data['TIME']):
			raise ValueError("Lightcurves are not time-aligned. Rebin or align them first.")

		# Add count rates and errors
		combined_rate = lc0_data['COUNT_RATE'] + lc1_data['COUNT_RATE']
		combined_error = np.sqrt(lc0_data['COUNT_RATE_ERR']**2 + lc1_data['COUNT_RATE_ERR']**2)

		# Create a new table with the combined data
		new_cols = [
			fits.Column(name='TIME', array=lc0_data['TIME'], format='D'),
			fits.Column(name='COUNTS', array=lc0_data['COUNTS'] + lc1_data['COUNTS'], format='D'),
			fits.Column(name='COUNT_RATE', array=combined_rate, format='E'),
			fits.Column(name='COUNT_RATE_ERR', array=combined_error, format='E'),
		]
		# Create a new table with the combined data and preserve header
		hdu = fits.BinTableHDU.from_columns(new_cols, header=lc0[1].header)
		hdu.name = 'LIGHTCURVE'

		# Copy GTI extension from lc0 or lc1 (they should be identical)
		gti_hdu = fits.BinTableHDU(data=lc0['GTI'].data, header=lc0['GTI'].header)
		gti_hdu.name = 'GTI'

		# Primary HDU (required for valid FITS format)
		primary = fits.PrimaryHDU()

		# Write out all HDUs
		hdul = fits.HDUList([primary, hdu, gti_hdu])
		hdul.writeto(combined_file, overwrite=True)

	#Copies events used in Sgr A* lightcurve to new file.

	# Step 1: Extract tg_m=0 from order0 region
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(order0.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}][tg_m=0]" outfile="{observationID}_sgra_order0_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes option=all', shell=True, cwd=repro_wd)

	# Step 2: Extract tg_m=±1 from order1 region
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(order1.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}][tg_m=-1,1]" outfile="{observationID}_sgra_order1_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes option=all', shell=True, cwd=repro_wd)

	# Step 3: Merge both event files
	subprocess.call('punlearn dmmerge', shell=True, cwd=repro_wd)
	subprocess.call(f'dmmerge infile="{observationID}_sgra_order0_{erange[0]}-{erange[1]}keV_evt.fits, {observationID}_sgra_order1_{erange[0]}-{erange[1]}keV_evt.fits" outfile="{observationID}_sgra_{erange[0]}-{erange[1]}keV_evt_unsorted.fits" clobber=yes', shell=True, cwd=repro_wd)

	# Step 4: Sort merged events by time
	subprocess.call(f'dmsort "{observationID}_sgra_{erange[0]}-{erange[1]}keV_evt_unsorted.fits[EVENTS]" {observationID}_sgra_{erange[0]}-{erange[1]}keV_evt.fits keys=TIME', shell=True, cwd=repro_wd)
