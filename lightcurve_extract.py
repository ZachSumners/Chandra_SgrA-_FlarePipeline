import subprocess
import os
import numpy as np
from astropy.io import fits

def general_lightcurve_extraction(infile, outfile, bkg, repro_wd):
	'''The commands we use to extract a lightcurve are generally the same. This function are those commands, we just need to change file names.'''

	#Extract a general lightcurve given some file names.
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile={infile}', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile={outfile}', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract bkg={bkg}', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)

def extract_lightcurve_magnetar(observationID, repro_wd, erange, fileName):
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
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(sgra.reg)][bin time=::300]"', f'"{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Magnetar lightcurve extraction as given by the Guide to Analyzing Flares.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(mag.reg)][bin time=::300]"', f'"{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300_magnetar.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Contamination region lightcurve extraction as given by the Guide to Analyzing Flares.
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(contam.reg)][bin time=::300]"', f'"{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300_contam.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Copies events used in Sgr A* lightcurve to new file.
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_{erange[0]}-{erange[1]}keV_evt.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)	


def extract_lightcurve(observationID, repro_wd, erange, fileName):
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
	general_lightcurve_extraction(f'"acisf{observationID}_{fileName}_evt2.fits[energy={int(erange[0])*1000}:{int(erange[1])*1000},sky=region(sgra.reg)][bin time=::300]"', f'"{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300.fits"', f'"acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', repro_wd)
	
	#Copies events used in Sgr A* lightcurve to new file.
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_{erange[0]}-{erange[1]}keV_evt.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)
