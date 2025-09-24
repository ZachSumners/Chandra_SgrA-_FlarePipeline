import subprocess
import numpy as np
from astropy.io import fits
import glob
import os

def count_sources(regfile_path):
	'''This function counts how many sources were detected'''
	count = 0
	with open(regfile_path) as f:
		for line in f:
			line = line.strip()
			if not line:
				continue
			count += 1
	return count

def wcs_correct(fp, observationID, repro_wd, erange, fileName):
	'''This function attempts to do a WCS precision correction to improve the astrometry of our observations. This is a 
 	method outlined by the ciao documentation, and it just follows their instructions.'''

	#Check if pcadf is 000N or 00N (Chandra why do you do this to me)
	pattern = os.path.join(repro_wd, f"pcadf{observationID}_*_asol1.fits")
	matches = glob.glob(pattern)
	asol_path = matches[0]

	num_sources = count_sources(f'{repro_wd}/src.reg')

	#Match the sources in it with a catalogue.
	if num_sources > 3:
		subprocess.call(f'cp {fp}/csc_catalogue.fits {repro_wd}/csc_catalogue.fits', shell=True, cwd=repro_wd)
		subprocess.call(f'wcs_match infile="src.fits" refsrcfile="csc_catalogue.fits[cols ra=RA,dec=DEC,ra_err=RA_ERR,dec_err=DEC_ERR]" outfile="{observationID}_xfm.fits" wcsfile={observationID}_broad_thresh.img radius=1 clob+', shell=True, cwd=repro_wd)
		
		corrected_pattern = os.path.join(repro_wd, f"pcadf{observationID}_*_asol1_corrected.fits")
		corrected_matches = glob.glob(corrected_pattern)
		if len(corrected_matches) == 0:
			print('Initial WCS correction failed. Attempting again with lower order solution and larger radius. The WCS solution will not be as accurate and you may want to manually place the region centroid.')
			subprocess.call(f'cp {fp}/csc_catalogue.fits {repro_wd}/csc_catalogue.fits', shell=True, cwd=repro_wd)
			subprocess.call(f'wcs_match infile="src.fits" refsrcfile="csc_catalogue.fits[cols ra=RA,dec=DEC,ra_err=RA_ERR,dec_err=DEC_ERR]" outfile="{observationID}_xfm.fits" wcsfile={observationID}_broad_thresh.img radius=2 clob+', shell=True, cwd=repro_wd)
	else:
		subprocess.call(f'cp {fp}/csc_catalogue.fits {repro_wd}/csc_catalogue.fits', shell=True, cwd=repro_wd)
		subprocess.call(f'wcs_match infile="src.fits" refsrcfile="csc_catalogue.fits[cols ra=RA,dec=DEC,ra_err=RA_ERR,dec_err=DEC_ERR]" outfile="{observationID}_xfm.fits" wcsfile={observationID}_broad_thresh.img radius=1 method=trans clob+', shell=True, cwd=repro_wd)

	#With that matching, update the wcs information in various files.
	subprocess.call(f'wcs_update infile={asol_path} outfile=pcadf{observationID}_000N001_asol1_corrected.fits transformfile={observationID}_xfm.fits wcsfile={observationID}_broad_thresh.img', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile=acisf{observationID}_{fileName}_evt2.fits transformfile={observationID}_xfm.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile={observationID}_repro_{erange[0]}-{erange[1]}keV_img.fits transformfile={observationID}_xfm.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit acisf{observationID}_{fileName}_evt2.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_repro_{erange[0]}-{erange[1]}keV_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)

	#Make a copy of the image and update the WCS information in it.
	
	subprocess.call(f'wcs_update infile={observationID}_broad_thresh_img.fits transformfile={observationID}_xfm.fits outfile= clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_broad_thresh_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
