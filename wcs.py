import subprocess
import numpy as np
from astropy.io import fits

def wcs_correct(fp, observationID, repro_wd, erange, fileName):
	'''This function attempts to do a WCS precision correction to improve the astrometry of our observations. This is a 
 	method outlined by the ciao documentation, and it just follows their instructions.'''

	#Create new image and match the sources in it with a catalogue.
	
	subprocess.call(f'cp {fp}/csc_catalogue.fits {repro_wd}/csc_catalogue.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_match infile="src.fits" refsrcfile="csc_catalogue.fits[cols ra=RA,dec=DEC,ra_err=RA_ERR,dec_err=DEC_ERR]" outfile="{observationID}_xfm.fits" wcsfile={observationID}_broad_thresh.img radius=1 clob+', shell=True, cwd=repro_wd)

	#With that matching, update the wcs information in various files.
	subprocess.call(f'wcs_update infile=pcadf{observationID}_000N001_asol1.fits outfile=pcadf{observationID}_000N001_asol1_corrected.fits transformfile={observationID}_xfm.fits wcsfile={observationID}_broad_thresh.img', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile=acisf{observationID}_{fileName}_evt2.fits transformfile={observationID}_xfm.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile={observationID}_repro_{erange[0]}-{erange[1]}keV_img.fits transformfile={observationID}_xfm.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit acisf{observationID}_{fileName}_evt2.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_repro_{erange[0]}-{erange[1]}keV_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)

	#Make a copy of the image and update the WCS information in it.
	
	subprocess.call(f'wcs_update infile={observationID}_broad_thresh_img.fits transformfile={observationID}_xfm.fits outfile= clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_broad_thresh_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
