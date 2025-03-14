import subprocess
import numpy as np
from astropy.io import fits

def wcs_correct(fp, observationID, repro_wd, fileName):
	'''This function attempts to do a WCS precision correction to improve the astrometry of our observations. This is a 
 	method outlined by the ciao documentation, and it just follows their instructions.'''

	#Create new image and match the sources in it with a catalogue.
	subprocess.call(f'fluximage acisf{observationID}_{fileName}_evt2.fits {observationID} bin=1 band=broad', shell=True, cwd=repro_wd)
	subprocess.call(f'cp {fp}/sources_csc.tsv {repro_wd}/sources_csc.tsv', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_match infile="src.fits" refsrcfile="sources_csc.tsv[opt skip=1][cols ra=col1,dec=col2]" outfile="{observationID}_xfm.fits" wcsfile={observationID}_broad_thresh.img radius=1 clob+', shell=True, cwd=repro_wd)

	#With that matching, update the wcs information in various files.
	subprocess.call(f'wcs_update infile=pcadf{observationID}_000N001_asol1.fits outfile=pcadf{observationID}_000N001_asol1_corrected.fits transformfile={observationID}_xfm.fits wcsfile={observationID}_broad_thresh.img', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile=acisf{observationID}_{fileName}_evt2.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile={observationID}_repro_{erange[0]}-{erange[1]}keV_img.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit acisf{observationID}_{fileName}_evt2.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_repro_{erange[0]}-{erange[1]}keV_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)

	#Make a copy of the image and update the WCS information in it.
	subprocess.call(f'cp {repro_wd}/{observationID}_broad_thresh.img {repro_wd}/{observationID}_broad_thresh_img.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile={observationID}_broad_thresh_img.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_broad_thresh_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
