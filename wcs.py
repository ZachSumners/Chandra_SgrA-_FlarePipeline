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

def wcs_correct(fp, observationID, repro_wd, fileName):
	subprocess.call(f'fluximage acisf{observationID}_{fileName}_evt2.fits {observationID} bin=1 band=broad', shell=True, cwd=repro_wd)
	#subprocess.call(f'mkpsfmap {observationID}_broad_thresh.img {observationID}_broad_thresh.psfmap energy=2.3 ecf=0.9 mode=h', shell=True, cwd=repro_wd)
	
	#subprocess.call(f'wavdetect infile="{observationID}_broad_thresh.img" psffile="{observationID}_broad_thresh.psfmap" expfile="{observationID}_broad_thresh.expmap" scales="1 2 4 6 8 12 16 24 32" outfile="{observationID}_broad_src.fits" scell={observationID}_broad.cell imagefile={observationID}_broad.recon defnbkg={observationID}_broad.nbkg mode=h clob+', shell=True, cwd=repro_wd)
	
	#subprocess.call(f'ds9 {observationID}_broad_thresh.img -catalog csc', shell=True, cwd=repro_wd)
	subprocess.call(f'cp {fp}/sources_csc.tsv {repro_wd}/sources_csc.tsv', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_match infile="src.fits" refsrcfile="sources_csc.tsv[opt skip=1][cols ra=col1,dec=col2]" outfile="{observationID}_xfm.fits" wcsfile="{observationID}_broad_thresh.img" radius=1 clob+', shell=True, cwd=repro_wd)

	subprocess.call(f'wcs_update infile=pcadf{observationID}_000N001_asol1.fits outfile=pcadf{observationID}_000N001_asol1_corrected.fits transformfile={observationID}_xfm.fits wcsfile={observationID}_broad_thresh.img', shell=True, cwd=repro_wd)

	subprocess.call(f'wcs_update infile=acisf{observationID}_{fileName}_evt2.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'wcs_update infile={observationID}_repro_2-8keV_img.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit acisf{observationID}_{fileName}_evt2.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_repro_2-8keV_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
	
	subprocess.call(f'cp {repro_wd}/{observationID}_broad_thresh.img {repro_wd}/{observationID}_broad_thresh_img.fits', shell=True, cwd=repro_wd)
	
	subprocess.call(f'wcs_update infile={observationID}_broad_thresh_img.fits outfile=', shell=True, cwd=repro_wd)
	subprocess.call(f'dmhedit {observationID}_broad_thresh_img.fits file= op=add key=ASOLFILE value=pcadf{observationID}_000N001_asol1_corrected.fits', shell=True, cwd=repro_wd)
