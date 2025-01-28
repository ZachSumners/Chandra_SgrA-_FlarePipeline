import subprocess
import os
from astropy.io import fits

from barycenter import barycenter_corr
from wcs import wcs_correct
from regions import regions_search, regions_search_grating
from pileup import pileup_correction, pileup_correction_grating
from plotcurve import plot_lightcurve
from searchsources import find_sources
from lightcurve_extract import extract_lightcurve, extract_lightcurve_grating, extract_lightcurve_magnetar
from countOrders import grating_pileup
from magnetar import magnetar_correction, quiescent_correction, magnetar_extraction2

#Runs MEGA's "Guide to analyzing flares" (plus some other corrections). Prompts will show up in the terminal but you should just be able to hit enter each time to continue.

#Should just have to change Chandra observation ID and working directory. It should be the directory with the primary and secondary folders in it.
#Also need to extract the orbit eph file from the zip in the primary folder before use.

#Change here
observationID = 28230
barycentric = True
wcsCorrect = True
reprocess = True

if observationID > 14703 and observationID < 18731:
	magnetar = True
else:
	magnetar = False

#print(os.getcwd())
cwd = os.getcwd()
wd = f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}'

subprocess.call('punlearn ardlib', shell=True, cwd=wd)
#subprocess.call('gunzip ./primary/*gz ./secondary/*gz')
if reprocess == True:
	subprocess.call('chandra_repro', shell=True, cwd=wd)

subprocess.call('punlearn ardlib', shell=True, cwd=wd)
#Change here
repro_wd = f'{wd}/repro'

#Checks to see if a grating was used. Requires a slightly different treatment if so.
f_evt2 = fits.open(f'{repro_wd}/acisf{observationID}_repro_evt2.fits')
grating = f_evt2[1].header['GRATING'].strip()
if grating != 'NONE':
	grating_check = True
	magnetar = False
else:
	grating_check = False
	
#Barycentric corrrections
if barycentric == True:
	fileName = 'bary'
	barycenter_corr(observationID, repro_wd, fileName)
else:
	fileName='repro'

#Search for sources in the image
find_sources(observationID, repro_wd, fileName)

#WCS correction.
if wcsCorrect == True:
	wcs_correct(observationID, repro_wd, fileName)

#Identify SgrA*, background region and first order region (if using HETG grating).
if grating_check == False and magnetar == False:
	regions_search(observationID, repro_wd, fileName)
elif grating_check == False and magnetar == True:
	magnetar_extraction2(observationID, repro_wd, fileName)
else:
	regions_search_grating(observationID, repro_wd, fileName)

#Finds CCD in use and extracts light curve.
if grating_check == False and magnetar == False:
	extract_lightcurve(observationID, repro_wd, fileName)
elif grating_check == False and magnetar == True:
	extract_lightcurve_magnetar(observationID, repro_wd, fileName)
else:
	extract_lightcurve_grating(observationID, repro_wd, fileName)

#Find fraction of light that leaks into Sgr A* region.
if magnetar == True:
	leak_frac, q_mag = magnetar_correction(observationID, repro_wd, fileName)

#Pileup correction
if grating_check == False:
	pileup_correction(observationID, repro_wd, fileName)
else:
	pileup_correction_grating(observationID, repro_wd, fileName)

#Plots the light curve
plot_lightcurve(observationID, repro_wd, fileName)

#Runs bayesian blocks
subprocess.call(f'python3 RUN.py {observationID} False {grating_check}', shell=True)

#Fix the pileup if uses grating.
if grating_check == True:
	grating_pileup(observationID)

#Correct reported quiescent level if magnetar contaminates.
if magnetar == True:
	quiescent_correction(observationID, repro_wd, fileName, leak_frac, q_mag)

