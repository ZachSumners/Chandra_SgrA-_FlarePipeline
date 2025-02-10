'''
This pipeline by Zach Sumners and Nicole Ford was created to minimize the human intervention needed to produce light curves of Chandra Xray observations of Sagittarius A*.

It is inspired by the "Guide to analyzing flares" by McGill University's Extreme Gravity and Accretion group (MEGA) but adds additional treatment of
important precision calibrations such as photon pileup and barycentric timing. 

**To use the pipeline, you need to change:**
- The Chandra observation ID
- Your working directory (this will be the folder that the "primary" and "secondary" folders from a Chandra observation download live in).
- You need to unzip the orbit eph file from the primary folder before using the pipeline.

You can also change which processes you want the pipeline to run.
'''

#Libraries to open the fits files and run command line calls.
import subprocess
import os
import sys
from astropy.io import fits

#Import utility functions from the various pipeline scripts.
from barycenter import barycenter_corr
from wcs import wcs_correct
from regions import regions_search
from pileup import pileup_correction
from plotcurve import plot_lightcurve
from searchsources import find_sources
from lightcurve_extract import extract_lightcurve, extract_lightcurve_magnetar
from magnetar import magnetar_correction, quiescent_correction, magnetar_extraction2

#============================#
#Change the observation ID.
observationID = 28229
#Set directory filepath. Defaults to the current working directory.
fp = os.getcwd()
#Change your working directory to the observation subfolder.
wd = f'{fp}/{observationID}'
#============================#
barycentric = True
wcsCorrect = True
reprocess = True
#============================#

#Sgr A* observations need special treatment if the magnetar (SGR J1745-2900) was active. This occured in a time period and since observation ID's are sequential, 
#we can define a range of IDs that need magnetar care.
if observationID > 14703 and observationID < 18731:
	magnetar = True
else:
	magnetar = False

#The data needs to be reprocessed by the latest Chandra calibration files.
subprocess.call('punlearn ardlib', shell=True, cwd=wd)
if reprocess == True:
	subprocess.call('chandra_repro', shell=True, cwd=wd)

subprocess.call('punlearn ardlib', shell=True, cwd=wd)

#Whether the reprocessed data was just created or already exists, we need to work with the files in the reprocessed folder.
repro_wd = f'{wd}/repro'

#Checks to see if a grating was used. Requires a different treatment if so. We can check this in the fits header.
f_evt2 = fits.open(f'{repro_wd}/acisf{observationID}_repro_evt2.fits')
grating = f_evt2[1].header['GRATING'].strip()
if grating != 'NONE':
	print('Chandra observations with gratings are not supported at this time.')
	sys.exit()
	grating_check = True
	magnetar = False
else:
	grating_check = False
	
#We apply a barycenter timing correction. If we do so, the corrected files will have "bary" in the name.
if barycentric == True:
	fileName = 'bary'
	barycenter_corr(wd, observationID, repro_wd, fileName)
else:
	fileName='repro'

#Find all the sources in the image, and store a text file with a best fit ellipse for each one.
find_sources(observationID, repro_wd, fileName)

#Computes a WCS correction on our observations to improve the precision of our coordinate system. This is important when we define where the Sgr A*
#region should go.
if wcsCorrect == True:
	wcs_correct(fp, observationID, repro_wd, fileName)

#This step identifies the Sgr A* source region, defines a background region and finds the first order region if using a HETG grating.
if grating_check == False and magnetar == False:
	regions_search(observationID, repro_wd, fileName)
elif grating_check == False and magnetar == True:
	magnetar_extraction2(observationID, repro_wd, fileName)

#Finds the CCD in use and extracts a light curve based on the regions we just defined. We need to store the light curve of the zeroth and first order
#regions separately for pileup correction later on.
if grating_check == False and magnetar == False:
	extract_lightcurve(observationID, repro_wd, fileName)
elif grating_check == False and magnetar == True:
	extract_lightcurve_magnetar(observationID, repro_wd, fileName)

#If the magnetar is bright, this step finds fraction of light that leaks into Sgr A* region.
if magnetar == True:
	leak_frac, q_mag = magnetar_correction(observationID, repro_wd, fileName)

#This step comptues the pileup correction and scales the lightcurves appropriately.
if grating_check == False:
	pileup_correction(observationID, repro_wd, fileName)

#Plots the light curve
plot_lightcurve(observationID, repro_wd, fileName)

#Runs the bayesian blocks algorithm to determine whether a flare has occured and what parameters that flare has.
subprocess.call(f'python3 RUN.py {observationID} False {grating_check}', shell=True)

#Similarly, this step scales the quiescent bayesian blocks region if the magnetar has contaminated.
if magnetar == True:
	quiescent_correction(observationID, repro_wd, fileName, leak_frac, q_mag)

print('The lightcurve pipeline is complete.')