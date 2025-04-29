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
from regions import regions_search, regions_search_manual_select, regions_search_grating
from pileup import pileup_correction, pileup_correction_magnetar, pileup_correction_eff, pileup_correction_contam, pileup_correction_grating
from plotcurve import plot_lightcurve
from searchsources import find_sources
from lightcurve_extract import extract_lightcurve, extract_lightcurve_magnetar, extract_lightcurve_grating
from magnetar import magnetar_correction, quiescent_correction, magnetar_extraction2

#============================#
#Change the observation ID.
observationID = 14460
#Set directory filepath. Defaults to the current working directory.
fp = os.getcwd()
#Change your working directory to the observation subfolder.
wd = f'{fp}/{observationID}'
#============================#
barycentric = True
wcsCorrect = True
reprocess = True
search = False
#============================#
#defining other input variables
#Chandra energy range: list w/lower and upper limits
erange = [2,8]
#define the time bin size in seconds for lightcurve extraction
tbin = 300
#coordinates of the source in degrees 
src_coords = [266.41683708333333333, -29.007810555555556]
#coordinates for the inner and outer radii of the background selection region in pixels
bkg_coords = [28.455285,40.650407]

#Sgr A* observations need special treatment if the magnetar (SGR J1745-2900) was active. This occured in a time period and since observation ID's are sequential, 
#we can define a range of IDs that need magnetar care.
if observationID >= 14702 and observationID < 18731:
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
	grating_check = True
	if magnetar == True:
		print('Chandra observations with gratings and the magnetar are not supported at this time.')
		sys.exit()
else:
	grating_check = False
	
#We apply a barycenter timing correction. If we do so, the corrected files will have "bary" in the name.
if barycentric == True:
	fileName = 'bary'
	barycenter_corr(wd, observationID, repro_wd, fileName)
else:
	fileName='repro'

find_sources(observationID, repro_wd, erange, fileName)

#We need to create smoothed images of the observations to find the sources and reference catalogues in the WCS correction.
subprocess.call(f'fluximage acisf{observationID}_{fileName}_evt2.fits {observationID} bin=1 band=broad clobber=yes', shell=True, cwd=repro_wd)
subprocess.call(f'cp {repro_wd}/{observationID}_broad_thresh.img {repro_wd}/{observationID}_broad_thresh_img.fits', shell=True, cwd=repro_wd)

#Computes a WCS correction on our observations to improve the precision of our coordinate system. This is important when we define where the Sgr A*
#region should go.
if wcsCorrect == True:
	wcs_correct(fp, observationID, repro_wd, erange, fileName)

if grating_check == True:
	#Identifies zero and first order source regions for HETG grating observations and stores them in .reg files.
	regions_search_grating(observationID, repro_wd, src_coords, bkg_coords, fileName)
elif grating_check == False:
	if search == True:
		#Find all the sources in the image, and store a text file with a best fit ellipse for each one.
		regions_search_manual_select(observationID, repro_wd, erange, bkg_coords, fileName)
	else:
		#This step identifies the Sgr A* source region, defines a background region.
		if magnetar == False:
			regions_search(observationID, repro_wd, src_coords, bkg_coords, fileName)
		elif magnetar == True:
			#Identifies the special regions for magnetar observations (see Bouffard 2019)
			magnetar_extraction2(observationID, repro_wd, erange, src_coords, bkg_coords, fileName)
	

#Finds the CCD in use and extracts a light curve based on the regions we just defined. We need to store the light curve of the zeroth and first order
#regions separately for pileup correction later on.
if grating_check == False:
	if magnetar == False:
		extract_lightcurve(observationID, repro_wd, erange, tbin, fileName)
	elif magnetar == True:
		extract_lightcurve_magnetar(observationID, repro_wd, erange, tbin, fileName)
elif grating_check == True:
	extract_lightcurve_grating(observationID, repro_wd, erange, tbin, fileName)

#This step comptues the pileup correction and scales the lightcurves appropriately. Applies to 3 lightcurves if magnetar is present.
if grating_check == False:
	if magnetar == False:
		pileup_correction(observationID, repro_wd, erange, tbin, fileName)
	elif magnetar == True:
		pileup_correction_magnetar(observationID, repro_wd, erange, tbin, fileName)
		pileup_correction_eff(observationID, repro_wd, erange, tbin, fileName)
		pileup_correction_contam(observationID, repro_wd, erange, tbin, fileName)
elif grating_check == True:
	pileup_correction_grating(observationID, repro_wd, erange, tbin, fileName)

#Plots the light curve
#plot_lightcurve(observationID, repro_wd, erange, tbin, fileName)

if magnetar == True:
	leak_frac, q_mag = magnetar_correction(observationID, repro_wd, erange, tbin, fileName)
else:
	leak_frac = 0

#Runs the bayesian blocks algorithm to determine whether a flare has occured and what parameters that flare has.
subprocess.call(f'python3 RUN.py {observationID} {magnetar} {grating_check} {erange[0]} {erange[1]} {tbin} {leak_frac}', shell=True)

print('The lightcurve pipeline is complete.')