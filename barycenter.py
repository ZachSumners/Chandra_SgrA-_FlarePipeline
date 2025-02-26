import subprocess
import os
import glob
from astropy.io import fits

def barycenter_corr(wd, observationID, repro_wd, fileName):
	'''
 	This script applies a barycenter correction to a Chandra observation event file (evt2.fits). Correction parameters are stored in the eph files located in the "primary" folder.

   	The barycenter corrected event file (output of this function) will have "bary" in its name. Further calibration should use that file.

    	Steps of this correction are from the Chandra CIAO documentation (https://cxc.cfa.harvard.edu/ciao/threads/axbary/)
 	'''

	#Find the RA and Dec of the target from the fits header.
	coord = subprocess.run(f'dmlist acisf{observationID}_repro_evt2.fits header | egrep "(RA|DEC)_TARG"', shell=True, cwd=repro_wd, capture_output=True)
	coord_result = coord.stdout.decode("utf-8")
	ra = coord_result.split()[2]
	dec = coord_result.split()[11]

	#Run the barycenter correction on the command line.
	subprocess.call('punlearn axbary', shell=True, cwd=repro_wd)
	#Unzip the eph file if not done already.
	if not glob.glob(f"{wd}/primary/*eph1.fits"):
		subprocess.call(f'gunzip {wd}/primary/*eph1.fits.gz', shell=True, cwd=repro_wd)
	#Copy the eph file originally in the "primary" folder to the repro folder.
	subprocess.call(f'cp {wd}/primary/*eph1.fits orbit_eph0.fits', shell=True, cwd=repro_wd)

	#Run the barycenter correction with the evt2.fits infile and various other files describing the barycenter calculation.
	subprocess.call(f'axbary infile="acisf{observationID}_repro_evt2.fits" orbitfile="orbit_eph0.fits" outfile="acisf{observationID}_{fileName}_evt2.fits" ra={ra} dec={dec}', shell=True, cwd=repro_wd)
