import subprocess
import numpy as np
from astropy.io import fits
from crates_contrib.utils import *
	
def magnetar_extraction2(observationID, repro_wd, fileName):
	'''This function finds the neighbouring regions needed to correct for the magnetar just like we found regions in regions.py
 	How we calculate the contamination is outlined in Bouffard (2019).'''

	#Calculate the center of Sgr A* assuming appropriate WCS corrections from the image.
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', 266.41683708333333333, -29.007810555555556)

	#The radius of the Sgr A* region.
	sgra_rad = 2.5406504

	#The radius of the magnetar contamination regions
	mag_rad = 2.6422764

	#Calculate where the center of the magnetar and contamination regions are in our image.
	mag_ra_px, mag_dec_px = tr.convert('world', 'physical', 266.4173708, -29.0082889)
	con1_ra_px, con1_dec_px = tr.convert('world', 'physical', 266.4179046, -29.0087694)
	con2_ra_px, con2_dec_px = tr.convert('world', 'physical', 266.4168925, -29.0088226)
	con3_ra_px, con3_dec_px = tr.convert('world', 'physical', 266.4178492, -29.0077551)

	#Save the Sgr A* region.
	sgra_f = open(f'{repro_wd}/sgra.reg', 'w')
	sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	sgra_f.close()

	#Save the magnetar region.
	mag_f = open(f'{repro_wd}/mag.reg', 'w')
	mag_f.write(f'ellipse({mag_ra_px},{mag_dec_px},{mag_rad},{mag_rad},{0})')
	mag_f.close()

	#Save the background region
	bkg_f = open(f'{repro_wd}/bkg.reg', 'w')
	bkg_f.write(f'annulus({mag_ra_px},{mag_dec_px},12.703252,20.3252032)')
	bkg_f.close()

	#Save the contamination regions.
	contam_f = open(f'{repro_wd}/contam.reg', 'w')
	contam_f.write(f'ellipse({con1_ra_px},{con1_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	contam_f.write(f'ellipse({con2_ra_px},{con2_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	contam_f.write(f'ellipse({con3_ra_px},{con3_dec_px},{sgra_rad},{sgra_rad},{0})')
	contam_f.close()
	
	

def magnetar_correction(observationID, repro_wd, erange, fileName):
	'''This function calculates how much of the signal in the Sgr A* region is from the magnetar contamination.
 	It does this by analyzing the lightcurves. Also outlined in Bouffard (2019).'''

	#Open the Sgr A*, magnetar and contamination lightcurves.
	sgra = fits.open(f'{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300.fits')
	magnetar = fits.open(f'{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300_magnetar.fits')
	contam = fits.open(f'{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc300_contam.fits')

	#Calculate the fraction of flux that leaks out of the magnetar region into the contamination regions.
	mean_contam = np.mean(contam[1].data['NET_RATE'])
	mean_mag = np.mean(magnetar[1].data['NET_RATE'])
	leak_frac = mean_contam/mean_mag

	#Return this value for later.
	return leak_frac, mean_mag
	
def quiescent_correction(observationID, repro_wd, fileName, leak_frac, q_mag):
	'''This function, much like the pileup correction, corrects the Sgr A* quiescent rate based on the leak fraction.
 	This is reported in the RESULTS.txt file of the bayesian blocks output, so we need to correct that.'''

	#Open the bayesian blocks results summary txt file.
	table_res = "./" + str(observationID) + "/repro/" + "Results/"  + str(observationID) + "_SGRA_TABLE_RESULTS_pileupcorr.txt"
	#Find where it reports the quiescent level.
	with open(table_res, 'r') as f:
		data = f.readlines()
		old_quiescent = str(data[16][35:40]) + 'e-03'
		#Change that level to the corrected rate.
		new_quiescent = float(old_quiescent) - leak_frac/3*q_mag
		print(old_quiescent, new_quiescent)
		data[16] = f'Quiescent Count Rate (10^-3 ct/s): {np.around(new_quiescent/(10**(-3)), 3)} (MAGNETAR CORRECTED)\n'

	#Save this new rate.
	with open(table_res, 'w') as f:
		f.writelines(data)
	
