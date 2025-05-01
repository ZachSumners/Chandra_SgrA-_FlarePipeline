import subprocess
import numpy as np
from astropy.io import fits
from crates_contrib.utils import *
	
def magnetar_extraction2(observationID, repro_wd, erange, src_coords, bkg_coords, fileName):
	'''This function finds the neighbouring regions needed to correct for the magnetar just like we found regions in regions.py
 	How we calculate the contamination is outlined in Bouffard (2019).'''

	#Calculate the center of Sgr A* assuming appropriate WCS corrections from the image.
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', src_coords[0], src_coords[1])

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
	bkg_f.write(f'annulus({mag_ra_px},{mag_dec_px},{bkg_coords[0]},{bkg_coords[1]})')
	bkg_f.close()

	#Save the contamination regions.
	contam_f = open(f'{repro_wd}/contam.reg', 'w')
	contam_f.write(f'ellipse({con1_ra_px},{con1_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	contam_f.write(f'ellipse({con2_ra_px},{con2_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	contam_f.write(f'ellipse({con3_ra_px},{con3_dec_px},{sgra_rad},{sgra_rad},{0})')
	contam_f.close()

	#Also extract the events from the magnetar region
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(mag.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_magnetar_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_eff_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(contam.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_contam_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	
	
	

def magnetar_correction(observationID, repro_wd, erange, tbin, fileName):
	'''This function calculates how much of the signal in the Sgr A* region is from the magnetar contamination.
 	It does this by analyzing the lightcurves. Also outlined in Bouffard (2019).'''

	#Open the Sgr A*, magnetar and contamination lightcurves.
	eff = fits.open(f'{repro_wd}/{observationID}_eff_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')
	magnetar = fits.open(f'{repro_wd}/{observationID}_magnetar_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')
	contam = fits.open(f'{repro_wd}/{observationID}_contam_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')

	#Calculate the fraction of flux that leaks out of the magnetar region into the contamination regions.
	contamin = contam[1].data['RATE_PILEUP']
	contamin_err = contam[1].data['PILEUP_ERR']

	magn = magnetar[1].data['RATE_PILEUP']
	magn_err = magnetar[1].data['PILEUP_ERR']
	leak_frac = np.nanmean(contamin)/(3*np.nanmean(magn))
	print(f'Leak_frac is {leak_frac}')

	#Calculate the real lightcurve of Sgr A* based on contamination factor from magnetar (Bouffard 2019).
	effective = eff[1].data['RATE_PILEUP']
	print(f'Q_eff is {np.mean(effective)}')
	#Correction
	sgr_lightcurve = effective - leak_frac*magn

	# Open the original file
	with fits.open(f"{repro_wd}/{observationID}_eff_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits", mode='readonly') as hdul:
		# Make a copy of the HDU list in memory
		new_hdul = hdul.copy()

		# Modify a column (example: overwrite column 'FLUX' with random data)
		data = new_hdul[1].data  # Assuming the table is in the first extension (index 1)
		data['RATE_PILEUP'] = sgr_lightcurve  # Replace with your logic

		# Save to new file
		new_hdul.writeto(f"{repro_wd}/{observationID}_sgra_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits", overwrite=True)

	#Return this value for later.
	return np.nanmean(leak_frac), np.nanmean(magn)
	
def quiescent_correction(observationID, repro_wd, fileName, leak_frac, q_mag):
	'''This function, much like the pileup correction, corrects the Sgr A* quiescent rate based on the leak fraction.
 	This is reported in the RESULTS.txt file of the bayesian blocks output, so we need to correct that.'''

	#Open the bayesian blocks results summary txt file.
	table_res = "./" + str(observationID) + "/repro/" + "Results/"  + str(observationID) + "_SGRA_TABLE_RESULTS_pileupcorr.txt"
	#Find where it reports the quiescent level.
	with open(table_res, 'r') as f:
		data = f.readlines()
		old_quiescent = str(data[16][35:40]) + 'e-03'
		print(f'old quiescent rate is {old_quiescent}')
		#Change that level to the corrected rate.
		new_quiescent = float(old_quiescent) - leak_frac*q_mag
		print(old_quiescent, new_quiescent)
		data[16] = f'Quiescent Count Rate (10^-3 ct/s): {np.around(new_quiescent/(10**(-3)), 3)} (MAGNETAR CORRECTED)\n'

	#Save this new rate.
	with open(table_res, 'w') as f:
		f.writelines(data)
	
