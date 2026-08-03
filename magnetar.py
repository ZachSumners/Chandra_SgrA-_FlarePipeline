import subprocess
import numpy as np
from astropy.io import fits
from crates_contrib.utils import *
import os

def read_ellipse_centers(filename):
    centers = []
    with open(filename, "r") as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("ellipse"):
                # expect something like: ellipse(xc,yc,a,b,theta)
                inside = line[line.find("(")+1 : line.find(")")]
                parts = [p.strip() for p in inside.split(",")]
                if len(parts) < 2:
                    continue
                x, y = float(parts[0]), float(parts[1])
                centers.append((x, y, i, line))
    return centers

def magnetar_extraction2(observationID, repro_wd, erange, src_coords, bkg_coords, fileName):
	'''This function finds the neighbouring regions needed to correct for the magnetar just like we found regions in regions.py
 	How we calculate the contamination is outlined in Bouffard (2019).'''

	#Calculate the center of Sgr A* assuming appropriate WCS corrections from the image.
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', src_coords[0], src_coords[1])

	src_centers = read_ellipse_centers(f'{repro_wd}/src.reg')
	ref = np.array([sgra_ra_px, sgra_dec_px])

	distances = []
	for (x, y, line_no, raw_line) in src_centers:
		d = np.hypot(x - ref[0], y - ref[1])
		distances.append((d, x, y, line_no, raw_line))

	# find closest
	distances.sort(key=lambda t: t[0])

	if distances[0][0] < 2:
		dmin, x_closest, y_closest, line_no, raw_line = distances[0]
		sgra_ra_px = x_closest
		sgra_dec_px = y_closest

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
	if not os.path.exists(f"{repro_wd}/sgra.reg"):
		sgra_f = open(f'{repro_wd}/sgra.reg', 'w')
		sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
		sgra_f.close()

	#Save the magnetar region.
	if not os.path.exists(f"{repro_wd}/mag.reg"):
		mag_f = open(f'{repro_wd}/mag.reg', 'w')
		mag_f.write(f'ellipse({mag_ra_px},{mag_dec_px},{mag_rad},{mag_rad},{0})')
		mag_f.close()

	#Save the background region
	if not os.path.exists(f"{repro_wd}/bkg.reg"):
		bkg_f = open(f'{repro_wd}/bkg.reg', 'w')
		bkg_f.write(f'annulus({mag_ra_px},{mag_dec_px},{bkg_coords[0]},{bkg_coords[1]})')
		#bkg_f.write(f'ellipse(4118.0273,4102.5745,16.217534,16.217534,0)')
		bkg_f.close()

	#Save the contamination regions.
	if not os.path.exists(f"{repro_wd}/contam1.reg"):
		contam_f = open(f'{repro_wd}/contam1.reg', 'w')
		contam_f.write(f'ellipse({con1_ra_px},{con1_dec_px},{sgra_rad},{sgra_rad},{0})\n')
		contam_f.close()

	#Save the contamination regions.
	if not os.path.exists(f"{repro_wd}/contam2.reg"):
		contam_f = open(f'{repro_wd}/contam2.reg', 'w')
		contam_f.write(f'ellipse({con2_ra_px},{con2_dec_px},{sgra_rad},{sgra_rad},{0})\n')
		contam_f.close()

	#Save the contamination regions.
	if not os.path.exists(f"{repro_wd}/contam3.reg"):
		contam_f = open(f'{repro_wd}/contam3.reg', 'w')
		contam_f.write(f'ellipse({con3_ra_px},{con3_dec_px},{sgra_rad},{sgra_rad},{0})')
		contam_f.close()

	#Also extract the events from the magnetar region
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(mag.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_magnetar_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_eff_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(contam1.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_contam1_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(contam2.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_contam2_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(contam3.reg)][energy={int(erange[0])*1000}:{int(erange[1])*1000}]" outfile="acisf{observationID}_contam3_{erange[0]}-{erange[1]}keV_evt.fits" clobber=yes', shell=True, cwd=repro_wd)
	
	
	

def magnetar_correction(observationID, repro_wd, erange, tbin, fileName):
	'''This function calculates how much of the signal in the Sgr A* region is from the magnetar contamination.
 	It does this by analyzing the lightcurves. Also outlined in Bouffard (2019).'''

	#Open the Sgr A*, magnetar and contamination lightcurves.
	eff = fits.open(f'{repro_wd}/{observationID}_eff_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')
	magnetar = fits.open(f'{repro_wd}/{observationID}_magnetar_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')
	contam1 = fits.open(f'{repro_wd}/{observationID}_contam1_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')
	contam2 = fits.open(f'{repro_wd}/{observationID}_contam2_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')
	contam3 = fits.open(f'{repro_wd}/{observationID}_contam3_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits')

	#Calculate the fraction of flux that leaks out of the magnetar region into the contamination regions.
	contamin1 = contam1[1].data['RATE_PILEUP']
	contamin_err1 = contam1[1].data['PILEUP_ERR']
	contamin2 = contam2[1].data['RATE_PILEUP']
	contamin_err2 = contam2[1].data['PILEUP_ERR']
	contamin3 = contam3[1].data['RATE_PILEUP']
	contamin_err3 = contam3[1].data['PILEUP_ERR']

	contamin = (contamin1 + contamin2 + contamin3)/3
	contamin_err = np.sqrt((contamin_err1**2 + contamin_err2**2 + contamin_err3**2)/3)

	magn = magnetar[1].data['RATE_PILEUP']
	magn_err = magnetar[1].data['PILEUP_ERR']

	if np.mean(contamin) < 0.004:
		contamin1 = contam1[1].data['NET_RATE']
		contamin_err1 = contam1[1].data['ERR_RATE']
		contamin2 = contam2[1].data['NET_RATE']
		contamin_err2 = contam2[1].data['ERR_RATE']
		contamin3 = contam3[1].data['NET_RATE']
		contamin_err3 = contam3[1].data['ERR_RATE']

		contamin = (contamin1 + contamin2 + contamin3)/3
		contamin_err = np.sqrt((contamin_err1**2 + contamin_err2**2 + contamin_err3**2)/3)

		
	
	if np.mean(magn) < 0.004:
		magn = magnetar[1].data['NET_RATE']
		magn_err = magnetar[1].data['ERR_RATE']
		
	contamin_mean_err = np.std(contamin_err)/np.sqrt(len(contamin_err))
	
	magn_mean_err = np.std(magn_err)/np.sqrt(len(magn_err))

	leak_frac = np.nanmean(contamin)/(np.nanmean(magn))

	leak_frac_err = np.sqrt((1/np.nanmean(magn))**2 * contamin_mean_err**2 + (np.nanmean(contamin)/np.nanmean(magn)**2)**2 * magn_mean_err)
	print(f'Leak_frac is {leak_frac}')

	#Calculate the real lightcurve of Sgr A* based on contamination factor from magnetar (Bouffard 2019).
	effective = eff[1].data['RATE_PILEUP']
	print(f'Q_eff is {np.mean(effective)}')
	#Correction
	sgr_lightcurve = effective - leak_frac*magn

	with open(f'./{observationID}/repro/MagnetarProperties.txt', 'w') as f:
		f.write(f"Leak fraction is {leak_frac:.6f} +/- {leak_frac_err:.6f}\n")
		f.write(f"Q_eff is {np.mean(effective):.6f}\n")
		f.write(f"Q_contam is {np.mean(contamin):.6f}\n")
		f.write(f'Q_sgra is {np.mean(sgr_lightcurve)}\n')
		f.write(f'Q_mag is {np.nanmean(magn)} +/- {magn_mean_err}')
	
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
	
