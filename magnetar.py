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
from crates_contrib.utils import *

'''def magnetar_extraction(observationID, repro_wd, fileName):
	src_list = []
	f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/src.reg', 'r')
	for line in f:
		src_param = np.array(line[8:-2].split(','), dtype=np.float32)
		src_list.append(src_param)
		
	src_list_np = np.array(src_list)

	hdu_list = fits.open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/{observationID}_repro_2-8keV_img.fits')
	hdu_list.info()

	image_data = hdu_list[0].data

	plt.figure()
	ax = plt.gca()
	norm = simple_norm(image_data, 'log')

	for i in range(len(src_list_np)):
		ellipse = Ellipse(xy=(src_list_np[i][0]-3993.21, src_list_np[i][1]-3991.74), width=2*src_list_np[i][2], height=2*src_list_np[i][3], angle=src_list_np[i][4], edgecolor='r', fc='None', lw=2)
		ax.text(src_list_np[i][0]-3993.21, src_list_np[i][1]-3991.74, str(i))
		ax.add_patch(ellipse)

	ax.imshow(image_data, origin='lower', norm=norm)
	#ax.set_xlim(80,130)
	#ax.set_ylim(80, 130) 
	plt.show(block=False)

	source = int(input('What is the number of the Sgr A* source? (Enter a single integer): '))
	magnetar = int(input('What is the number of the magnetar? (Enter a single integer): '))
	plt.close()

	sgra = src_list_np[source]
	sgra = [4072.0911, 4108.5754, 2, 2, 90]
	sgra[2] = 2.5406504
	sgra[3] = 2.5406504
	
	magnetar = src_list_np[magnetar]
	magnetar = [4068.6601, 4105.0592, 1, 1, 90]
	magnetar_rad = 2.6422764
	
	delta = 5.1829268*np.cos(np.pi/4)
	print(delta)

	sgra_f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/sgra.reg', 'w')
	sgra_f.write(f'ellipse(17:45:40.0409, -29:00:28.118, 1.25", 1.25", 90)')
	#sgra_f.write(f'ellipse({sgra[0]},{sgra[1]},{sgra[2]},{sgra[3]},{sgra[4]})')
	sgra_f.close()
	
	mag_f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/mag.reg', 'w')
	mag_f.write(f'ellipse({magnetar[0]},{magnetar[1]},{2.6422764},{2.6422764},{magnetar[4]})')
	mag_f.close()

	bkg_f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/bkg.reg', 'w')
	bkg_f.write(f'annulus({magnetar[0]},{magnetar[1]},12.703252,20.3252032)')
	bkg_f.close()

	contam_f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/contam.reg', 'w')
	contam_f.write(f'ellipse({magnetar[0] + delta},{magnetar[1] - delta},{sgra[2]},{sgra[3]},{sgra[4]})\n')
	contam_f.write(f'ellipse({magnetar[0] - delta},{magnetar[1] + delta},{sgra[2]},{sgra[3]},{sgra[4]})\n')
	contam_f.write(f'ellipse({magnetar[0] - delta},{magnetar[1] - delta},{sgra[2]},{sgra[3]},{sgra[4]})')
	contam_f.close()'''
	
	
def magnetar_extraction2(observationID, repro_wd, fileName):
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', 266.41683708333333333, -29.007810555555556)
	print(sgra_ra_px, sgra_dec_px)
	sgra_rad = 2.5406504
	mag_ra_px, mag_dec_px = tr.convert('world', 'physical', 266.4173708, -29.0082889)
	mag_rad = 2.6422764
	#bkg_ra_px, bkg_dec_px = tr.convert('world', 'physical', 266.4170167, -29.0079722)
	con1_ra_px, con1_dec_px = tr.convert('world', 'physical', 266.4179046, -29.0087694)
	con2_ra_px, con2_dec_px = tr.convert('world', 'physical', 266.4168925, -29.0088226)
	con3_ra_px, con3_dec_px = tr.convert('world', 'physical', 266.4178492, -29.0077551)


	sgra_f = open(f'{repro_wd}/sgra.reg', 'w')
	#sgra_f.write(f'ellipse(17:45:40.0409, -29:00:28.118, 1.25", 1.25", 90)')
	sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	sgra_f.close()
	
	mag_f = open(f'{repro_wd}/mag.reg', 'w')
	mag_f.write(f'ellipse({mag_ra_px},{mag_dec_px},{mag_rad},{mag_rad},{0})')
	mag_f.close()

	bkg_f = open(f'{repro_wd}/bkg.reg', 'w')
	bkg_f.write(f'annulus({mag_ra_px},{mag_dec_px},12.703252,20.3252032)')
	bkg_f.close()

	contam_f = open(f'{repro_wd}/contam.reg', 'w')
	contam_f.write(f'ellipse({con1_ra_px},{con1_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	contam_f.write(f'ellipse({con2_ra_px},{con2_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	contam_f.write(f'ellipse({con3_ra_px},{con3_dec_px},{sgra_rad},{sgra_rad},{0})')
	contam_f.close()
	
	

def magnetar_correction(observationID, repro_wd, fileName):
	sgra = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300.fits')
	magnetar = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_magnetar.fits')
	contam = fits.open(f'{repro_wd}/{observationID}_sgra_2-8keV_lc300_contam.fits')
	
	mean_contam = np.mean(contam[1].data['NET_RATE'])
	mean_mag = np.mean(magnetar[1].data['NET_RATE'])
	leak_frac = mean_contam/mean_mag
	
	return leak_frac, mean_mag
	
def quiescent_correction(observationID, repro_wd, fileName, leak_frac, q_mag):
	table_res = "./" + str(observationID) + "/repro/" + "Results/"  + str(observationID) + "_SGRA_TABLE_RESULTS_pileupcorr.txt"
	with open(table_res, 'r') as f:
		data = f.readlines()
		old_quiescent = str(data[16][35:40]) + 'e-03'
		new_quiescent = float(old_quiescent) - leak_frac/3*q_mag
		print(old_quiescent, new_quiescent)
		data[16] = f'Quiescent Count Rate (10^-3 ct/s): {np.around(new_quiescent/(10**(-3)), 3)} (MAGNETAR CORRECTED)\n'
	
	with open(table_res, 'w') as f:
		f.writelines(data)
	
