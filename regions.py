import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
from astropy.visualization import simple_norm
from pycrates import read_file
from astropy.time import Time
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from astropy.table import Table
import pandas as pd
from crates_contrib.utils import *

def regions_search_grating(observationID, repro_wd, fileName):
	#src_list = []
	#f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/src.reg', 'r')
	#for line in f:
	#	src_param = np.array(line[8:-2].split(','), dtype=np.float32)
	#	src_list.append(src_param)
		
	#src_list_np = np.array(src_list)

	#hdu_list = fits.open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/{observationID}_repro_2-8keV_img.fits')
	#hdu_list.info()

	#image_data = hdu_list[0].data

	#plt.figure()
	#ax = plt.gca()
	#norm = simple_norm(image_data, 'log')

	#for i in range(len(src_list_np)):
	#	ellipse = Ellipse(xy=(src_list_np[i][0]-3993.21, src_list_np[i][1]-3991.74), width=2*src_list_np[i][2], height=2*src_list_np[i][3], angle=src_list_np[i][4], edgecolor='r', fc='None', lw=2)
	#	ax.text(src_list_np[i][0]-3993.21, src_list_np[i][1]-3991.74, str(i))
	#	ax.add_patch(ellipse)

	#ax.imshow(image_data, origin='lower', norm=norm)
	#ax.set_xlim(80,130)
	#ax.set_ylim(80, 130) 
	#plt.show(block=False)

	#source = int(input('What is the number of the Sgr A* source? (Enter a single integer): '))
	#plt.close()

	#sgra = src_list_np[source]
	#sgra[2] = 2.5406504
	#sgra[3] = 2.5406504
	
	tr = SimpleCoordTransform(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', 266.41683708333333333, -29.007810555555556)
	sgra_rad = 2.5406504

	hetg_region_file = fits.open(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/acisf{observationID}_tgmask.fits')
	regions = hetg_region_file[1].data

	sgra_f = open(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/sgra.reg', 'w')
	sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	sgra_f.close()

	bkg_f = open(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/bkg.reg', 'w')
	bkg_f.write(f'annulus({sgra_ra_px},{sgra_dec_px},28.455285,40.650407)')
	bkg_f.close()

	box_f = open(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/box.reg', 'w')
	box_f.write(f'rotbox({sgra_ra_px},{sgra_dec_px},1020,5.0813008,{regions[1][5]})\n')
	box_f.write(f'rotbox({sgra_ra_px},{sgra_dec_px},1020,5.0813008,{regions[2][5]})')
	box_f.close()
	
	
	
	
def regions_search(observationID, repro_wd, fileName):
	'''src_list = []
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
	plt.close()

	sgra = src_list_np[source]
	sgra[2] = 2.5406504
	sgra[3] = 2.5406504

	sgra_f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/sgra.reg', 'w')
	sgra_f.write(f'ellipse({sgra[0]},{sgra[1]},{sgra[2]},{sgra[3]},{sgra[4]})')
	sgra_f.close()

	bkg_f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/bkg.reg', 'w')
	bkg_f.write(f'annulus({sgra[0]},{sgra[1]},28.455285,40.650407)')
	bkg_f.close()'''
	
	
	tr = SimpleCoordTransform(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', 266.41683708333333333, -29.007810555555556)
	sgra_rad = 2.5406504


	sgra_f = open(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/sgra.reg', 'w')
	#sgra_f.write(f'ellipse(17:45:40.0409, -29:00:28.118, 1.25", 1.25", 90)')
	sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	sgra_f.close()

	bkg_f = open(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/bkg.reg', 'w')
	bkg_f.write(f'annulus({sgra_ra_px},{sgra_dec_px},12.703252,20.3252032)')
	bkg_f.close()

