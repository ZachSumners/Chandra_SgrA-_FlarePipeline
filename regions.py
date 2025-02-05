import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
from astropy.visualization import simple_norm
	
def regions_search(observationID, repro_wd, fileName):
	'''This function selects which region found by searchsources.py is Sgr A*. This can be done with a manual selection or automatic selection
 	assuming the WCS correction is appropriate.
  
  	It also creates regions that the lightcurve extraction commands will use. In other words, it defines which events will be used for lightcurve
   	extraction.'''
	
	#Open the image created earlier and selects which pixel is the center of Sgr A* based on its literature position.
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', 266.41683708333333333, -29.007810555555556)

	#Radius of Sgr A* region in pixels based on the resolution of Chandra.
	sgra_rad = 2.5406504

	#Creates the Sgr A* region.
	sgra_f = open(f'{repro_wd}/sgra.reg', 'w')
	sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	sgra_f.close()

	#Creates the background region.
	bkg_f = open(f'{repro_wd}/bkg.reg', 'w')
	bkg_f.write(f'annulus({sgra_ra_px},{sgra_dec_px},12.703252,20.3252032)')
	bkg_f.close()

