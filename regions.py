import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
from astropy.visualization import simple_norm
from crates_contrib.utils import *
import subprocess
from astropy.wcs import WCS
	
def regions_search(observationID, repro_wd, src_coords, bkg_coords, fileName):
	'''This function selects which region found by searchsources.py is Sgr A*. This can be done with a manual selection or automatic selection
 	assuming the WCS correction is appropriate.
  
  	It also creates regions that the lightcurve extraction commands will use. In other words, it defines which events will be used for lightcurve
   	extraction.'''
	
	#Open the image created earlier and selects which pixel is the center of Sgr A* based on its literature position.
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	#converts coordinates in degrees to pixel coordinates
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', src_coords[0], src_coords[1])

	#Radius of Sgr A* region in pixels based on the resolution of Chandra.
	sgra_rad = 2.5406504

	#Creates the Sgr A* region.
	sgra_f = open(f'{repro_wd}/sgra.reg', 'w')
	sgra_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	sgra_f.close()

	#Creates the background region.
	bkg_f = open(f'{repro_wd}/bkg.reg', 'w')
	
	#bkg_f.write(f'annulus({sgra_ra_px},{sgra_dec_px},{bkg_coords[0]},{bkg_coords[1]})')
	bkg_f.write(f'annulus({sgra_ra_px},{sgra_dec_px},15.455285,25.650407)')
	bkg_f.close()


def regions_search_manual_select(observationID, repro_wd, erange, bkg_coords, fileName):
	'''
	This function allows the user to visually select which region corresponds to Sgr A* and extracts the Sgr A* region at the center of the 
	CIAO located source, not the absolute literature coordinates. This allows for slight offsets in WCS coordinates.
	'''

	#Open the file with all the source regions CIAO found.
	src_list = []
	f = open(f'{repro_wd}/src.reg', 'r')
	for line in f:
		src_param = np.array(line[8:-2].split(','), dtype=np.float32)
		src_list.append(src_param)
	
	#Store these in an array
	src_list_np = np.array(src_list)

	#Create a smoothed image of the source extraction region for visualization
	subprocess.call(f'csmooth {repro_wd}/{observationID}_{fileName}_{erange[0]}-{erange[1]}keV_cropped.fits outfile={repro_wd}/{observationID}_smoothed_image.fits outsigfile={repro_wd}/{observationID}_smoothed_sig.fits outsclfile={repro_wd}/{observationID}_smoothed_scl.fits sigmin=2 clobber=yes sclmap= conmeth=fft conkerneltype=gauss sclmode=compute sigmax=5 sclmin=INDEF sclmax=INDEF', shell=True, cwd=repro_wd)

	#Open the smoothed fits file.
	hdu_list = fits.open(f'{repro_wd}/{observationID}_smoothed_image.fits')

	#Calculate the WCS coordinates of this image.
	world_coords = WCS(hdu_list[0].header)
	image_data = hdu_list[0].data

	#Create a visualization that has the smoothed image and detected regions (labeled with a number).
	fig = plt.figure(figsize=(8, 6))
	ax = fig.add_subplot(111, projection=world_coords)

	norm = simple_norm(image_data, 'log')

	#Plot all the regions on the image.
	for i in range(len(src_list_np)):
		ellipse = Ellipse(xy=(src_list_np[i][0]-3993.21, src_list_np[i][1]-3991.74), width=2*src_list_np[i][2], height=2*src_list_np[i][3], angle=src_list_np[i][4], edgecolor='r', fc='None', lw=2)
		ax.text(src_list_np[i][0]-3993.21, src_list_np[i][1]-3991.74, str(i))
		ax.add_patch(ellipse)

	#Show the image.
	ax.imshow(image_data, origin='lower', norm=norm) 
	plt.show(block=False)

	#Get the user to select which region corresponds to Sgr A*.
	source = int(input('What is the number of the Sgr A* source? (Enter a single integer): '))
	plt.close()

	#Save the coordinates of the selected region.
	sgra = src_list_np[source]
	sgra[2] = 2.5406504
	sgra[3] = 2.5406504

	#Save the central coordinates into sgra.reg 
	sgra_f = open(f'{repro_wd}/sgra.reg', 'w')
	sgra_f.write(f'ellipse({sgra[0]},{sgra[1]},{sgra[2]},{sgra[3]},{sgra[4]})')
	sgra_f.close()

	#Save the background region with the central coordinates selected as well.
	bkg_f = open(f'{repro_wd}/bkg.reg', 'w')
	#bkg_f.write(f'annulus({sgra[0]},{sgra[1]},28.455285,40.650407)')
	bkg_f.write(f'annulus({sgra[0]},{sgra[1]},{bkg_coords[0]},{bkg_coords[1]})')
	bkg_f.close()




def regions_search_grating(observationID, repro_wd, src_coords, bkg_coords, fileName):
	'''This function defines the zeroth and first order regions for Chandra HETG grating observations, as given by the region slice of the evt2 HDUL fits file.'''

	#Open the image created earlier and selects which pixel is the center of Sgr A* based on its literature position.
	tr = SimpleCoordTransform(f'{repro_wd}/{observationID}_broad_thresh_img.fits')
	#converts coordinates in degrees to pixel coordinates
	sgra_ra_px, sgra_dec_px = tr.convert('world', 'physical', src_coords[0], src_coords[1])


	#Open file that gives the zeroth and first order region from the reprocessing.
	hetg_region_file = fits.open(f'{repro_wd}/acisf{observationID}_tgmask.fits')
	regions = hetg_region_file[1].data

	#Radius of Sgr A* (zeroth order) region.
	sgra_rad = 2.5406504
	
	#Create the zeroth order region.
	order0_f = open(f'{repro_wd}/order0.reg', 'w')
	order0_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})')
	#order0_f.write(f'ellipse({float(regions[0][2])},{float(regions[0][3])},{sgra_rad},{sgra_rad},{0})')
	order0_f.close()

	#Write the first order regions.
	order1_f = open(f'{repro_wd}/order1.reg', 'w')
	order1_f.write(f'rotbox({sgra_ra_px},{sgra_dec_px},{float(regions[1][4][0])},5,{regions[1][5]})\n')
	order1_f.write(f'rotbox({sgra_ra_px},{sgra_dec_px},{float(regions[2][4][0])},5,{regions[2][5]})')
	#order1_f.write(f'rotbox({float(regions[1][2])},{float(regions[1][3])},{float(regions[1][4][0])},5,{regions[1][5]})\n')
	#order1_f.write(f'rotbox({float(regions[2][2])},{float(regions[2][3])},{float(regions[2][4][0])},5,{regions[2][5]})')
	order1_f.close()

	#Write a region file with both orders
	order1and0_f = open(f'{repro_wd}/order1and0.reg', 'w')
	order1and0_f.write(f'ellipse({sgra_ra_px},{sgra_dec_px},{sgra_rad},{sgra_rad},{0})\n')
	order1and0_f.write(f'rotbox({sgra_ra_px},{sgra_dec_px},{float(regions[1][4][0])},5,{regions[1][5]})\n')
	order1and0_f.write(f'rotbox({sgra_ra_px},{sgra_dec_px},{float(regions[2][4][0])},5,{regions[2][5]})')
	#order1and0_f.write(f'ellipse({float(regions[0][2])},{float(regions[0][3])},{sgra_rad},{sgra_rad},{0})\n')
	#order1and0_f.write(f'rotbox({float(regions[1][2])},{float(regions[1][3])},{float(regions[1][4][0])},5,{regions[1][5]})\n')
	#order1and0_f.write(f'rotbox({float(regions[2][2])},{float(regions[2][3])},{float(regions[2][4][0])},5,{regions[2][5]})')
	order1and0_f.close()
