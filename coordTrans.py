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
#from crates_contrib.utils import *

observationID = 15045
#new_hdulist = fits.HDUList()

#with fits.open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/acisf15045_bary_evt2.fits') as ff:
	#new_hdulist.append(ff[0])
	
	#plt.imshow(ff[1].data)
	#plt.show()

	#table = Table(ff[1].data)
	#header = ff[1].header
	
	#selection = table[['time', 'expno', 'ccd_id', 'node_id', 'chipx', 'chipy', 'tdetx', 'tdety', 'detx', 'dety', 'x', 'y', 'pha', 'pha_ro', 'energy', 'pi', 'fltgrade', 'grade']]
	#img_hdu = fits.ImageHDU(data=selection.to_pandas().values, header=header)
	#new_hdulist.append(img_hdu)

#print('here')

#new_hdulist.writeto(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/acisf15045_bary_evt2_img.fits', overwrite=True)
#print('here2')
#new_img = fits.open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/acisf15045_bary_evt2_img.fits')
#print(new_img.info())

tr = SimpleCoordTransform(f'/Users/zachsumners/Desktop/Research/Chandra/Pipeline/{observationID}/repro/acisf15045_bary_evt2_img.fits')
ra, dec = tr.convert('world', 'physical', 266.4168371, -29.0078106)
print(ra, dec)
