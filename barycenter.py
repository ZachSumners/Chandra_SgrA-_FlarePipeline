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

def barycenter_corr(observationID, repro_wd, fileName):
	coord = subprocess.run(f'dmlist acisf{observationID}_repro_evt2.fits header | egrep "(RA|DEC)_TARG"', shell=True, cwd=repro_wd, capture_output=True)
	coord_result = coord.stdout.decode("utf-8")
	ra = coord_result.split()[2]
	dec = coord_result.split()[11]

	subprocess.call('punlearn axbary', shell=True, cwd=repro_wd)
	subprocess.call(f'pset axbary infile="acisf{observationID}_repro_evt2.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'cp ~/Desktop/Research/Chandra/Pipeline/{observationID}/primary/*eph1.fits ~/Desktop/Research/Chandra/Pipeline/{observationID}/repro/orbit_eph0.fits', shell=True, cwd=repro_wd)
	subprocess.call(f'pset axbary orbitfile="orbit_eph0.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset axbary outfile="acisf{observationID}_{fileName}_evt2.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset axbary ra={ra} dec={dec}', shell=True, cwd=repro_wd)
	subprocess.call(f'axbary', shell=True, cwd=repro_wd)
