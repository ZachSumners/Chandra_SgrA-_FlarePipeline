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

def find_sources(observationID, repro_wd, fileName):
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][bin x=3992.71:4174.71:1,y=3991.24:4173.24:1][energy=2000:8000]" outfile="{observationID}_{fileName}_2-8keV_cropped.fits" clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call('punlearn mkpsfmap', shell=True, cwd=repro_wd)
	subprocess.call(f'mkpsfmap infile="{observationID}_{fileName}_2-8keV_cropped.fits" outfile="{observationID}_repro_2-8keV_psfmap.fits" energy=3.8 ecf=0.393 clobber=yes', shell=True, cwd=repro_wd)

	#Search for sources in image.
	subprocess.call('punlearn wavdetect', shell=True, cwd=repro_wd)
	subprocess.call(f'pset wavdetect infile="{observationID}_{fileName}_2-8keV_cropped.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset wavdetect psffile="{observationID}_repro_2-8keV_psfmap.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset wavdetect outfile="src.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset wavdetect scellfile="{observationID}_repro_2-8keV_scell.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset wavdetect imagefile="{observationID}_repro_2-8keV_img.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset wavdetect defnbkgfile="{observationID}_repro_2-8keV_nbkg.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset wavdetect regfile="src.reg"', shell=True, cwd=repro_wd)
	subprocess.call('pset wavdetect scales="1.0 2.0 4.0 8.0 16.0"', shell=True, cwd=repro_wd)
	subprocess.call('pset wavdetect sigthresh=1.e-06', shell=True, cwd=repro_wd)
	subprocess.call('pset wavdetect clobber=yes', shell=True, cwd=repro_wd)
	subprocess.call('wavdetect', shell=True, cwd=repro_wd)
