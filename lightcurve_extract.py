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

def extract_lightcurve_grating(observationID, repro_wd, fileName):
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	p = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(sgra.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	sgra_ccd_id = result[16]

	#pb = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(bkg.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	#result_bkg = pb.stdout.decode("utf-8")
	#bkg_ccd_id = result_bkg[16]

	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[(energy=2000:8000,sky=region(sgra.reg),tg_m=0) || (energy=2000:8000,sky=region(box.reg),tg_m=-1)  || (energy=2000:8000,sky=region(box.reg),tg_m=1)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300.fits"', shell=True, cwd=repro_wd)
	#subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)
	
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[(energy=2000:8000,sky=region(sgra.reg),tg_m=0)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300_order0.fits"', shell=True, cwd=repro_wd)
	#subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)
	
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[(energy=2000:8000,sky=region(box.reg),tg_m=-1)  || (energy=2000:8000,sky=region(box.reg),tg_m=1)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300_order1.fits"', shell=True, cwd=repro_wd)
	#subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)

	#Copies events used in light curve to new file.
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][(sky=region(sgra.reg), tg_m=0) || (sky=region(box.reg), tg_m=1) || (sky=region(box.reg), tg_m=-1)][energy=2000:8000]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_2-8keV_evt.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)
	
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg), tg_m=0][energy=2000:8000]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_2-8keV_evt_order0.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)
	
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][(sky=region(box.reg), tg_m=1) || (sky=region(box.reg), tg_m=-1)][energy=2000:8000]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_2-8keV_evt_order1.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)




	
def extract_lightcurve_magnetar(observationID, repro_wd, fileName):
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	p = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(sgra.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	sgra_ccd_id = 7#result[16]

	pb = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(bkg.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result_bkg = pb.stdout.decode("utf-8")
	bkg_ccd_id = 7#result_bkg[16]

	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[energy=2000:8000,sky=region(sgra.reg)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)
	
	
	
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[energy=2000:8000,sky=region(mag.reg)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300_magnetar.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)
	
	
	
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[energy=2000:8000,sky=region(contam.reg)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300_contam.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)



	#Copies events used in light curve to new file.
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy=2000:8000]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_2-8keV_evt.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)	





def extract_lightcurve(observationID, repro_wd, fileName):
	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	p = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(sgra.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	sgra_ccd_id = result[16]

	pb = subprocess.run(f'dmstat "acisf{observationID}_{fileName}_evt2.fits[sky=region(bkg.reg)][cols ccd_id]"', shell=True, cwd=repro_wd, capture_output=True)
	result_bkg = pb.stdout.decode("utf-8")
	bkg_ccd_id = result_bkg[16]

	subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract infile="acisf{observationID}_{fileName}_evt2.fits[energy=2000:8000,sky=region(sgra.reg)][bin time=::300]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300.fits"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmextract bkg="acisf{observationID}_{fileName}_evt2.fits[ccd_id={bkg_ccd_id},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('dmextract', shell=True, cwd=repro_wd)

	#Copies events used in light curve to new file.
	subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy infile="acisf{observationID}_{fileName}_evt2.fits[EVENTS][sky=region(sgra.reg)][energy=2000:8000]"', shell=True, cwd=repro_wd)
	subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_2-8keV_evt.fits"', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
	subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
	subprocess.call('dmcopy', shell=True, cwd=repro_wd)
