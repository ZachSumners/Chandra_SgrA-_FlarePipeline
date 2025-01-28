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

def plot_lightcurve(observationID, repro_wd, fileName):
	fig = plt.figure(figsize=(14, 8))
	ax = plt.gca()

	tab = read_file(f"{repro_wd}/{observationID}_sgra_2-8keV_lc300_pileup.fits")

	chandra_time = tab.get_column("time").values
	myFmt = mdates.DateFormatter('%H:%M:%S')

	chandra_mjd = (50814 + (0 + chandra_time)/86400)
	chandra_times = Time(chandra_mjd, format='mjd')
	chandra_datetimes = chandra_times.datetime

	rate = tab.get_column("rate_pileup").values
	erate = tab.get_column("pileup_err").values
	ax.errorbar(chandra_datetimes, rate, yerr=abs(erate), marker="o", color="red", mfc="black",mec="black", ecolor="grey")
	ax.xaxis.set_major_formatter(myFmt)
	ax.set_xlabel("Time - hours (UTC)")
	ax.set_ylabel("Net Count Rate (counts/sec)")

	obsDate = f'{chandra_datetimes[0].month}/{chandra_datetimes[0].day}/{chandra_datetimes[0].year}'

	ax.set_title(f'SgrA* Light Curve - ObsID {observationID} - Obs Start on {obsDate}')
	plt.savefig(f'{repro_wd}/{observationID}_300lc_2-8keV.png')
