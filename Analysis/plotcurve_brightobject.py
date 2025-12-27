import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pycrates import read_file
from astropy.time import Time
import os
import matplotlib.dates as mdates

def plot_lightcurve(observationID, repro_wd, erange, tbin, fileName):
	'''This function plots the lightcurve that was extracted in lightcurve_extract.py'''

	#Create a figure.
	fig = plt.figure(figsize=(14, 8))
	ax = plt.gca()

	#Open the pileup corrected lightcurve file.
	tab = read_file(f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/31020/repro/31020_bright_object_lc.fits")

	#Convert the chandra MJD time to a readable time.
	chandra_time = tab.get_column("time").values
	myFmt = mdates.DateFormatter('%H:%M:%S')
	chandra_mjd = (50814 + (0 + chandra_time)/86400)
	chandra_times = Time(chandra_mjd, format='mjd')
	chandra_datetimes = chandra_times.datetime

	#Get the flux rates and errors.
	rate = tab.get_column("NET_RATE").values
	erate = tab.get_column("ERR_RATE").values

	#Plot everything together.
	ax.errorbar(chandra_datetimes, rate, yerr=abs(erate), marker="o", color="red", mfc="black",mec="black", ecolor="grey")
	ax.xaxis.set_major_formatter(myFmt)
	ax.set_xlabel("Time - hours (UTC)")
	ax.set_ylabel("Net Count Rate (counts/sec)")

	#Get the date of the observation.
	obsDate = f'{chandra_datetimes[0].month}/{chandra_datetimes[0].day}/{chandra_datetimes[0].year}'

	#Save this plot.
	ax.set_title(f'Bright Object Light Curve - ObsID {observationID} - Obs Start on {obsDate}')
	plt.savefig(f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/31020/repro/{observationID}_bright_object.png')

fp = os.getcwd()
#Change your working directory to the observation subfolder.
wd = f'{fp}/31020'
repro_wd = f'{wd}/repro'
plot_lightcurve('31020', repro_wd, [2,8], 300, 'bary')