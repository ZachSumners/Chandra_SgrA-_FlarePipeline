import os
import subprocess
from pycrates import read_file
import matplotlib.pylab as plt
import numpy as np
from astropy.time import Time

#Runs the second part of MEGA's "Guide to analyzing flares" (the part after the DS9 stuff but before the Bayesian blocks). Prompts will show up in the terminal but you should just be able to hit enter each time to continue.

#Should just have to change Chandra observation ID, chip ID and working directory.

#Change here
observationID = 28229
chipID = 7
repro_wd = f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro'

subprocess.call('punlearn dmextract', shell=True, cwd=repro_wd)
subprocess.call(f'pset dmextract infile="acisf{observationID}_repro_evt2.fits[ccd_id={chipID}, energy=2000:8000,sky=region(sgra_large.reg)][bin time=::300]"', shell=True, cwd=repro_wd)
subprocess.call(f'pset dmextract outfile="{observationID}_sgra_2-8keV_lc300_large.fits"', shell=True, cwd=repro_wd)
subprocess.call(f'pset dmextract bkg="acisf{observationID}_repro_evt2.fits[ccd_id={chipID},sky=region(bkg.reg)]"', shell=True, cwd=repro_wd)
subprocess.call('pset dmextract opt="ltc1"', shell=True, cwd=repro_wd)
subprocess.call('pset dmextract clobber = yes', shell=True, cwd=repro_wd)
subprocess.call('dmextract', shell=True, cwd=repro_wd)

subprocess.call('punlearn dmcopy', shell=True, cwd=repro_wd)
subprocess.call(f'pset dmcopy infile="acisf{observationID}_repro_evt2.fits[EVENTS][sky=region(sgra_large.reg)][energy=2000:8000]"', shell=True, cwd=repro_wd)
subprocess.call(f'pset dmcopy outfile="{observationID}_sgra_2-8keV_evt_large.fits"', shell=True, cwd=repro_wd)
subprocess.call('pset dmcopy clobber = yes', shell=True, cwd=repro_wd)
subprocess.call('pset dmcopy option="all"', shell=True, cwd=repro_wd)
subprocess.call('dmcopy', shell=True, cwd=repro_wd)

#Plots the light curve
tab = read_file(f"{repro_wd}/{observationID}_sgra_2-8keV_lc300_large.fits")

chandra_time = tab.get_column("time").values
dt = ((50814 + (0 + chandra_time)/86400) - 60404)*24
rate = tab.get_column("net_rate").values
erate = tab.get_column("err_rate").values
plt.errorbar(dt, rate, yerr=erate, marker="o", color="red", mfc="black",mec="black", ecolor="grey")
plt.xlabel("Time - hours (UTC)")
plt.ylabel("Net Count Rate (counts/sec)")
plt.title('10" SgrA* Light Curve')
plt.show()

#print(dt[0])

#t = Time(times, format='isot')
#textstr = '\n'.join(())
#ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

#plt.title(f"{observationID}_sgra_2-8keV_lc300.fits")
#plt.show()
