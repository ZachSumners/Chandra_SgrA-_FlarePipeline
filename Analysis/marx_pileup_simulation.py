import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pycrates import read_file
from astropy.time import Time
import matplotlib.dates as mdates

wd = os.getcwd()
differences = []
true_fluxes = []

fluxes = []

for i in range(1, 301, 10):
    flux = 0.00001622793002*i
    fluxes.append(flux)

    def plot_lightcurve(observationID, name):
        '''This function plots the lightcurve that was extracted in lightcurve_extract.py'''

        #Create a figure.
        fig = plt.figure(figsize=(14, 8))
        ax = plt.gca()

        #Open the pileup corrected lightcurve file.
        tab = read_file(f"/Users/zachsumners/Desktop/15043/repro/{name}.fits")

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
        ax.set_title(f'SgrA* Light Curve - ObsID {observationID} - Obs Start on {obsDate}')
        plt.savefig(f'/Users/zachsumners/Desktop/15043/repro/{name}.png')



    def difference(name1, name2, flux):
        plt.clf()
        tab1 = read_file(f"/Users/zachsumners/Desktop/15043/repro/{name1}.fits")
        tab2 = read_file(f"/Users/zachsumners/Desktop/15043/repro/{name2}.fits")

        #Convert the chandra MJD time to a readable time.
        chandra_time = tab1.get_column("time").values
        myFmt = mdates.DateFormatter('%H:%M:%S')
        chandra_mjd = (50814 + (0 + chandra_time)/86400)
        chandra_times = Time(chandra_mjd, format='mjd')
        chandra_datetimes = chandra_times.datetime

        #Get the flux rates and errors.
        rate1 = tab1.get_column("NET_RATE").values
        erate1 = tab1.get_column("ERR_RATE").values
        rate2 = tab2.get_column("NET_RATE").values
        erate2 = tab2.get_column("ERR_RATE").values

        difference = rate1 - rate2
        difference_mean = np.mean(difference)
        differences_oneflux.append(difference_mean)

        #true_fluxes.append(np.mean(rate1))
        #observed_fluxes.append(np.mean(rate2))


        plt.plot(chandra_datetimes, difference, marker="o", color="red", mfc="black",mec="black")

        #Get the date of the observation.
        obsDate = f'{chandra_datetimes[0].month}/{chandra_datetimes[0].day}/{chandra_datetimes[0].year}'

        #Save this plot.
        plt.title(f'SgrA* Light Curve - ObsID 15043 Difference - Obs Start on {obsDate}')
        plt.savefig(f'/Users/zachsumners/Desktop/15043/repro/{flux}_DIFFERENCE.png')

    differences_oneflux = []
    for i in range(5):
        subprocess.call(f'pset marx.par SourceFlux={flux}', shell=True, cwd=wd)
        subprocess.call('marx', shell=True, cwd=wd)
        subprocess.call(f'marx2fits point source_{flux}_evt2.fits', shell=True, cwd=wd)
        subprocess.call('marxpileup MarxOutputDir=point', shell=True, cwd=wd)
        subprocess.call(f'marx2fits --pileup point/pileup source_{flux}_pileup_evt2.fits', shell=True, cwd=wd)

        subprocess.call('punlearn dmextract', shell=True, cwd=wd)
        subprocess.call(f'dmextract infile="source_{flux}_evt2.fits[energy=2000:8000,sky=region(sgra.reg)][bin time=::300]" outfile="marx_{flux}_sgra_2-8keV_lc300.fits" bkg="source_{flux}_evt2.fits[ccd_id=7,sky=region(bkg.reg)]" opt="ltc1" clobber=yes', shell=True, cwd=wd)
        subprocess.call(f'dmextract infile="source_{flux}_pileup_evt2.fits[energy=2000:8000,sky=region(sgra.reg)][bin time=::300]" outfile="marx_{flux}_sgra_2-8keV_lc300_pileup.fits" bkg="source_{flux}_pileup_evt2.fits[ccd_id=7,sky=region(bkg.reg)]" opt="ltc1" clobber=yes', shell=True, cwd=wd)

        plot_lightcurve("15043", f"marx_{flux}_sgra_2-8keV_lc300")
        plot_lightcurve("15043", f"marx_{flux}_sgra_2-8keV_lc300_pileup")

        difference(f"marx_{flux}_sgra_2-8keV_lc300", f"marx_{flux}_sgra_2-8keV_lc300_pileup", flux)

    differences.append(np.mean(np.array(differences_oneflux)))
    #print(differences_oneflux, np.mean(np.array(differences_oneflux)))



fluxes_ct = np.array(fluxes)*308.111


observed_analytical_fluxes = np.linspace(0.002, 1.15, 300)
frame_fluxes = observed_analytical_fluxes*0.44104
fd = 1-(((np.exp(frame_fluxes)-1)*np.exp(-frame_fluxes))/(frame_fluxes))
true_analytical_fluxes = frame_fluxes/(1-(fd))*1/0.44014
analytical_difference = true_analytical_fluxes - observed_analytical_fluxes

plt.clf()
plt.plot(fluxes_ct, np.array(differences), label='MARX Simulations')
plt.plot(true_analytical_fluxes, analytical_difference, label='Analytical Correction')
plt.title('Pileup Difference in ct/s')
plt.xlabel('True Flux (ct/s)')
plt.ylabel('Pileup Difference (ct/s)')
plt.legend()
plt.grid()
plt.savefig('PILEUP_DIFFERENCES.png')