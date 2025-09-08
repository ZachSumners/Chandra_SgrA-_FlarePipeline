import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pycrates import read_file
from astropy.time import Time
import matplotlib.dates as mdates
import glob
from scipy.interpolate import CubicSpline
from astropy.table import Table	
from astropy.io.fits import BinTableHDU

def difference(name1, name2, flux, differences_oneflux, trueflux_oneflux, observedflux_oneflux, wd):
    marx_folder = os.getcwd()
    
    #Open MARX fits files.
    tab1 = read_file(f"{marx_folder}/MARX_SupportFiles/{name1}.fits")
    tab2 = read_file(f"{marx_folder}/MARX_SupportFiles/{name2}.fits")

    #Convert the chandra MJD time to a readable time.
    chandra_time = tab1.get_column("time").values
    myFmt = mdates.DateFormatter('%H:%M:%S')
    chandra_mjd = (50814 + (0 + chandra_time)/86400)
    chandra_times = Time(chandra_mjd, format='mjd')
    chandra_datetimes = chandra_times.datetime

    #Get the flux rates and errors.
    try:
        rate1 = tab1.get_column("NET_RATE").values
        erate1 = tab1.get_column("ERR_RATE").values
        rate2 = tab2.get_column("NET_RATE").values
        erate2 = tab2.get_column("ERR_RATE").values
    except:
        rate1 = tab1.get_column("COUNT_RATE").values
        erate1 = tab1.get_column("COUNT_RATE_ERR").values
        rate2 = tab2.get_column("COUNT_RATE").values
        erate2 = tab2.get_column("COUNT_RATE_ERR").values

    #Calculate difference between piled and unpiled
    difference = rate1 - rate2
    difference_mean = np.mean(difference)
    differences_oneflux.append(difference_mean)

    #Save results, remove intermediate files.
    trueflux_oneflux.append(np.mean(rate1))
    observedflux_oneflux.append(np.mean(rate2))

    subprocess.call(f'rm {marx_folder}/MARX_SupportFiles/{name1}.fits', shell=True, cwd=wd)
    subprocess.call(f'rm {marx_folder}/MARX_SupportFiles/{name2}.fits', shell=True, cwd=wd)
    
    return differences_oneflux, trueflux_oneflux, observedflux_oneflux
    
def param_set(flux, tstart, exposuretime, gratingtype, detectortype, ditherfile, acis_exposure_time, sourcera, sourcedec, ra_nom, dec_nom, roll_nom, detoffsetx, detoffsety, detoffsetz, wd):
    #This is all just setting the params of MARX to match the ObservationID 
    marx_folder = os.getcwd()
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par SourceFlux={flux}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par TStart={tstart}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par ExposureTime={exposuretime}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par GratingType={gratingtype}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par DetectorType={detectortype}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par DitherFile={ditherfile}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par ACIS_Exposure_Time={acis_exposure_time}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par SourceRA={sourcera}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par SourceDEC={sourcedec}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par RA_Nom={ra_nom}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par Dec_Nom={dec_nom}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par Roll_Nom={roll_nom}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par DetOffsetX={detoffsetx}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par DetOffsetY={detoffsety}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par DetOffsetZ={detoffsetz}', shell=True, cwd=wd)
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marx.par OutputDir={marx_folder}/MARX_SupportFiles/output', shell=True, cwd=wd)
    
    #Set the pileup frame time
    subprocess.call(f'pset {marx_folder}/MARX_SupportFiles/marxpileup.par FrameTime={acis_exposure_time}', shell=True, cwd=wd)

def run_marx(flux, wd):
    marx_folder = os.getcwd()
    #Run the piled and unpiled versions of marx. Save result in fits file.
    subprocess.call(f'marx @@/{marx_folder}/MARX_SupportFiles/marx.par Verbose=no', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.call(f'marx2fits {marx_folder}/MARX_SupportFiles/output {marx_folder}/MARX_SupportFiles/source_{flux}_evt2.fits', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.call(f'marxpileup @@/{marx_folder}/MARX_SupportFiles/marxpileup.par MarxOutputDir={marx_folder}/MARX_SupportFiles/output Verbose=0', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.call(f'marx2fits --pileup {marx_folder}/MARX_SupportFiles/output/pileup {marx_folder}/MARX_SupportFiles/source_{flux}_pileup_evt2.fits', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def bouffard_version(acis_exposure_time):
    #Calculate and plot the Bouffard (2019)/Nowak (2012) pileup perscription
    alpha = 0.5
    true_analytical_fluxes = np.linspace(0.005, 2, 100)
    fraction = 1 - ((np.exp(alpha*true_analytical_fluxes*acis_exposure_time) - 1)*np.exp(-true_analytical_fluxes*acis_exposure_time))/(alpha*true_analytical_fluxes*acis_exposure_time)
    observed_analytical_fluxes = true_analytical_fluxes * (1-fraction)
    plt.plot(observed_analytical_fluxes, true_analytical_fluxes, label='Bouffard (2019) Correction (from Nowak 2012)')

def ponti_version(a, b, c, d, observed_analytical_fluxes):
    #Calculate and plot the Ponti (2015) pileup method
    ponti_true_flux = a* observed_analytical_fluxes**b + c* observed_analytical_fluxes**d
    plt.plot(observed_analytical_fluxes, ponti_true_flux, label='Ponti (2015) Correction')
    return ponti_true_flux

def neilsen_version(acis_exposure_time):
    #Calculate the Neilsen grating pileup method
    true_analytical_fluxes = np.linspace(0.005, 1, 300)
    alpha = 0.5
    observed_analytical_fluxes = ((np.array(true_analytical_fluxes)*acis_exposure_time) * (1 + 0.94/(np.array(true_analytical_fluxes)*alpha*acis_exposure_time) * (np.exp(alpha*np.array(true_analytical_fluxes)*acis_exposure_time) - 1)*np.exp(-np.array(true_analytical_fluxes)*acis_exposure_time))) * 0.52 * (1/acis_exposure_time)#(np.array(observedfluxes_order1)*acis_exposure_time)/(np.array(true_fluxes)*acis_exposure_time)) * (1/acis_exposure_time)
    plt.plot(observed_analytical_fluxes, true_analytical_fluxes, label='Neilsen (2013) Correction')

def marx_pileup_estimation(observationID, repro_wd):
    differences = []
    true_fluxes = []
    observed_fluxes = []
    fluxes = []
    wd = repro_wd

    #Open the event file of the ObsID to get observation information for MARX
    observationID_5digit = str(observationID).zfill(5)
    evt_hdul = fits.open(f'{repro_wd}/acisf{observationID_5digit}_bary_evt2.fits')
    header = evt_hdul[1].header

    marx_folder = os.getcwd()

    #Get observations params
    tstart = Time((50814 + (0 + header['TSTART'])/86400), format='mjd').decimalyear
    exposuretime = float(header['TSTOP']) - float(header['TSTART'])
    gratingtype = header['GRATING']
    detectortype = header['DETNAM']
    
    if ("0" in detectortype) or ("1" in detectortype) or ("2" in detectortype) or ("3" in detectortype):
        detectortype = 'ACIS-I'
    elif ("4" in detectortype) or ("5" in detectortype) or ("6" in detectortype) or ("7" in detectortype) or ("8" in detectortype) or ("9" in detectortype):
        detectortype = 'ACIS-S'

    folder = repro_wd
    pattern = os.path.join(folder, "*corrected.fits")
    files = glob.glob(pattern)
    ditherfile = files[0]

    acis_exposure_time = header['EXPTIME']
    sourcera = header['RA_TARG']
    sourcedec = header['DEC_TARG']
    ra_nom = header['RA_NOM']
    dec_nom = header['DEC_NOM']
    roll_nom = header['ROLL_NOM']

    if detectortype == 'ACIS-I':
        detoffsetx = -0.78234819833843 - float(header['SIM_X'])
        detoffsety = 0 - float(header['SIM_Y'])
        detoffsetz = -233.5924630914 - float(header['SIM_Z'])
    elif detectortype == 'ACIS-S':
        detoffsetx = -0.68426746699586 - float(header['SIM_X'])
        detoffsety = 0 - float(header['SIM_Y'])
        detoffsetz = -190.1325231040 - float(header['SIM_Z'])
    else:
        raise ValueError('Detector type not valid.')

    #This is the main loop that runs the simulation. For varying flux, find how much pileup there is.
    for i in range(1, 602, 20):
        #Set the flux (0.00005 is ~Sgr A* quiescense)
        flux = 0.00002*i
        fluxes.append(flux)

        differences_oneflux = []
        trueflux_oneflux = []
        observedflux_oneflux = []
        #Run the simulation multiple times for a given flux to average out the poisson noise that the simulation generates.
        for i in range(2):
            #Set params for MARX
            param_set(flux, tstart, exposuretime, gratingtype, detectortype, ditherfile, acis_exposure_time, sourcera, sourcedec, ra_nom, dec_nom, roll_nom, detoffsetx, detoffsety, detoffsetz, wd)
            #Run marx simulation
            run_marx(flux, wd)

            #Extract count rates from the piled and unpiled versions of the simulation. Grating has no background so needs own method to extract.
            subprocess.call('punlearn dmextract', shell=True, cwd=wd)
            if gratingtype == 'NONE':
                subprocess.call(f'dmextract infile="{marx_folder}/MARX_SupportFiles/source_{flux}_evt2.fits[energy=2000:8000,sky=region({repro_wd}/sgra.reg)][bin time=::300]" outfile="{marx_folder}/MARX_SupportFiles/marx_{flux}_sgra_2-8keV_lc300.fits" bkg="{marx_folder}/MARX_SupportFiles/source_{flux}_evt2.fits[sky=region({repro_wd}/bkg.reg)]" opt="ltc1" clobber=yes', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                subprocess.call(f'dmextract infile="{marx_folder}/MARX_SupportFiles/source_{flux}_pileup_evt2.fits[energy=2000:8000,sky=region({repro_wd}/sgra.reg)][bin time=::300]" outfile="{marx_folder}/MARX_SupportFiles/marx_{flux}_sgra_2-8keV_lc300_pileup.fits" bkg="{marx_folder}/MARX_SupportFiles/source_{flux}_pileup_evt2.fits[sky=region({repro_wd}/bkg.reg)]" opt="ltc1" clobber=yes', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            elif gratingtype == 'HETG':
                subprocess.call(f'dmextract infile="{marx_folder}/MARX_SupportFiles/source_{flux}_evt2.fits[energy=2000:8000,sky=region({repro_wd}/order1and0.reg)][bin time=::300]" outfile="{marx_folder}/MARX_SupportFiles/marx_{flux}_sgra_2-8keV_lc300.fits" opt="ltc1" clobber=yes', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                subprocess.call(f'dmextract infile="{marx_folder}/MARX_SupportFiles/source_{flux}_pileup_evt2.fits[energy=2000:8000,sky=region({repro_wd}/order1and0.reg)][bin time=::300]" outfile="{marx_folder}/MARX_SupportFiles/marx_{flux}_sgra_2-8keV_lc300_pileup.fits" opt="ltc1" clobber=yes', shell=True, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            #Calculate the difference between the true and observed count rate.
            differences_oneflux, trueflux_oneflux, observedflux_oneflux = difference(f"marx_{flux}_sgra_2-8keV_lc300", f"marx_{flux}_sgra_2-8keV_lc300_pileup", flux, differences_oneflux, trueflux_oneflux, observedflux_oneflux, wd)
            #Remove intermediate files
            subprocess.call(f'rm {marx_folder}/MARX_SupportFiles/source_{flux}_evt2.fits', shell=True, cwd=wd)
            subprocess.call(f'rm {marx_folder}/MARX_SupportFiles/source_{flux}_pileup_evt2.fits', shell=True, cwd=wd)

        #Keeps track of count rates.
        differences.append(np.mean(np.array(differences_oneflux)))
        true_fluxes.append(np.mean(np.array(trueflux_oneflux)))
        observed_fluxes.append(np.mean(np.array(observedflux_oneflux)))

    observed_analytical_fluxes = np.linspace(0.005, 1.25, 300)
    
    #Plot MARX results
    plt.figure(figsize=(14, 10))
    plt.plot(np.array(observed_fluxes), np.array(true_fluxes), label='MARX Simulations')

    #Plot analytical solution for comparison
    if gratingtype == 'NONE':
        bouffard_version(acis_exposure_time)
    elif gratingtype == 'HETG':
        neilsen_version(acis_exposure_time)

    #Set params for and plot Ponti (2015) version of pileup
    if detectortype == 'ACIS-I' and gratingtype == 'NONE':
        a = 1.563
        b = 1.099
        c = 1185
        d = 4.866
    elif detectortype == 'ACIS-S' and gratingtype == 'NONE':
        a = 0.2366
        b = 6.936
        c = 1.393
        d = 1.179
    elif detectortype == 'ACIS-S' and gratingtype == 'HETG':
        a = 802
        b = 4.743
        c = 1.599
        d = 1.11
    ponti_true_flux = ponti_version(a, b, c, d, observed_analytical_fluxes)

    #Plot WebPIMMS version of pileup for a given few observations (the brightest ObsID of the observing mode)
    if observationID == 23739:
        plt.plot([4.9E-3, 9.7E-3, 4.78E-2, 9.35E-2, 1.787E-1, 2.561E-1, 3.263E-1, 3.9E-1, 4.476E-1, 4.996E-1, 5.464E-1, 5.885E-1, 6.262E-1, 6.599E-1, 6.901E-1, 7.168E-1, 7.405E-1, 7.615E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5], label='WebPIMMS')
    elif observationID == 26760:
        plt.plot([0.005, 0.01, 0.049, 9.62E-2, 1.849E-1, 2.668E-1, 3.423E-1, 4.118E-1, 4.756E-1, 5.342E-1, 5.88E-1, 6.372E-1, 6.821E-1, 7.231E-1, 7.604E-1, 7.944E-1, 8.251E-1, 8.53E-1, 8.781E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6], label='WebPIMMS')
    elif observationID == 242:
        plt.plot([4.8E-3, 9.6E-3, 4.33E-2, 7.68E-2, 1.217E-1, 1.462E-1, 1.582E-1, 1.628E-1, 1.636E-1, 1.627E-1, 1.615E-1, 1.607E-1, 1.608E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], label='WebPIMMS')
    elif observationID == 1561:
        plt.plot([4.8E-3, 9.6E-3, 4.33E-2, 7.68E-2, 1.217E-1, 1.462E-1, 1.582E-1, 1.628E-1, 1.636E-1, 1.627E-1, 1.615E-1, 1.607E-1, 1.608E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], label='WebPIMMS')
    elif observationID == 10556:
        plt.plot([4.8E-3, 9.5E-3, 4.33E-2, 7.72E-2, 1.232E-1, 1.49E-1, 1.621E-1, 1.676E-1, 1.689E-1, 1.683E-1, 1.67E-1, 1.66E-1, 1.658E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], label='WebPIMMS')
    elif observationID == 5950:
        plt.plot([4.8E-3, 9.5E-3, 4.34E-2, 7.73E-2, 1.233E-1, 1.491E-1, 1.622E-1, 1.676E-1, 1.689E-1, 1.682E-1, 1.67E-1, 1.66E-1, 1.658E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], label='WebPIMMS')
    elif observationID == 9169:
        plt.plot([4.8E-3, 9.5E-3, 4.32E-2, 7.66E-2, 1.215E-1, 1.46E-1, 1.58E-1, 1.628E-1, 1.636E-1, 1.627E-1, 1.615E-1, 1.607E-1, 1.608E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], label='WebPIMMS')
    elif observationID == 15043:
        plt.plot([4.9E-3, 9.8E-3, 4.82E-2, 9.5E-2, 1.842E-1, 2.68E-1, 3.466E-1, 4.204E-1, 4.895E-1, 5.542E-1, 6.148E-1, 6.714E-1, 7.242E-1, 7.735E-1, 8.195E-1, 8.624E-1, 9.023E-1, 9.394E-1, 9.738E-1, 1.0057E0, 1.0352E0, 1.0625E0, 1.0877E0, 1.1109E0, 1.1323E0, 1.1519E0, 1.1697E0, 1.1861E0], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5], label='WebPIMMS')
    elif observationID == 14941:
        plt.plot([4.8E-3, 9.5E-3, 4.33E-2, 7.71E-2, 1.231E-1, 1.489E-1, 1.621E-1, 1.676E-1, 1.689E-1, 1.683E-1, 1.67E-1, 1.661E-1, 1.658E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], label='WebPIMMS')
    elif observationID == 13839:
        plt.plot([4.8E-3, 9.4E-3, 4.3E-2, 7.66E-2, 1.225E-1, 1.485E-1, 1.618E-1, 1.675E-1, 1.689E-1, 1.683E-1, 1.671E-1, 1.661E-1, 1.658E-1], [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], label='WebPIMMS')

    #Plot where the brightest count rate of the ObsID is to give sense of maxima.
    #obs_lc_hdul = fits.open(f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{observationID}/repro/{observationID_5digit}_sgra_2-8keV_lc300_pileup.fits")
    #lc = obs_lc_hdul[1].data

    #try:
    #    plt.axvline(x=np.max(lc['NET_RATE']), color='black', linestyle='--', linewidth=2, label='Max Count Rate of ObsID')
    #except:
    #    plt.axvline(x=np.max(lc['COUNT_RATE']), color='black', linestyle='--', linewidth=2, label='Max Count Rate of ObsID')

    #Plotting all the solutions together
    plt.title(f'Pileup for ObsID {observationID} in ct/s')
    plt.xlabel('Observed Flux (ct/s)')
    plt.ylabel('True Flux (ct/s)')
    plt.legend()
    plt.grid()

    #Setting plot bounds
    if gratingtype == 'NONE':
        if np.max(np.array(true_fluxes)) > 5:
            plt.ylim(-0.05, 5)
        if  np.max(ponti_true_flux) > 5:
            plt.ylim(-0.05, np.max(np.array(true_fluxes)) + 0.1)
        if np.max(np.array(observed_fluxes)) > 0.5:
            plt.xlim(-0.05, 1.25)
        if np.max(np.array(observed_fluxes)) < 0.4:
            plt.xlim(-0.01, 0.2)
    elif gratingtype == 'HETG':
        plt.ylim(-0.02, 0.8)
        plt.xlim(-0.02, 0.4)

    pileup_conversion = np.column_stack((np.array(observed_fluxes), np.array(true_fluxes)))
    np.savetxt(f"{repro_wd}/marx_pileup_conversion.txt", pileup_conversion, fmt="%.8f", delimiter="\t")

    plt.savefig(f'PILEUP_DIFFERENCES_{observationID}.png')

    return np.array(observed_fluxes), np.array(true_fluxes)

def marx_pileup_interpolation(marx_observed_flux, marx_true_flux, observationID, erange, tbin, fileName, region_name, repro_wd):
    #Open the lightcurve.

    f = fits.open(f'{repro_wd}/{observationID}_{region_name}_{erange[0]}-{erange[1]}keV_lc{tbin}.fits')
    table = Table(f[1].data)

	#Open the event files to get the exposure time of the observation.
    f_evt = fits.open(f'{repro_wd}/acisf{observationID}_{fileName}_evt2.fits')
    exptime = float(f_evt[1].header['EXPTIME'])

    count_rate = table['COUNT_RATE']
    count_error = table['COUNT_RATE_ERR']
    
    #The curve is possibly non monotomic (turn over at burnout) so we need to parameterize the interpolation splines
    true_flux = []
    t = np.arange(len(marx_observed_flux))
    spl_x = CubicSpline(t, marx_observed_flux)
    spl_y = CubicSpline(t, marx_true_flux)

    ts = np.linspace(t[0], t[-1], 10000)
    
    for flux in count_rate:
        find_t = spl_x(ts) - flux
        t0 = np.where(find_t == find_t[find_t > 0].min())[0][0]
        
        x0 = spl_x(ts[t0])
        y0 = spl_y(ts[t0])

        true_flux.append(y0)


    table.add_column(np.asarray(true_flux), index=21, name='RATE_PILEUP')
    table.add_column(count_error, index=22, name='PILEUP_ERR')

	#Save this data into a new fits file.
	# Create new HDU with updated table
    hdu1 = BinTableHDU(data=table, header=f[1].header)
    hdu1.header['MJDREF'] = 5.08140000000000E+04  # Add/update MJDREF

	# Assemble final file with original PRIMARY and GTI extensions
    hdul = fits.HDUList([f[0], hdu1, f[2]])
    hdul.writeto(f'{repro_wd}/{observationID}_{region_name}_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits', overwrite=True)


def marx_pileup_interpolation_block(marx_observed_flux, marx_true_flux, count_rate):  
    #The curve is possibly non monotomic (turn over at burnout) so we need to parameterize the interpolation splines
    true_flux = []
    t = np.arange(len(marx_observed_flux))
    spl_x = CubicSpline(t, marx_observed_flux)
    spl_y = CubicSpline(t, marx_true_flux)

    ts = np.linspace(t[0], t[-1], 10000)
    for flux in count_rate:
        #print(flux)
        find_t = spl_x(ts) - flux
        try:
            t0 = np.where(find_t == find_t[find_t > 0].min())[0][0]
            x0 = spl_x(ts[t0])
            y0 = spl_y(ts[t0])
        except:
            y0 = flux

        true_flux.append(y0)

    return np.asarray(true_flux)