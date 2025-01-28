import numpy as np
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
import subprocess
import re

observationID = 28229
data_dir = f'/home/zach/Desktop/McGill-MSc/Chandra/data/Barycentric/{observationID}/{observationID}_sgra_2-8keV_lc300.fits'
pimms_dir = f'/home/zach/Desktop/Software/pimms/pimms'
script_directory = f'/home/zach/Desktop/McGill-MSc/Chandra/data/'
data_output = f'/home/zach/Desktop/McGill-MSc/Chandra/data/BB/{observationID}'

f = fits.open(data_dir)
table = Table(f[1].data)

def extract_fluxes(file_path):
	fluxes = []
	with open(file_path, 'r') as file:
		for line in file:
			# Match the line with flux values
			match = re.search(r'PIMMS predicts a flux.*?of\s+([\d\.E\-\+]+)\s+ergs/cm/cm/s', line)
			if match:
				flux_value = float(match.group(1))
				fluxes.append(flux_value)
	return np.array(fluxes)

def counts_to_flux(pimms_directory, script_directory, observationID, counts_np):
	with open(f'{script_directory}/pimms_commands_chandra.txt', 'w') as f:
		f.write("MODEL PL 2 14.2e+22\n") # See what to put here
		f.write("FROM CHANDRA ACIS-S\n")
		f.write("INSTRUMENT FLUX ERGS 2.0-10.0\n")
	    
		for count in counts_np:
			f.write(f"GO {count}\n")

	print("pimms_commands_chandra.txt file created successfully.")

	subprocess.call(f'{pimms_directory} < {script_directory}/pimms_commands_chandra.txt > {script_directory}/pimms_output_chandra.txt', shell = True)

	# Function to extract flux values from PIMMS output

	# Path to your PIMMS output file
	output_file_path = f'{script_directory}/pimms_output_chandra.txt'

	# Extract flux values and store in a numpy array
	flux = extract_fluxes(output_file_path)
	return flux

counts = table['NET_RATE']
counts_np = np.array(counts)
flux = counts_to_flux(pimms_dir, script_directory, observationID, counts_np)

counts_error = table['ERR_RATE']
counts_error_np = np.array(counts_error)
counts_error_np = np.nan_to_num(counts_error_np)
flux_error = counts_to_flux(pimms_dir, script_directory, observationID, counts_error_np)

luminosity = flux * 4 * np.pi * 2.523E22**2

r = 2.532E22
ur = 2.16E20
luminosity_error = abs(luminosity*np.sqrt((flux_error/flux)**2 + (ur/r)**2)) #According to the CompleteBB script.

mjdref = f[0].header['mjdref']
timezero = f[0].header['timezero']

time = table['TIME']
utc_time = mjdref + (timezero + time)/86400

print(np.mean(luminosity))

table.add_column(utc_time, index=4, name='TIME_MJD')
table.add_column(flux, index=22, name='FLUX')
table.add_column(flux_error, index=23, name='FLUX_ERR')
table.add_column(luminosity, index=24, name='LUMINOSITY')
table.add_column(luminosity_error, index=25, name='LUMINOSITY_ERR')

hdu = table.write(f'{data_output}/{observationID}_sgra_2-8keV_lc300_TABLE_lum.fits', overwrite=True, format='fits')
ft = fits.open(f'{data_output}/{observationID}_sgra_2-8keV_lc300_TABLE_lum.fits')

hdul = fits.HDUList([f[0], ft[1], f[2]])
hdul.writeto(f'{data_output}/{observationID}_sgra_2-8keV_lc300_mjd_lum.fits', overwrite=True)


