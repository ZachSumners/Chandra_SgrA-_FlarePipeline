'''sources = fits.open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/src.fits')
sources_table = Table(sources[1].data)
sources_table.sort('NET_COUNTS')
sources_table.reverse()

ra_catalog = []
dec_catalog = []
ra_reported = []
dec_reported = []
ra_reported_error = []
dec_reported_error = []

for i in range(3):
	p = subprocess.run(f'search_csc pos="{sources_table["RA"][i+1]}, {sources_table["DEC"][i+1]}" radius="0.02" outfile="{observationID}_wcsLookup{i}"', shell=True, cwd=repro_wd, capture_output=True)
	result = p.stdout.decode("utf-8")
	
	f = open(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/{observationID}_wcsLookup{i}', 'r')
	lines = f.readlines()

	ra_catalog.append(lines[65].split()[3])
	dec_catalog.append(lines[65].split()[4])
	ra_reported.append(sources_table["RA"][i+1])
	ra_reported_error.append(sources_table["RA_ERR"][i+1])
	dec_reported.append(sources_table["DEC"][i+1])
	dec_reported_error.append(sources_table["DEC_ERR"][i+1])
	
catalog_wcs = {'ra': ra_catalog, 'dec': dec_catalog}
df = pd.DataFrame(data=catalog_wcs)
df.to_csv(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/{observationID}_src_catalog_wcs.csv', index=False, sep='\t')

reported_wcs = {'ra': ra_reported, 'RA_ERR': ra_reported_error, 'dec': dec_reported, 'DEC_ERR': dec_reported_error}
df = pd.DataFrame(data=reported_wcs)
df.to_csv(f'/home/zach/Desktop/McGill-MSc/Chandra/data/{observationID}/repro/{observationID}_src_reported_wcs.csv', index=False, sep='\t')

print(reported_wcs, catalog_wcs)

col1 = fits.Column(name='ra', array=ra_catalog, format='E')
col2 = fits.Column(name='dec', array=dec_catalog, format='E')
hdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1, col2]))

hdu.writeto(f'{repro_wd}/{observationID}_src_catalog_wcs.fits', overwrite=True)
'''
