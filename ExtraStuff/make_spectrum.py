import subprocess
import pandas as pd 

df = pd.read_csv('flare_properties_nov19.csv')

print(df['fluence_ct'].dtype)
df = df[df['obs_id'] == 28231]

ids = df['obs_id'].values
start_flare = df['start_mjd'].values
end_flare = df['end_mjd'].values 

grating_ids = ['14392',
'13851',
'13843',
'13839',
'13845',
'14432',
'13842',
'13838',
'13849',
'13842']

flare_counter = 1
old_flare_id = 0
for i, flare_id in enumerate(ids):
     
    start = (start_flare[i] - 50814)*86400
    end = (end_flare[i] - 50814)*86400

    id_5digit = str(flare_id).zfill(5)

    if flare_id == old_flare_id:
        flare_counter += 1
    else:
        flare_counter = 1
    
    print(flare_id, flare_counter, start, end)
    #subprocess.call(f'rm flare_spectrum*', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    #subprocess.call(f'rm *_flare*_evt2.fits', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    
    subprocess.call(f'dmcopy infile="acisf{id_5digit}_bary_evt2.fits[EVENTS][time={start}:{end}]" outfile="{id_5digit}_flare{i}_evt2.fits" clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    if str(flare_id) in grating_ids:
        subprocess.call(f'specextract infile="{id_5digit}_flare{i}_evt2.fits[sky=region(order0.reg)]" outroot="flare_spectrum{flare_counter}" weight=no correctpsf=yes clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    else:
        subprocess.call(f'specextract infile="{id_5digit}_flare{i}_evt2.fits[sky=region(sgra.reg)]" bkgfile="{id_5digit}_flare{i}_evt2.fits[sky=region(bkg.reg)]" outroot="flare_spectrum{flare_counter}" grouptype=NONE binspec=NONE weight=no correctpsf=yes clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    
    old_flare_id = flare_id
    




import subprocess
import pandas as pd 

'''ids = [15043,
20751,
16218,
14392,
1561,
22595,
13851,
20346,
11843,
23739,
22594,
13843,
13839,
3393,
13845,
20041,
18055,
14432,
15042,
13842,
18055,
13838,
16966,
3393,
10556,
6363,
26760,
21454,
13849,
22592,
13842,
22595,
22937,
4684,
5953,
3663]

ids = [3393, 3393]

grating_ids = ['14392',
'13851',
'13843',
'13839',
'13845',
'14432',
'13842',
'13838',
'13849',
'13842']

flare_counter = 1
old_flare_id = 0

for i, flare_id in enumerate(ids):
    flare_id = str(flare_id)

    id_5digit = str(flare_id).zfill(5)

    if flare_id == old_flare_id:
        flare_counter += 1
    else:
        flare_counter = 1
    
    print(flare_id, flare_counter)

    subprocess.call(f'dmcopy infile="acisf{id_5digit}_bary_evt2.fits[EVENTS]" outfile="{id_5digit}_flare{i}_evt2.fits" clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    if str(flare_id) in grating_ids:
        subprocess.call(f'specextract infile="{id_5digit}_flare{i}_evt2.fits[sky=region(order0.reg)]" outroot="flare_spectrum{flare_counter}" weight=no correctpsf=yes clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    else:
        subprocess.call(f'specextract infile="{id_5digit}_flare{i}_evt2.fits[sky=region(sgra.reg)]" bkgfile="{id_5digit}_flare{i}_evt2.fits[sky=region(bkg.reg)]" outroot="flare_spectrum{flare_counter}" weight=no correctpsf=yes clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    
    old_flare_id = flare_id

'''   

'''
mag_ids = [15043, 16218, 15042, 16966]

for i, flare_id in enumerate(mag_ids):
    flare_id = str(flare_id)

    id_5digit = str(flare_id).zfill(5)

    if flare_id == old_flare_id:
        flare_counter += 1
    else:
        flare_counter = 1
    
    print(flare_id, flare_counter)

    subprocess.call(f'dmcopy infile="acisf{id_5digit}_bary_evt2.fits[EVENTS]" outfile="{id_5digit}_flare{i}_evt2.fits" clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    if str(flare_id) in grating_ids:
        subprocess.call(f'specextract infile="{id_5digit}_flare{i}_evt2.fits[sky=region(order0.reg)]" outroot="flare_spectrum{flare_counter}" weight=no correctpsf=yes clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    else:
        subprocess.call(f'specextract infile="{id_5digit}_flare{i}_evt2.fits[sky=region(mag.reg)]" bkgfile="{id_5digit}_flare{i}_evt2.fits[sky=region(bkg.reg)]" outroot="mag_spectrum{flare_counter}" weight=no correctpsf=yes clobber=yes', shell=True, cwd=f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_id}/repro')
    
    old_flare_id = flare_id
'''