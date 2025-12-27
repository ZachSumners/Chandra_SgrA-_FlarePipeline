import subprocess

obsIDs = [15043,
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
9173,
3663]

print(obsIDs)

for ids in obsIDs:
    if ids == 9173:
        continue
    id_5digit = str(ids).zfill(5)
    subprocess.call(f'dmkeypar /Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{ids}/repro/acisf{id_5digit}_bary_evt2.fits LIVETIME echo+', shell=True)