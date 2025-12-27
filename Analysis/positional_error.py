import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
from astropy.visualization import simple_norm
from crates_contrib.utils import *
import subprocess
from astropy.wcs import WCS
import os
import pandas as pd

obs_ids = [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(),d)) and d.isdigit()]
ra_coords = []
dec_coords = []

for ids in obs_ids:
    if ids == str(1561) or float(ids) < 14702:
        continue
    repro_wd = f'/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{ids}/repro'
    src_coords = [266.4086250, -29.0061944]

    observationID_5digit = str(ids).zfill(5)

    tr = SimpleCoordTransform(f'{repro_wd}/{observationID_5digit}_broad_thresh_img.fits')
    #converts coordinates in degrees to pixel coordinates
    ra_px, dec_px = tr.convert('world', 'physical', src_coords[0], src_coords[1])

    # Read the file
    with open(f"{repro_wd}/src.reg") as f:
        rows = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Remove "ellipse(" and ")"
            values = line.replace("ellipse(", "").replace(")", "")
            # Split by comma and convert to float
            nums = [float(v) for v in values.split(",")]
            rows.append(nums)

    # Turn into a DataFrame with column names
    df = pd.DataFrame(rows, columns=["x", "y", "a", "b", "angle"])

    distance = np.sqrt((df['x'].to_numpy() - ra_px)**2 + (df['y'].to_numpy() - dec_px)**2)
    src_idx = np.argmin(distance)
    if distance[src_idx] > 2:
        continue
    
    image_ra = df['x'][src_idx]
    image_dec = df['y'][src_idx]

    ra_diff = abs(image_ra - ra_px)
    dec_diff = abs(image_dec - dec_px)

    ra_coords.append(ra_diff)
    dec_coords.append(dec_diff)
print(len(ra_coords))
print(np.std(ra_coords), np.std(dec_coords))
plt.hist(ra_coords, bins=np.linspace(0, 10, 100))
plt.show()