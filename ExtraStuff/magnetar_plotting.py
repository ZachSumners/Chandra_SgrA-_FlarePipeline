#!/usr/bin/env python3
import os
import re
import csv
from astropy.time import Time   # pip install astropy

import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd

# change this to the directory that contains all your observation folders
ROOT = "/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub"

# pattern for files like 16218_MAGNETAR_TABLE_RESULTS.txt
RESULTS_FILENAME_SUFFIX = "_MAGNETAR_TABLE_RESULTS.txt"

out_rows = []

for dirpath, dirnames, filenames in os.walk(ROOT):
    # only care about directories that are called Results
    if os.path.basename(dirpath) != "Results":
        continue

    for fname in filenames:
        if not fname.endswith(RESULTS_FILENAME_SUFFIX):
            continue

        fullpath = os.path.join(dirpath, fname)
        with open(fullpath, "r") as f:
            text = f.read()

        # Obs_ID (optional but handy)
        m_id = re.search(r"Obs_ID:\s*(\S+)", text)
        obs_id = m_id.group(1) if m_id else "UNKNOWN"

        # Obs Date (ISO string)
        m_date = re.search(r"Obs Date:\s*(\S+)", text)
        if not m_date:
            print(f"Warning: no Obs Date in {fullpath}")
            continue
        iso_date = m_date.group(1)

        # convert to MJD
        t = Time(iso_date, format="isot", scale="utc")
        mjd = t.mjd

        # Quiescent Count Rate line:
        # e.g. "Quiescent Count Rate (10^-3 ct/s): 64.044 +/- 1.264"
        m_rate = re.search(
            r"Quiescent Count Rate.*:\s*([0-9.eE+-]+)\s*\+/-\s*([0-9.eE+-]+)",
            text,
        )
        if not m_rate:
            print(f"Warning: no Quiescent Count Rate in {fullpath}")
            continue

        rate = float(m_rate.group(1))       # in 10^-3 ct/s, as in the file
        rate_err = float(m_rate.group(2))

        out_rows.append({
            "obs_id": obs_id,
            "mjd": mjd,
            "rate_1e-3_ct_s": rate,
            "rate_err_1e-3_ct_s": rate_err,
            "file": fullpath,
        })

# write CSV
out_csv = "magnetar_fluxes.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "obs_id", "mjd", "rate_1e-3_ct_s", "rate_err_1e-3_ct_s", "file"
    ])
    writer.writeheader()
    writer.writerows(out_rows)

print(f"Wrote {len(out_rows)} rows to {out_csv}")


df = pd.read_csv('magnetar_fluxes.csv')

df = df.sort_values(by=['mjd'])

date = df['mjd'] - 56407
rate = df['rate_1e-3_ct_s']
rate_err = df['rate_err_1e-3_ct_s']

markersize = 6
capsize = 2
capthick = 0.6
elinewidth= 0.6
edgewidth = 0.6

fontsize = 12

plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.titlesize'] = 10


plt.rcParams['xtick.major.size'] = 8  # length of major ticks
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.minor.size'] = 5

plt.rcParams.update({'text.usetex' : True, 'font.family' : 'Computer Modern Roman'})

fontsize = 14  # or whatever you use


date = date.to_numpy()
rate = rate.to_numpy()
rate_err = rate_err.to_numpy()

import numpy as np
import matplotlib.pyplot as plt

radius = 2.523467125e22

# count-rate -> flux conversion
def rate_to_flux(r):
    # r is in ct/s already here (your left axis units)
    L = r * (10**34 / 0.0118766)          # erg/s
    return L / (4 * np.pi * radius**2) /1e-12    # erg/cm^2/s

def flux_to_rate(f):
    L = f * (4 * np.pi * radius**2)
    return L / (10**34 / 0.0118766)

fig, ax = plt.subplots(figsize=(5, 4))

# left axis: count rate in ct/s
ax.errorbar(date, rate/1000.0, yerr=rate_err/1000.0,
            fmt='-o', mec='black', mfc='steelblue')
ax.set_xlabel('Days Since Outburst')
ax.set_ylabel('Count Rate (ct/s)')
ax.set_ylim(bottom=0)   # avoid negative counts/ticks

# right axis: flux
ax2 = ax.secondary_yaxis('right',
                         functions=(rate_to_flux, flux_to_rate))
ax2.set_ylabel(r'Flux (10$^{-12}$ erg cm$^{-2}$ s$^{-1}$)')

plt.savefig("MagnetarFlux.pdf", dpi=300, bbox_inches='tight')