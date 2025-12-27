import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Example data — replace these with your actual lightcurve arrays
lc1 = fits.open('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/28229/repro/28229_sgra_2-8keV_lc300_pileup.fits')[1].data
lc2 = fits.open('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/28230/repro/28230_sgra_2-8keV_lc300_pileup.fits')[1].data
lc3 = fits.open('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/28231/repro/28231_sgra_2-8keV_lc300_pileup.fits')[1].data
lc4 = fits.open('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/28232/repro/28232_sgra_2-8keV_lc300_pileup.fits')[1].data

lc1['RATE_PILEUP'][2] = 0
lc1['PILEUP_ERR'][2] = 0

# 2x2 subplot grid
fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=False, sharey=True)

# Plotting each lightcurve
axes[0, 0].errorbar((lc1['TIME']/86400 + 50814 - 60404)*24, lc1['RATE_PILEUP'], yerr=lc1['PILEUP_ERR'], label='April 4, 2024', fmt='o-', ecolor='black')
axes[0, 1].errorbar((lc2['TIME']/86400 + 50814 - 60406)*24, lc2['RATE_PILEUP'], yerr=lc2['PILEUP_ERR'], label='April 6, 2024', fmt='o-', ecolor='black')
axes[1, 0].errorbar((lc3['TIME']/86400 + 50814 - 60408)*24, lc3['RATE_PILEUP'], yerr=lc3['PILEUP_ERR'], label='April 8, 2024', fmt='o-', ecolor='black')
axes[1, 1].errorbar((lc4['TIME']/86400 + 50814 - 60409)*24, lc4['RATE_PILEUP'], yerr=lc4['PILEUP_ERR'], label='April 9, 2024', fmt='o-', ecolor='black')



# Optional: add legends
for ax in axes.flat:
    ax.legend()
    ax.set_xlabel('Hour (UTC)')
    ax.set_ylabel('Count Rate (ct/s)')

plt.tight_layout()
plt.show()