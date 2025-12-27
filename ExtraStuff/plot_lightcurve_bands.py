import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import pandas as pd
import os
from astropy.io import fits

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


def plot_lightcurve_with_ring_bands(t, y, yerr, peak_indices=None,
                                    bands=(5, 10, 15), obsid=None, shape=None):
    """
    Plot a lightcurve and add non-overlapping full-height coloured bands
    corresponding to ±N points around each peak.

    The smallest band (e.g. ±5) is the central coloured bar.
    The next (±10) colours only the time range between ±5 and ±10
    (two side strips).
    The largest (±15) colours only the time range between ±10 and ±15
    (outer side strips).

    Parameters
    ----------
    t : array-like
        Time array.
    y : array-like
        Lightcurve values.
    peak_indices : array-like of int, optional
        Indices of flare peaks. If None, peaks are found automatically.
    bands : iterable of int
        Numbers of points to include on each side of the peak, e.g. (5,10,15).
        Assumed to be increasing.
    """
    t = np.asarray(t)
    y = np.asarray(y)

    if peak_indices is None:
        peak_indices, _ = find_peaks(y)
    peak_indices = np.asarray(peak_indices, dtype=int)

    bands = sorted(bands)              # ensure ascending
    nb = len(bands)

    # colours for each band radius (inner → outer)
    colors = ["#b6d7a8", "#9fc5e8", "#f9cb9c"]  # green, blue, orange
    while len(colors) < nb:
        colors.append(colors[-1])      # just in case

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(t, y, yerr=yerr, marker='.', lw=1, label="Lightcurve")

    ymin, ymax = ax.get_ylim()

    for p in peak_indices:
        # mark the peak
        ax.plot(t[p], y[p], "ko", marker='.', ms=6, zorder=5)

        # pre-compute left/right times for each band radius
        lefts = []
        rights = []
        for n_side in bands:
            i_left = max(0, p - n_side)
            i_right = min(len(t) - 1, p + n_side)
            lefts.append(t[i_left])
            rights.append(t[i_right])

        # build ring-like, non-overlapping intervals
        prev_left = None
        prev_right = None
        for i, n_side in enumerate(bands):
            left = lefts[i]
            right = rights[i]

            col = colors[i]

            if prev_left is None:
                # innermost band: just [left, right]
                intervals = [(left, right)]
                label = f"±{n_side*300} seconds"
            else:
                # only shade the region between this band and the previous one
                intervals = []
                if left < prev_left:
                    intervals.append((left, prev_left))
                if prev_right < right:
                    intervals.append((prev_right, right))
                label = f"±{n_side*300} seconds"

            for (L, R) in intervals:
                if L >= R:
                    continue
                ax.axvspan(
                    L, R,
                    ymin, 5,
                    color=col,
                    alpha=0.8,
                    #label=label if p == peak_indices[0] else None,
                )

            prev_left = left
            prev_right = right

    ax.axhspan(0, 0.005*3, xmin=0, xmax=1, color='gray', alpha=0.4, label='Quiescence')
    
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(t[p] - 6000, t[p] + 6000)
    ax.set_xlabel("Time")
    ax.set_ylabel("Count rate")
    ax.set_title(f"Lightcurve ObsID: {obsid}, Shape: {shape}")
    ax.legend(loc="upper right")

    

    fig.tight_layout()
    return fig, ax


df = pd.read_csv('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/flare_properties_dec2.csv')
df = df.sort_values(by="obs_id")


flare_ids = df['obs_id'].to_numpy()
flare_numbers = df['flare_number'].to_numpy()
start = df['start_mjd'].to_numpy()
end = df['end_mjd'].to_numpy()
mean_rate = df['rate_mean'].to_numpy()
mean_rate_err = df['rate_mean_err'].to_numpy()
max_rate = df['rate_max'].to_numpy()
max_rate_err = df['rate_max_err'].to_numpy()
duration = df['duration_s'].to_numpy()
fluence = df['fluence_ct'].to_numpy()
#shape= df['Shape'].to_numpy()
#hr = df['hardness_ratio'].to_numpy()

wd = os.getcwd()

'''for i in range(len(flare_ids)):
    observationID_5digit = str(flare_ids[i]).zfill(5)
    hdul = fits.open(f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/{flare_ids[i]}/repro/{observationID_5digit}_sgra_2-8keV_lc300_pileup.fits")
    data = hdul[1].data 
    
    mask = ((50814 + data['time']/86400) > float(start[i])) & ((50814 + data['time']/86400) < float(end[i]))
    flare = data[mask]['RATE_PILEUP']

    max_arg = np.argmax(flare) + list(mask).index(True)

    plot_lightcurve_with_ring_bands(data['time'], data['RATE_PILEUP'], peak_indices=[max_arg], bands=(15, 10, 5), obsid=flare_ids[i], shape=shape[i])

    plt.show()
    '''

#highlighted = [13854, 16966, 13849, 20346, 6363, 4684, 11843, 13854, 13839, 20041, 13852, 22594, 18055, 13845, 14432, 22595, 15570, 13843, 3393, 16218, 25297, 23739, 15045, 5953, 15042, 1561, 14392, 15043, 20751, 13851]

highlighted = [
13849, 13851, 14463
]

from matplotlib.backends.backend_pdf import PdfPages

with PdfPages("new_lightcurves.pdf") as pdf:
    for i in range(len(flare_ids)):
        if flare_ids[i] not in highlighted:
            continue
        observationID_5digit = str(flare_ids[i]).zfill(5)
        hdul = fits.open(
            f"/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/"
            f"{flare_ids[i]}/repro/{observationID_5digit}_sgra_2-8keV_lc300_pileup.fits"
        )
        data = hdul[1].data 
        
        mask = ((50814 + data['time']/86400) > float(start[i])) & \
            ((50814 + data['time']/86400) < float(end[i]))
        flare = data[mask]['RATE_PILEUP']
        max_arg = np.argmax(flare) + list(mask).index(True)

        fig, ax = plot_lightcurve_with_ring_bands(
            data['time'], data['RATE_PILEUP'], data['PILEUP_ERR'],
            peak_indices=[max_arg],
            bands=(15, 10, 5),
            obsid=flare_ids[i],
            #shape=shape[i],
        )

        outname = f"lightcurve_{flare_ids[i]}_flare{flare_numbers[i]}.pdf"
        fig.savefig(outname, bbox_inches="tight")   # <-- saves as PDF
        plt.close(fig)   # optional, keeps memory down

        fig, ax = plot_lightcurve_with_ring_bands(
            data['time'], data['RATE_PILEUP'], data['PILEUP_ERR'],
            peak_indices=[max_arg],
            bands=(15, 10, 5),
            obsid=flare_ids[i],
            #shape=shape[i],
        )

        pdf.savefig(fig, bbox_inches="tight")   # add page to PDF
        plt.close(fig)