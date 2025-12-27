import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from pycrates import read_file
import matplotlib.ticker as mtick

# ----------------------------
# Load one light curve (MJD + data)
# ----------------------------
def load_lc(file):
    tab = read_file(file)
    t_sec  = tab.get_column("time").values         # seconds since Chandra epoch
    mjd    = 50814.0 + (t_sec / 86400.0)           # -> MJD (days)
    rate   = tab.get_column("NET_RATE").values
    erate  = np.abs(tab.get_column("ERR_RATE").values)
    return mjd, rate, erate

# ----------------------------
# Compress gaps across *all* series
# thresh_min: gaps larger than this are compressed
# keep_min: small gap left in place (visual separator)
# ----------------------------
def compress_mjd_lists(mjd_lists, thresh_min=20.0, keep_min=2.0):
    thresh_d = thresh_min / (24.0 * 60.0)  # minutes -> days
    keep_d   = keep_min   / (24.0 * 60.0)

    # Concatenate while remembering segment membership
    seg_id = []
    all_mjd = []
    for i, arr in enumerate(mjd_lists):
        all_mjd.append(arr)
        seg_id.append(np.full(arr.size, i, dtype=int))
    all_mjd = np.concatenate(all_mjd)
    seg_id  = np.concatenate(seg_id)

    # Sort by time (preserve original order within equal times)
    order = np.argsort(all_mjd, kind="mergesort")
    mjd_sorted = all_mjd[order]

    # Build compression offsets
    offset = 0.0
    offsets = np.zeros_like(mjd_sorted)
    for k in range(1, mjd_sorted.size):
        dt = mjd_sorted[k] - mjd_sorted[k-1]
        if dt > thresh_d:                 # big gap detected
            offset += (dt - keep_d)       # remove everything except keep_d
        offsets[k] = offset
    mjd_compressed_sorted = mjd_sorted - offsets

    # Unsort back to original order
    unsort = np.empty_like(order)
    unsort[order] = np.arange(order.size)
    mjd_compressed = mjd_compressed_sorted[unsort]

    # Split back per segment
    out = []
    for i in range(len(mjd_lists)):
        out.append(mjd_compressed[seg_id == i])
    return out

# ----------------------------
# Pretty HH:MM formatter for float hours
# ----------------------------
def hhmm_formatter(x, pos):
    h = int(np.floor(x))
    m = int(round((x - h) * 60))
    if m == 60:
        h += 1
        m = 0
    return f"{h:02d}:{m:02d}"

# ----------------------------
# Main plotting
# ----------------------------
files = [
    "/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/31019/repro/bright_object_lc.fits",
    "/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/31020/repro/bright_object_lc.fits",
    "/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/31021/repro/bright_object_lc.fits",
    "/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/31022/repro/bright_object_lc.fits",
]

# Load all curves
mjds, rates, erates = [], [], []
for f in files:
    m, r, e = load_lc(f)
    mjds.append(m)
    rates.append(r)
    erates.append(e)

# Compress gaps (e.g., compress gaps > 20 min down to 2 min)
mjds_c = compress_mjd_lists(mjds, thresh_min=20.0, keep_min=50.0)

# Convert to elapsed hours from the very first (compressed) time
t0 = min(m[0] for m in mjds_c)
hours_lists = [ (m - t0) * 24.0 for m in mjds_c ]

# Plot
fig, ax = plt.subplots(figsize=(14, 8))
for h, r, e in zip(hours_lists, rates, erates):
    ax.errorbar(h, r, yerr=e, marker="o", linestyle="none",
                color="red", mfc="black", mec="black", ecolor="grey", alpha=0.9)

ax.set_xlabel("Time (compressed) — HH:MM")
ax.set_ylabel("Net Count Rate (counts/s)")
ax.xaxis.set_major_formatter(mtick.FuncFormatter(hhmm_formatter))
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.show()