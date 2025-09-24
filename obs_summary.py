import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import LogNorm, SymLogNorm
import matplotlib.patches as patches
import matplotlib.image as mpimg
import textwrap
import os
import re
import math
import CompleteBB as bb

def set_column_widths(tbl, col_widths):
    """
    col_widths: dict {col_index: width_fraction}
    Works with matplotlib.table.Table even with header row at -1.
    """
    for (r, c), cell in tbl.get_celld().items():
        if c in col_widths:
            cell.set_width(col_widths[c])

    # optional: tighten padding and freeze font size
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    for cell in tbl.get_celld().values():
        cell.PAD = 0.05

def observation_summary_figure(observationID, repro_wd, obsid_5digit, erange, tbin, grating_check, flags):
    fig = plt.figure(figsize=(15, 10), dpi=150)
    gs  = GridSpec(2, 2, height_ratios=[1.5, 1], width_ratios=[1, 2], hspace=0.1, wspace=0.2)

    ax_img = fig.add_subplot(gs[0, 0])
    image_hdul = fits.open(f'{repro_wd}/{obsid_5digit}_bary_2-8keV_cropped.fits')
    image = image_hdul[0].data

    if grating_check == False:
        region_name = 'sgra'
    else:
        region_name = 'order0'
    with open(f'{repro_wd}/{region_name}.reg', "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("ellipse"):
                contents = line[len("ellipse("):-1]
                parts = contents.split(",")
                x, y = float(parts[0]), float(parts[1])

    physical_x = x - 3992.21
    physical_y = y - 3990.74

    cropped_stamp = image[round(physical_y-20):round(physical_y+20), round(physical_x-20):round(physical_x+20)]

    plt.imshow(cropped_stamp+1, origin='lower', norm = LogNorm(vmin=0.9, vmax=cropped_stamp.max()))
    circle = patches.Circle((physical_x - round(physical_x-20)-1, physical_y - round(physical_y-20)-1), 2.54, edgecolor="red", facecolor="none", linewidth=1)
    ax_img.add_patch(circle)
    ax_img.axis('off')

    if len(flags) > 0:
        ax_img.text(0.02, -0.07, f"FLAGS: {flags[0]}", ha='left', va='top', transform=ax_img.transAxes, fontsize=12, fontweight="bold")
    else:
        ax_img.text(0.02, -0.07, "FLAGS: NONE", ha='left', va='top', transform=ax_img.transAxes, fontsize=12, fontweight="bold")


    ax_lc = fig.add_subplot(gs[0, 1])
    bb.plot_bb("./"  + str(observationID) + "/repro/" + "Results/" + str(obsid_5digit) + "_sgra_bayesianBlocks_info.txt", ax_lc) 
    bb.plot_lc("./" +  str(observationID) + "/repro/" + (str(obsid_5digit) + f"_sgra_{erange[0]}-{erange[1]}keV_lc{tbin}_pileup.fits"), 'RATE_PILEUP', 'PILEUP_ERR', ax_lc) 
    ax_lc.set_xlabel("Time (UTC)")
    ax_lc.set_ylabel("Count Rate")
    #ax_lc.plot()

    #ax_lc = fig.add_subplot(gs[0, 1])
    #bb_img = mpimg.imread(f"{repro_wd}/Results/{observationID}_sgra_PLOT.png")
    #ax_lc.imshow(bb_img)
    #ax_lc.axis('off')

    ax_flare = fig.add_subplot(gs[1, :])
    ax_flare.axis('off')
    with open(f'{repro_wd}/Results/{obsid_5digit}_SGRA_TABLE_RESULTS.txt', "r") as f:
        raw = f.read()

    m = re.search(r'^\s*FLARE NUMBER', raw, flags=re.M)
    basic = raw[:m.start()].strip() if m else raw.strip()

    # Remove decorative lines/headers if present
    basic = re.sub(r'^-+\s*$', '', basic, flags=re.M)
    basic = re.sub(r'^\s*BASIC INFORMATION\s*$', '', basic, flags=re.M).strip()

    pairs = []
    for line in basic.splitlines():
        line = line.strip()
        if not line:
            continue
        mm = re.match(r'^([^:]+):\s*(.+)$', line)
        if mm:
            key = mm.group(1).strip()
            val = mm.group(2).strip()
            pairs.append([key, val])

    # Optional: enforce a tidy order if those keys exist
    preferred = [
        "Obs_ID", "Obs Date", "Telescope", "Instrument",
        "Exposure (ks)", "Quiescent Count Rate (10^-3 ct/s)"
    ]
    if pairs:
        # index by key
        d = {k: v for k, v in pairs}
        ordered = [[k, d[k]] for k in preferred if k in d]
        # add any remaining keys not in preferred (in original order)
        seen = set(preferred)
        rest = [[k, v] for k, v in pairs if k not in seen]
        pairs = ordered + rest

    pairs = np.asarray(pairs)

    # Split on lines that start with 'FLARE NUMBER'
    blocks = re.split(r'^\s*(?=FLARE NUMBER\s+\d+)', raw, flags=re.M)
    flares = []
    for b in blocks:
        b = b.strip()
        if not b.startswith("FLARE NUMBER"):
            continue
        d = {}

        # First line: FLARE NUMBER X
        m = re.match(r'^FLARE NUMBER\s+(\d+)', b)
        if m:
            d['FLARE NUMBER'] = m.group(1)

        # Simple "Key: value" lines (stop at empty lines)
        for line in b.splitlines()[1:]:
            line = line.strip()
            if not line:
                continue
            m = re.match(r'^([^:]+):\s*(.+)$', line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip()
                d[key] = val
        flares.append(d)

    right_bottom = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, :], height_ratios=[1, 4], hspace=0.02)
    if len(flares) != 0:
        metrics = ["Start Time", "End Time", "Duration", "Count Rate (mean)", "Count Rate (max)","Energy", "Luminosity", "Flux", "Fluence"]
        rows = []
        for m in metrics:
            row = [m]
            for d in flares:
                row.append(d.get(m, "—"))
            rows.append(row)

        col_labels = ["Metric"] + [f"Flare {d.get('FLARE NUMBER','?')}" for d in flares]

        ax_flare_tbl = fig.add_subplot(right_bottom[1, 0])
        tbl = ax_flare_tbl.table(cellText=rows, colLabels=col_labels, cellLoc="left", loc="center")
        ax_flare_tbl.axis('off')

        for key, cell in tbl.get_celld().items():
            row, col = key
            if row == 0:  # header row
                cell.set_text_props(weight="bold")       # make text bold
                cell.set_facecolor("#d9d9d9")            # light gray background

        # Colour the "Metric" column (col 0, including header and all rows)
        for key, cell in tbl.get_celld().items():
            row, col = key
            if col == 0:
                cell.set_facecolor((5/256, 122/256, 168/256, 0.3))            # light purple
                cell.set_text_props(weight="bold")  

        # Make "Metric" wide, each flare column narrower
        flare_col_widths = {0: 0.2}  # metric column
        for j in range(1, len(col_labels)):      # remaining flare columns
            flare_col_widths[j] = (0.80 / (len(col_labels) - 1))
        set_column_widths(tbl, flare_col_widths)

        tbl.scale(1, 1.2)
    else:
        ax_flare_area = fig.add_subplot(right_bottom[1, 0])
        ax_flare_area.axis('off')
        ax_flare_area.text(0.5, 0.5, "No flares in observation", ha="center", va="center", fontsize=14, fontweight="bold")

    ax_basic_top = fig.add_subplot(right_bottom[0, 0])
    tbl_basic = ax_basic_top.table(cellText=[pairs[:, 1].tolist()], colLabels=pairs[:, 0], cellLoc="left", loc="center")
    ax_basic_top.axis("off")

    for key, cell in tbl_basic.get_celld().items():
        row, col = key
        if row == 0:  # header row
            cell.set_text_props(weight="bold")       # make text bold
            cell.set_facecolor("#d9d9d9")            # light gray background

    # e.g., give the first few columns more space than the tiny numeric ones
    basic_col_widths = {
        0: 0.1,   # Obs_ID
        1: 0.20,   # Obs Date
        2: 0.10,   # Telescope
        3: 0.10,   # Instrument
        4: 0.20,   # Exposure (ks)
        5: 0.30,   # Quiescent Count Rate
    }
    set_column_widths(tbl_basic, basic_col_widths)

    tbl_basic.scale(1, 1.2)

    plt.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.01, wspace=0, hspace=0.25)

    y = 0.93
    fig.lines.append(plt.Line2D([0.01, 0.99], [y, y], transform=fig.transFigure,color="black", linewidth=1.2))
    fig.lines.append(plt.Line2D([0.01, 0.99], [y-0.01, y-0.01], transform=fig.transFigure, color="black", linewidth=1.2))
    
    fig.suptitle(f'Observation {observationID} Summary', fontsize=20)
    fig.savefig(f'{repro_wd}/Results/{observationID}_SUMMARY.png')