#!/usr/bin/env python3
"""
scan_sgra_flares_to_csv.py

Recursively searches a directory for files named
   *_SGRA_TABLE_RESULTS.txt
   *_SGRA_TABLE_RESULTS_pileupcorr.txt

For each file it finds, it extracts:
  • Obs_ID
  • For each “FLARE NUMBER N” block:
      – start_mjd
      – end_mjd
      – duration_s
      – rate_mean ± rate_mean_err
      – rate_max  ± rate_max_err
      – energy_1e37erg ± energy_err
      – luminosity_1e34 ± luminosity_err
      – flux_1e_minus12 ± flux_err
      – fluence_ct

If an observation has no flares (including “No Flares in Observation”), it writes no rows.

Usage:
    python scan_sgra_flares_to_csv.py [--rootdir PATH] [--outfile PATH]

Defaults:
    rootdir = current directory (“.”)
    outfile = "flare_properties.csv"
"""

import os
import re
import csv
import argparse
import sys

# ----------------------------
# REGEX PATTERNS (adjusted to match "+/-")
# ----------------------------
RE_OBS_ID       = re.compile(r"^\s*Obs_ID:\s*(\d+)\s*$")
RE_NO_FLARES    = re.compile(r"^\s*No Flares in Observation\s*$", re.IGNORECASE)
RE_FLARE_HEADER = re.compile(r"^\s*FLARE NUMBER\s+(\d+)\s*$")

RE_START_TIME   = re.compile(r"^\s*Start Time:\s*([\d\.]+)\s*\(MJD\)\s*$")
RE_END_TIME     = re.compile(r"^\s*End Time:\s*([\d\.]+)\s*\(MJD\)\s*$")
RE_DURATION     = re.compile(r"^\s*Duration:\s*([\d\.]+)\s*\(s\)\s*$")

# Count Rate (mean):  0.0219 +/- 0.0029 (ct/s)
RE_RATE_MEAN    = re.compile(
    r"^\s*Count Rate \(mean\):\s*([0-9Ee\.\-]+)\s*\+/-\s*([0-9Ee\.\-]+)\s*\(ct/s\)\s*$"
)
# Count Rate (max):   0.0445 +/- 0.0134 (ct/s)
RE_RATE_MAX     = re.compile(
    r"^\s*Count Rate \(max\):\s*([0-9Ee\.\-]+)\s*\+/-\s*([0-9Ee\.\-]+)\s*\(ct/s\)\s*$"
)

# Energy:       3.3376 +/- 0.578 10^37 ergs
RE_ENERGY       = re.compile(
    r"^\s*Energy:\s*([0-9Ee\.\-]+)\s*\+/-\s*([0-9Ee\.\-]+)\s*10\^37\s*ergs\s*$"
)
# Luminosity:   1.3019 +/- 0.2255 10^34 erg/s
RE_LUMINOSITY   = re.compile(
    r"^\s*Luminosity:\s*([0-9Ee\.\-]+)\s*\+/-\s*([0-9Ee\.\-]+)\s*10\^34\s*erg/s\s*$"
)
# Flux:         1.726  +/- 0.2989 10^-12 erg/s/cm2
RE_FLUX         = re.compile(
    r"^\s*Flux:\s*([0-9Ee\.\-]+)\s*\+/-\s*([0-9Ee\.\-]+)\s*10\^-\d+\s*erg/s/cm2\s*$"
)
# Fluence:      56.2442 ct
RE_FLUENCE      = re.compile(r"^\s*Fluence:\s*([0-9Ee\.\-]+)\s*ct\s*$")

# Only match files ending in *_SGRA_TABLE_RESULTS.txt or *_SGRA_TABLE_RESULTS_pileupcorr.txt
FILE_PATTERN = re.compile(r".*_SGRA_TABLE_RESULTS(_pileupcorr)?\.txt$")

# ----------------------------
# PARSING FUNCTION
# ----------------------------
def parse_sgra_table_results(path):
    """
    Parse a single SGRA_TABLE_RESULTS*.txt file and return:
      obs_id        (string or None)
      has_no_flares (bool)
      flares        (list of dicts, one per flare)

    Each flare dict contains keys:
      'flare_number'       : int
      'start_mjd'          : float
      'end_mjd'            : float
      'duration_s'         : float
      'rate_mean'          : float
      'rate_mean_err'      : float
      'rate_max'           : float
      'rate_max_err'       : float
      'energy_1e37erg'     : float
      'energy_err'         : float
      'luminosity_1e34'    : float
      'luminosity_err'     : float
      'flux_1e_minus12'    : float
      'flux_err'           : float
      'fluence_ct'         : float

    If "No Flares in Observation" is found, has_no_flares=True.
    If no FLARE NUMBER sections appear, flares will be [].
    """
    obs_id = None
    has_no_flares = False
    flares = []
    current = None

    with open(path, 'r') as f:
        for raw in f:
            line = raw.rstrip()

            # 1) Obs_ID:
            m = RE_OBS_ID.match(line)
            if m:
                obs_id = m.group(1)
                continue

            # 2) "No Flares in Observation"
            if RE_NO_FLARES.match(line):
                has_no_flares = True
                # don’t collect further flare blocks; still scan for Obs_ID if needed
                continue

            # 3) FLARE HEADER
            m = RE_FLARE_HEADER.match(line)
            if m and not has_no_flares:
                # close previous flare
                if current is not None:
                    flares.append(current)

                current = {
                    'flare_number': int(m.group(1)),
                    'start_mjd': None,
                    'end_mjd': None,
                    'duration_s': None,
                    'rate_mean': None,
                    'rate_mean_err': None,
                    'rate_max': None,
                    'rate_max_err': None,
                    'energy_1e37erg': None,
                    'energy_err': None,
                    'luminosity_1e34': None,
                    'luminosity_err': None,
                    'flux_1e_minus12': None,
                    'flux_err': None,
                    'fluence_ct': None
                }
                continue

            # 4) inside a flare block, fill fields
            if current is not None and not has_no_flares:
                m = RE_START_TIME.match(line)
                if m:
                    current['start_mjd'] = float(m.group(1))
                    continue

                m = RE_END_TIME.match(line)
                if m:
                    current['end_mjd'] = float(m.group(1))
                    continue

                m = RE_DURATION.match(line)
                if m:
                    current['duration_s'] = float(m.group(1))
                    continue

                m = RE_RATE_MEAN.match(line)
                if m:
                    current['rate_mean'] = float(m.group(1))
                    current['rate_mean_err'] = float(m.group(2))
                    continue

                m = RE_RATE_MAX.match(line)
                if m:
                    current['rate_max'] = float(m.group(1))
                    current['rate_max_err'] = float(m.group(2))
                    continue

                m = RE_ENERGY.match(line)
                if m:
                    current['energy_1e37erg'] = float(m.group(1))
                    current['energy_err'] = float(m.group(2))
                    continue

                m = RE_LUMINOSITY.match(line)
                if m:
                    current['luminosity_1e34'] = float(m.group(1))
                    current['luminosity_err'] = float(m.group(2))
                    continue

                m = RE_FLUX.match(line)
                if m:
                    current['flux_1e_minus12'] = float(m.group(1))
                    current['flux_err'] = float(m.group(2))
                    continue

                m = RE_FLUENCE.match(line)
                if m:
                    current['fluence_ct'] = float(m.group(1))
                    continue

    # append last flare if present
    if current is not None and not has_no_flares:
        flares.append(current)

    return obs_id, has_no_flares, flares


# ----------------------------
# MAIN SCRIPT
# ----------------------------
def main(rootdir, outfile):
    """
    Scan for *_SGRA_TABLE_RESULTS*.txt under `rootdir` and write one CSV row per real flare.
    Skip any file that has no flares (including “No Flares in Observation”).
    """
    fieldnames = [
        'obs_id',
        'filename',
        'flare_number',
        'start_mjd',
        'end_mjd',
        'duration_s',
        'rate_mean',
        'rate_mean_err',
        'rate_max',
        'rate_max_err',
        'energy_1e37erg',
        'energy_err',
        'luminosity_1e34',
        'luminosity_err',
        'flux_1e_minus12',
        'flux_err',
        'fluence_ct'
    ]

    with open(outfile, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        total_files = 0
        skipped_no_flares = 0
        total_obs_with_flares = 0
        total_flares = 0

        print(f"Scanning '{rootdir}' for *_SGRA_TABLE_RESULTS*.txt files…\n")

        for dirpath, dirnames, filenames in os.walk(rootdir):
            for fn in filenames:
                if not FILE_PATTERN.match(fn):
                    continue

                total_files += 1
                fullpath = os.path.join(dirpath, fn)

                obs_id, has_no_flares, flares = parse_sgra_table_results(fullpath)
                print(obs_id)

                # Skip if no flares
                if has_no_flares or not flares:
                    skipped_no_flares += 1
                    continue

                # Otherwise write one CSV row per flare
                total_obs_with_flares += 1
                for fl in flares:
                    row = {
                        'obs_id': obs_id or '',
                        'filename': fn,
                        'flare_number': fl.get('flare_number', ''),
                        'start_mjd': fl.get('start_mjd', ''),
                        'end_mjd': fl.get('end_mjd', ''),
                        'duration_s': fl.get('duration_s', ''),
                        'rate_mean': fl.get('rate_mean', ''),
                        'rate_mean_err': fl.get('rate_mean_err', ''),
                        'rate_max': fl.get('rate_max', ''),
                        'rate_max_err': fl.get('rate_max_err', ''),
                        'energy_1e37erg': fl.get('energy_1e37erg', ''),
                        'energy_err': fl.get('energy_err', ''),
                        'luminosity_1e34': fl.get('luminosity_1e34', ''),
                        'luminosity_err': fl.get('luminosity_err', ''),
                        'flux_1e_minus12': fl.get('flux_1e_minus12', ''),
                        'flux_err': fl.get('flux_err', ''),
                        'fluence_ct': fl.get('fluence_ct', '')
                    }
                    writer.writerow(row)
                    total_flares += 1

        # Summary
        print(f"Total matching files scanned:          {total_files}")
        print(f"Observations skipped (no flares):      {skipped_no_flares}")
        print(f"Observations with ≥1 flare(s):         {total_obs_with_flares}")
        print(f"Grand total flares written to CSV:     {total_flares}")
        print(f"\nAll flare properties have been saved to '{outfile}'.\n")


# ----------------------------
# ENTRY POINT
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract flares from *_SGRA_TABLE_RESULTS*.txt files into a CSV."
    )
    parser.add_argument(
        "--rootdir",
        default=".",
        help="Top‐level directory to scan (default: current directory)"
    )
    parser.add_argument(
        "--outfile",
        default="flare_properties_dec2_test.csv",
        help="CSV file to write (default: flare_properties.csv)"
    )
    args = parser.parse_args()

    main(args.rootdir, args.outfile)