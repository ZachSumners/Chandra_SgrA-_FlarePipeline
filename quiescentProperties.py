#!/usr/bin/env python3
"""
scan_sgra_quiescent_rates.py

Recursively scans a directory for files ending in:
    *_SGRA_TABLE_RESULTS.txt
    *_SGRA_TABLE_RESULTS_pileupcorr.txt

From each such file, it extracts the “BASIC INFORMATION” block:
  - Obs_ID
  - Obs Date
  - Telescope
  - Instrument
  - Exposure (ks)
  - Quiescent Count Rate (10^-3 ct/s) ± uncertainty

And writes one CSV row per observation, with columns:
  obs_id, obs_date, telescope, instrument, exposure_ks, quiescent_rate, quiescent_rate_err

Usage:
    python scan_sgra_quiescent_rates.py [--rootdir PATH] [--outfile PATH]

Defaults:
    rootdir = current directory (“.”)
    outfile = "quiescent_rates.csv"
"""

import os
import re
import csv
import argparse

# ----------------------------
# REGEX PATTERNS FOR PARSING “BASIC INFORMATION”
# ----------------------------
RE_OBS_ID       = re.compile(r"^\s*Obs_ID\s*:\s*(\d+)\s*$")
RE_OBS_DATE     = re.compile(r"^\s*Obs Date\s*:\s*([\d\-T:\.]+)\s*$")
RE_TELESCOPE    = re.compile(r"^\s*Telescope\s*:\s*(.+)\s*$")
RE_INSTRUMENT   = re.compile(r"^\s*Instrument\s*:\s*(.+)\s*$")
RE_EXPOSURE     = re.compile(r"^\s*Exposure\s*\(ks\)\s*:\s*([\d\.Ee\+\-]+)\s*$")
# Quiescent Count Rate (10^-3 ct/s):  5.015 +/- 0.263
RE_QUIESCENT    = re.compile(
    r"^\s*Quiescent Count Rate.*:\s*([0-9Ee\.\-]+)\s*\+/-\s*([0-9Ee\.\-]+)\s*$"
)

# Only match files ending in *_SGRA_TABLE_RESULTS.txt or *_SGRA_TABLE_RESULTS_pileupcorr.txt
FILE_PATTERN = re.compile(r".*_SGRA_TABLE_RESULTS(_pileupcorr)?\.txt$")

# ----------------------------
# PARSE ONE FILE'S BASIC INFO
# ----------------------------
def parse_quiescent_block(path):
    """
    Open a single *_SGRA_TABLE_RESULTS*.txt file at `path`, scan line by line,
    and extract:
      - obs_id                (string)
      - obs_date              (string, e.g. "2016-07-18T12:08:57")
      - telescope             (string, e.g. "CHANDRA")
      - instrument            (string, e.g. "ACIS")
      - exposure_ks           (float)
      - quiescent_rate        (float)  in units of 10^-3 ct/s
      - quiescent_rate_err    (float)

    Returns a dict with those keys. If any field is not found, its value will be None.
    """
    info = {
        'obs_id':            None,
        'obs_date':          None,
        'telescope':         None,
        'instrument':        None,
        'exposure_ks':       None,
        'quiescent_rate':    None,
        'quiescent_rate_err':None
    }

    with open(path, 'r') as f:
        for raw in f:
            line = raw.rstrip()

            m = RE_OBS_ID.match(line)
            if m:
                info['obs_id'] = m.group(1)
                continue

            m = RE_OBS_DATE.match(line)
            if m:
                info['obs_date'] = m.group(1)
                continue

            m = RE_TELESCOPE.match(line)
            if m:
                info['telescope'] = m.group(1)
                continue

            m = RE_INSTRUMENT.match(line)
            if m:
                info['instrument'] = m.group(1)
                continue

            m = RE_EXPOSURE.match(line)
            if m:
                info['exposure_ks'] = float(m.group(1))
                continue

            m = RE_QUIESCENT.match(line)
            if m:
                info['quiescent_rate']     = float(m.group(1))
                info['quiescent_rate_err'] = float(m.group(2))
                continue

            # Stop scanning once we've collected everything
            if (info['obs_id'] is not None and
                info['obs_date'] is not None and
                info['telescope'] is not None and
                info['instrument'] is not None and
                info['exposure_ks'] is not None and
                info['quiescent_rate'] is not None and
                info['quiescent_rate_err'] is not None):
                break

    return info

# ----------------------------
# MAIN ROUTINE
# ----------------------------
def main(rootdir, outfile):
    """
    Walk `rootdir` recursively, find all *_SGRA_TABLE_RESULTS*.txt files,
    extract BASIC INFORMATION (quiescent rate, etc.), and write one CSV row
    per observation.
    """
    fieldnames = [
        'obs_id',
        'obs_date',
        'telescope',
        'instrument',
        'exposure_ks',
        'quiescent_rate',
        'quiescent_rate_err'
    ]

    with open(outfile, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        total_files = 0
        total_written = 0

        print(f"Scanning '{rootdir}' for *_SGRA_TABLE_RESULTS*.txt files…\n")

        for dirpath, dirnames, filenames in os.walk(rootdir):
            for fn in filenames:
                if not FILE_PATTERN.match(fn):
                    continue

                total_files += 1
                fullpath = os.path.join(dirpath, fn)

                info = parse_quiescent_block(fullpath)

                # If we did not find Obs_ID or quiescent_rate, still write the row,
                # but those fields will be empty (None → written as blank).
                row = {
                    'obs_id':            info['obs_id'] or '',
                    'obs_date':          info['obs_date'] or '',
                    'telescope':         info['telescope'] or '',
                    'instrument':        info['instrument'] or '',
                    'exposure_ks':       (f"{info['exposure_ks']:.6f}" if info['exposure_ks'] is not None else ''),
                    'quiescent_rate':    (f"{info['quiescent_rate']:.6f}" if info['quiescent_rate'] is not None else ''),
                    'quiescent_rate_err':(f"{info['quiescent_rate_err']:.6f}" if info['quiescent_rate_err'] is not None else '')
                }
                writer.writerow(row)
                total_written += 1

        print(f"Total matching files scanned:  {total_files}")
        print(f"Total rows (observations) written: {total_written}")
        print(f"\nQuiescent rates have been saved to '{outfile}'.\n")


# ----------------------------
# ENTRY POINT
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract quiescent count rate (and other basic info) from *_SGRA_TABLE_RESULTS*.txt files."
    )
    parser.add_argument(
        "--rootdir",
        default=".",
        help="Top‐level directory to scan (default: current directory)"
    )
    parser.add_argument(
        "--outfile",
        default="quiescent_rates.csv",
        help="Name of CSV file to write (default: quiescent_rates.csv)"
    )
    args = parser.parse_args()

    main(args.rootdir, args.outfile)