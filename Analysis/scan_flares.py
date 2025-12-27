#!/usr/bin/env python3
"""
scan_sgra_flare_counts.py

Recursively walks a directory tree looking only for files whose names end in
   • *_SGRA_TABLE_RESULTS.txt
   • *_SGRA_TABLE_RESULTS_pileupcorr.txt

For each such file, it parses:
  • “Obs_ID: ####”                → logs the observation ID
  • “No Flares in Observation”   → counts as zero flares
  • each “FLARE NUMBER X” line   → each occurrence increments flare_count

At the end, prints:
  • Total number of matching SGRA_TABLE_RESULTS*.txt files scanned
  • How many had 0 flares
  • How many had exactly 1 flare
  • How many had >1 flares
  • Grand total number of flares across all files
"""

import os
import re
import argparse

# ----------------------------
# REGEX PATTERNS FOR PARSING
# ----------------------------
RE_OBS_ID       = re.compile(r"^\s*Obs_ID:\s*(\d+)\s*$")
RE_NO_FLARES    = re.compile(r"^\s*No Flares in Observation\s*$", re.IGNORECASE)
RE_FLARE_HEADER = re.compile(r"^\s*FLARE NUMBER\s+(\d+)\s*$")

# Only match filenames ending in:
#    _SGRA_TABLE_RESULTS.txt
#    _SGRA_TABLE_RESULTS_pileupcorr.txt
FILE_PATTERN = re.compile(r".*_SGRA_TABLE_RESULTS(_pileupcorr)?\.txt$")

# ----------------------------
# PARSING FUNCTION
# ----------------------------
def parse_sgra_table_results(path):
    """
    Parse one “*_SGRA_TABLE_RESULTS*.txt” file at `path` and count flares.

    Returns:
      (obs_id_str, num_flares, no_flares_flag)

    - obs_id_str    = the Obs_ID (string) if found (else None)
    - num_flares    = how many lines “FLARE NUMBER X” appear
    - no_flares_flag = True if "No Flares in Observation" appears anywhere
    """
    obs_id = None
    no_flares = False
    flare_count = 0

    with open(path, 'r') as f:
        for raw in f:
            line = raw.rstrip()

            # Capture Obs_ID if present
            m = RE_OBS_ID.match(line)
            if m:
                obs_id = m.group(1)
                continue

            # “No Flares in Observation”
            if RE_NO_FLARES.match(line):
                no_flares = True
                # Override any counted flares
                continue

            # Each “FLARE NUMBER X” line → increment count
            m = RE_FLARE_HEADER.match(line)
            if m:
                flare_count += 1
                continue

    return (obs_id, flare_count, no_flares)


# ----------------------------
# MAIN SCRIPT
# ----------------------------
def main(rootdir):
    total_files = 0
    obs_zero = 0
    obs_one  = 0
    obs_multi= 0
    total_flares = 0
    flare_ids = []

    print(f"Scanning '{rootdir}' for matching SGRA_TABLE_RESULTS*.txt files…\n")

    for dirpath, dirnames, filenames in os.walk(rootdir):
        for fn in filenames:
            # Only pick up filenames matching our SGRA pattern
            if not FILE_PATTERN.match(fn):
                continue

            fullpath = os.path.join(dirpath, fn)
            total_files += 1

            obs_id, flare_count, no_flares_flag = parse_sgra_table_results(fullpath)

            # Categorize based on what we found
            if no_flares_flag:
                obs_zero += 1
                print(f"[Obs {obs_id or 'unknown'}] {fn}: No flares detected.")
            else:
                if flare_count == 0:
                    # Neither “No Flares” nor any “FLARE NUMBER” lines
                    obs_zero += 1
                    print(f"[Obs {obs_id or 'unknown'}] {fn}: 0 flares (no FLARE NUMBER).")
                    
                elif flare_count == 1:
                    obs_one += 1
                    total_flares += 1
                    print(f"[Obs {obs_id or 'unknown'}] {fn}: 1 flare.")
                    flare_ids.append(obs_id)
                else:
                    obs_multi += 1
                    total_flares += flare_count
                    print(f"[Obs {obs_id or 'unknown'}] {fn}: {flare_count} flares.")
                    flare_ids.append(obs_id)

    # Final tally
    print("\n" + "-"*60)
    print(f"Total matching SGRA_TABLE_RESULTS*.txt files scanned: {total_files}")
    print(f"Observations with 0 flares:                            {obs_zero}")
    print(f"Observations with exactly 1 flare:                      {obs_one}")
    print(f"Observations with >1 flares:                            {obs_multi}")
    print(f"--------------------------------------------------------------")
    print(f"Grand total number of flares found:                     {total_flares}")
    print("-"*60 + "\n")

    print(flare_ids, len(flare_ids))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count flares in *_SGRA_TABLE_RESULTS*.txt files."
    )
    parser.add_argument(
        "rootdir",
        nargs="?",
        default=".",
        help="Top‐level directory to scan (default: current directory)"
    )
    args = parser.parse_args()
    main(args.rootdir)