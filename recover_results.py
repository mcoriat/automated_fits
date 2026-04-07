#!/usr/bin/env python3
"""
recover_results.py — Recover BXA fit results from per-source log files.

After the pipeline run of 06 Mar 2026, all 30 chunks completed their
fits but crashed at the FITS export step due to a bug in _write_rows_to_fits
(astropy could not auto-detect format from a .tmp extension).  The fit
results (median, p16, p84) were logged to each source's per-source log
file.  This script scans those log files and rebuilds the merged FITS
catalogue.

Usage (on bell):
    python3 recover_results.py /mnt/sas/pipeline_output \
        --model powerlaw \
        --output fit_results_all.fits

This will:
  1. Walk all SRCID subdirectories under the output dir
  2. Parse {SRCID}_process_log_{model}.txt for "Fit completed successfully"
  3. Extract parameter medians, p16, p84
  4. Write a single merged FITS table
"""

import argparse
import os
import re
import sys
import time
from collections import OrderedDict

import numpy as np
from astropy.table import Table


# Regex to match:  "   nH (median [percentiles 16,84]): 3.904e-01 [2.885e-01,5.414e-01]"
PARAM_RE = re.compile(
    r"^\s+(\S+)\s+\(median \[percentiles 16,84\]\):\s+"
    r"([^\s]+)\s+\[([^\s,]+),([^\s\]]+)\]"
)


def parse_log_file(log_path):
    """Parse a per-source log file and extract fit results.

    Returns a dict of {param_name: (median, p16, p84)} if the fit
    succeeded, or None if no successful fit was found.
    """
    params = OrderedDict()
    found_success = False

    try:
        with open(log_path, "r", errors="replace") as f:
            for line in f:
                if "Fit completed successfully" in line:
                    found_success = True
                    params.clear()  # reset in case of multiple fits
                    continue

                if found_success:
                    m = PARAM_RE.match(line)
                    if m:
                        name = m.group(1)
                        median = float(m.group(2))
                        p16 = float(m.group(3))
                        p84 = float(m.group(4))
                        params[name] = (median, p16, p84)
                    elif line.strip() and not line.strip().startswith("Cleaned up"):
                        # End of parameter block
                        found_success = False
    except Exception as e:
        print(f"  WARNING: Could not read {log_path}: {e}",
              file=sys.stderr)
        return None

    return params if params else None


def main():
    parser = argparse.ArgumentParser(
        description="Recover BXA fit results from per-source log files")
    parser.add_argument(
        "output_dir",
        help="Pipeline output directory containing SRCID subdirs")
    parser.add_argument(
        "--model", default="powerlaw",
        help="Model name used in log filenames (default: powerlaw)")
    parser.add_argument(
        "--output", default="fit_results_all.fits",
        help="Output FITS filename (default: fit_results_all.fits)")
    parser.add_argument(
        "--dry_run", action="store_true",
        help="Just count recoverable sources, don't write FITS")
    args = parser.parse_args()

    output_dir = args.output_dir
    model = args.model
    log_suffix = f"_process_log_{model}.txt"

    if not os.path.isdir(output_dir):
        print(f"ERROR: {output_dir} is not a directory")
        sys.exit(1)

    # Scan for SRCID directories (numeric names)
    print(f"Scanning {output_dir} for source directories...")
    t0 = time.time()

    entries = os.listdir(output_dir)
    srcid_dirs = [e for e in entries if e.isdigit()]
    srcid_dirs.sort()

    print(f"  Found {len(srcid_dirs)} source directories "
          f"({time.time() - t0:.1f}s)")

    # Parse log files
    rows = []
    n_no_log = 0
    n_no_fit = 0
    n_recovered = 0

    print(f"Parsing log files...")
    t1 = time.time()

    for i, srcid_str in enumerate(srcid_dirs):
        if (i + 1) % 10000 == 0:
            elapsed = time.time() - t1
            rate = (i + 1) / elapsed
            eta = (len(srcid_dirs) - i - 1) / rate
            print(f"  [{i+1}/{len(srcid_dirs)}] "
                  f"{n_recovered} recovered, "
                  f"{n_no_fit} no fit, "
                  f"{n_no_log} no log  "
                  f"({rate:.0f} src/s, ETA {eta:.0f}s)")

        srcid = int(srcid_str)
        log_path = os.path.join(
            output_dir, srcid_str, f"{srcid_str}{log_suffix}")

        if not os.path.exists(log_path):
            n_no_log += 1
            continue

        params = parse_log_file(log_path)
        if params is None:
            n_no_fit += 1
            continue

        row = {"SRCID": srcid}
        for pname, (median, p16, p84) in params.items():
            row[f"{pname}_median"] = median
            row[f"{pname}_p16"] = p16
            row[f"{pname}_p84"] = p84

        rows.append(row)
        n_recovered += 1

    elapsed = time.time() - t1
    print(f"\nDone parsing ({elapsed:.1f}s)")
    print(f"  Recovered:  {n_recovered}")
    print(f"  No fit:     {n_no_fit}")
    print(f"  No log:     {n_no_log}")
    print(f"  Total dirs: {len(srcid_dirs)}")

    if args.dry_run:
        print("\n[DRY RUN] No FITS file written.")
        return

    if not rows:
        print("\nNo results to write.")
        return

    # Build table
    print(f"\nBuilding FITS table with {len(rows)} rows...")

    # Collect all column names
    all_cols = OrderedDict()
    all_cols["SRCID"] = True
    for row in rows:
        for k in row:
            all_cols[k] = True

    # Build column arrays
    col_data = {}
    for col in all_cols:
        if col == "SRCID":
            col_data[col] = np.array(
                [r.get(col, 0) for r in rows], dtype=np.int64)
        else:
            col_data[col] = np.array(
                [r.get(col, np.nan) for r in rows], dtype=np.float64)

    table = Table(col_data)

    # Write
    fits_path = os.path.join(output_dir, args.output)
    table.write(fits_path, format="fits", overwrite=True)
    print(f"  Written: {fits_path}")
    print(f"  Rows:    {len(table)}")
    print(f"  Columns: {table.colnames}")


if __name__ == "__main__":
    main()
