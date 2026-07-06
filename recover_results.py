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
import multiprocessing as mp
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

# Pre-fit chi2/dof:  "   Pre-fit completed: chi2/dof = 245.6/253"
PREFIT_RE = re.compile(
    r"Pre-fit completed:\s+chi2/dof\s*=\s*([0-9.eE+-]+)/(\d+)")

# logZ:  "  logZ = -128.6 +- 0.1658"
LOGZ_RE = re.compile(
    r"logZ\s*=\s*([0-9.eE+-]+)\s*\+-\s*([0-9.eE+-]+)")

# Background p-value:  "   Background check: p=0.4316 → OK for ..."
BKG_PVAL_RE = re.compile(
    r"Background check:\s+p=([0-9.eE+-]+)")

# Goodness-of-fit (new pipeline):
#   "   Goodness-of-fit: cstat/dof = 128.3/253, KS p-value = 0.4321"
GOF_RE = re.compile(
    r"Goodness-of-fit:\s+cstat/dof\s*=\s*"
    r"([0-9.eE+-]+)/([0-9-]+),\s*"
    r"KS p-value\s*=\s*([0-9.eE+-]+)")

# BXA timeout/failure
BXA_TIMEOUT_RE = re.compile(r"BXA fit timed out")
BXA_FAIL_RE = re.compile(r"BXA solver failed")


def parse_log_file(log_path):
    """Parse a per-source log file and extract fit results.

    Returns a dict with:
      - "params": OrderedDict {name: (median, p16, p84)} or None
      - "prefit_chi2": float or NaN
      - "prefit_dof": int or -1
      - "prefit_pvalue": float or NaN  (from chi2/dof via scipy)
      - "logZ": float or NaN
      - "logZ_err": float or NaN
      - "bkg_pvalue": float or NaN
      - "flag": int  (0=OK, see flag definitions)
    Returns None only if the log file cannot be read at all.
    """
    params = OrderedDict()
    found_success = False
    prefit_chi2 = np.nan
    prefit_dof = -1
    logZ = np.nan
    logZ_err = np.nan
    bkg_pvalue = np.nan
    gof_cstat = np.nan
    gof_dof = -1
    gof_ks_pvalue = np.nan
    bxa_timed_out = False
    bxa_failed = False

    try:
        with open(log_path, "r", errors="replace") as f:
            for line in f:
                # --- Pre-fit chi2/dof ---
                m = PREFIT_RE.search(line)
                if m:
                    prefit_chi2 = float(m.group(1))
                    prefit_dof = int(m.group(2))
                    continue

                # --- logZ (take the first occurrence) ---
                if np.isnan(logZ):
                    m = LOGZ_RE.search(line)
                    if m:
                        logZ = float(m.group(1))
                        logZ_err = float(m.group(2))
                        continue

                # --- Background p-value ---
                m = BKG_PVAL_RE.search(line)
                if m:
                    bkg_pvalue = float(m.group(1))
                    continue

                # --- BXA timeout/failure ---
                if BXA_TIMEOUT_RE.search(line):
                    bxa_timed_out = True
                    continue
                if BXA_FAIL_RE.search(line):
                    bxa_failed = True
                    continue

                # --- Goodness-of-fit (new pipeline) ---
                m = GOF_RE.search(line)
                if m:
                    gof_cstat = float(m.group(1))
                    gof_dof = int(m.group(2))
                    gof_ks_pvalue = float(m.group(3))
                    continue

                # --- Fit completed + parameter block ---
                if "Fit completed successfully" in line:
                    found_success = True
                    params.clear()
                    continue

                if found_success:
                    m = PARAM_RE.match(line)
                    if m:
                        name = m.group(1)
                        median = float(m.group(2))
                        p16 = float(m.group(3))
                        p84 = float(m.group(4))
                        params[name] = (median, p16, p84)
                    elif line.strip() and not line.strip().startswith(
                            "Cleaned up"):
                        found_success = False

    except Exception as e:
        print(f"  WARNING: Could not read {log_path}: {e}",
              file=sys.stderr)
        return None

    if not params:
        return None

    # Compute pre-fit p-value from chi2/dof
    prefit_pvalue = np.nan
    if np.isfinite(prefit_chi2) and prefit_dof > 0:
        try:
            from scipy.stats import chi2 as chi2_dist
            prefit_pvalue = 1.0 - chi2_dist.cdf(
                prefit_chi2, prefit_dof)
        except Exception:
            pass

    # --- Flag logic ---
    # 0 = no issues
    # 1 = poor goodness-of-fit: KS p-value < 0.01 if available,
    #     otherwise pre-fit p-value < 0.01
    # 2 = PhoIndex pegged at prior boundary
    #     (median within 0.05 of hard limits 1.0 or 3.0)
    # 4 = high nH (median >= 100 in 1e22 units = 1e24 cm^-2,
    #     Compton-thick). Matches NH_PEG_THRESHOLD in
    #     spectral_fitting_bxa_adapted.py. The prior cap is
    #     1000, so this is a physical selector rather than a
    #     true boundary peg. (An earlier version used 9.5,
    #     from when the prior cap was 10.)
    # Flags are combined as a bitmask.
    flag = 0

    # Prefer KS p-value (proper GoF) over pre-fit p-value
    if np.isfinite(gof_ks_pvalue):
        if gof_ks_pvalue < 0.01:
            flag |= 1
    elif np.isfinite(prefit_pvalue) and prefit_pvalue < 0.01:
        flag |= 1

    if "PhoIndex" in params:
        ph_med = params["PhoIndex"][0]
        if ph_med <= 1.05 or ph_med >= 2.95:
            flag |= 2

    if "nH" in params:
        nh_med = params["nH"][0]
        if nh_med >= 100.0:
            flag |= 4

    return {
        "params": params,
        "prefit_chi2": prefit_chi2,
        "prefit_dof": prefit_dof,
        "prefit_pvalue": prefit_pvalue,
        "logZ": logZ,
        "logZ_err": logZ_err,
        "bkg_pvalue": bkg_pvalue,
        "cstat": gof_cstat,
        "cstat_dof": gof_dof,
        "ks_pvalue": gof_ks_pvalue,
        "flag": flag,
    }


def _parse_one_source(args_tuple):
    """Worker function for multiprocessing.

    Takes (srcid_str, log_path) and returns
    (srcid_str, "ok", row_dict) or (srcid_str, "no_log"|"no_fit", None).
    Must be a top-level function for pickling.
    """
    srcid_str, log_path = args_tuple

    if not os.path.exists(log_path):
        return (srcid_str, "no_log", None)

    result = parse_log_file(log_path)
    if result is None:
        return (srcid_str, "no_fit", None)

    row = {"SRCID": int(srcid_str)}
    for pname, (median, p16, p84) in result["params"].items():
        row[f"{pname}_median"] = median
        row[f"{pname}_p16"] = p16
        row[f"{pname}_p84"] = p84

    row["prefit_chi2"] = result["prefit_chi2"]
    row["prefit_dof"] = result["prefit_dof"]
    row["prefit_pvalue"] = result["prefit_pvalue"]
    row["logZ"] = result["logZ"]
    row["logZ_err"] = result["logZ_err"]
    row["bkg_pvalue"] = result["bkg_pvalue"]
    row["cstat"] = result["cstat"]
    row["cstat_dof"] = result["cstat_dof"]
    row["ks_pvalue"] = result["ks_pvalue"]
    row["flag"] = result["flag"]

    return (srcid_str, "ok", row)


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
        "--nworkers", type=int, default=0,
        help="Number of parallel workers (default: all CPUs)")
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

    # Parse log files in parallel
    rows = []
    n_no_log = 0
    n_no_fit = 0
    n_recovered = 0

    nworkers = args.nworkers if args.nworkers > 0 else mp.cpu_count()
    print(f"Parsing log files with {nworkers} workers...")
    t1 = time.time()

    # Build work items: (srcid_str, log_path) tuples
    work_items = [
        (s, os.path.join(output_dir, s, f"{s}{log_suffix}"))
        for s in srcid_dirs
    ]

    n_total = len(work_items)
    n_done = 0

    with mp.Pool(nworkers) as pool:
        for srcid_str, status, row in pool.imap_unordered(
                _parse_one_source, work_items, chunksize=500):
            n_done += 1
            if status == "ok":
                rows.append(row)
                n_recovered += 1
            elif status == "no_log":
                n_no_log += 1
            else:
                n_no_fit += 1

            if n_done % 50000 == 0:
                elapsed = time.time() - t1
                rate = n_done / elapsed
                eta = (n_total - n_done) / rate
                print(f"  [{n_done}/{n_total}] "
                      f"{n_recovered} recovered, "
                      f"{n_no_fit} no fit, "
                      f"{n_no_log} no log  "
                      f"({rate:.0f} src/s, ETA {eta:.0f}s)")

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

    # Integer columns (use specific defaults, not NaN)
    INT_COLS = {"SRCID": (np.int64, 0),
                "prefit_dof": (np.int32, -1),
                "cstat_dof": (np.int32, -1),
                "flag": (np.int16, -1)}

    # Build column arrays
    col_data = {}
    for col in all_cols:
        if col in INT_COLS:
            dtype, default = INT_COLS[col]
            col_data[col] = np.array(
                [r.get(col, default) for r in rows], dtype=dtype)
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

    # Flag summary
    flags = table["flag"]
    print(f"\n  Flag summary (bitmask):")
    print(f"    flag=0 (clean):          "
          f"{np.sum(flags == 0)}")
    print(f"    bit 1 (prefit p<0.01):   "
          f"{np.sum((flags & 1) > 0)}")
    print(f"    bit 2 (PhoIndex pegged): "
          f"{np.sum((flags & 2) > 0)}")
    print(f"    bit 4 (nH >= 1e24 cm^-2): "
          f"{np.sum((flags & 4) > 0)}")


if __name__ == "__main__":
    main()
