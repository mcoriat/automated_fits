#!/usr/bin/env python3
"""
backfill_instruments.py — Add 'instruments' and 'n_spectra' columns
to FITS result files produced before these columns were introduced.

Parses the chunk log files to recover the exact instrument order
used during each fit, then writes the columns into the FITS tables.

Usage:
    python3 backfill_instruments.py /data/scratch/pipeline_output

    # Dry run (show what would be done, don't modify files):
    python3 backfill_instruments.py /data/scratch/pipeline_output --dry-run

    # Also backfill the merged file:
    python3 backfill_instruments.py /data/scratch/pipeline_output --include-merged
"""

import argparse
import os
import re
import sys
from collections import defaultdict

from astropy.table import Table


def parse_logs_for_instruments(log_dir):
    """Parse chunk log files to extract instrument order per SRCID.

    Looks for log patterns like:
        Performing simultaneous BXA fit for source 123456 using 3 spectra:
           Spectrum file: /path/to/file.FTZ (pn)
           Spectrum file: /path/to/file.FTZ (MOS1)
           Spectrum file: /path/to/file.FTZ (MOS2)

    Also handles single-spectrum fits:
        Starting BXA fit on 1 spectrum file(s)

    Returns:
        dict: {srcid (int): {"instruments": "pn,MOS1,MOS2",
                              "n_spectra": 3}}
    """
    srcid_info = {}

    # Pattern for "Performing simultaneous BXA fit for source XXXXX using N spectra:"
    re_simul = re.compile(
        r"Performing simultaneous BXA fit for source (\d+) "
        r"using (\d+) spectra")
    # Pattern for "Spectrum file: /path/to/file (INSTRUMENT)"
    re_specfile = re.compile(
        r"Spectrum file:\s+\S+\s+\((\w+)\)")
    # Pattern for single-spectrum: "[N/M] Processing SRCID XXXXX"
    re_processing = re.compile(
        r"\] Processing SRCID (\d+)")
    # Pattern for "Starting BXA fit on 1 spectrum file(s)"
    re_single = re.compile(
        r"Starting BXA fit on (\d+) spectrum file")
    # Instrument encoded in spectrum filename:
    # P<obsid>PNS... → pn, P<obsid>M1S... → MOS1, P<obsid>M2S... → MOS2
    re_inst_from_file = re.compile(
        r"P\d{10}(PN|M1|M2)")

    log_files = sorted([
        os.path.join(log_dir, f) for f in os.listdir(log_dir)
        if f.startswith("chunk_") and f.endswith(".log")
    ])

    if not log_files:
        print(f"  WARNING: No chunk log files found in {log_dir}")
        return srcid_info

    print(f"  Parsing {len(log_files)} log files...")

    for log_file in log_files:
        current_srcid = None
        expected_nspec = 0
        collecting_instruments = []

        with open(log_file, "r", errors="replace") as fh:
            for line in fh:
                # Track current SRCID being processed
                m_proc = re_processing.search(line)
                if m_proc:
                    # Save previous source if we were collecting
                    if (current_srcid is not None
                            and collecting_instruments):
                        _save_instruments(
                            srcid_info, current_srcid,
                            collecting_instruments)
                    current_srcid = int(m_proc.group(1))
                    expected_nspec = 0
                    collecting_instruments = []

                # Multi-spectrum fit header
                m_simul = re_simul.search(line)
                if m_simul:
                    current_srcid = int(m_simul.group(1))
                    expected_nspec = int(m_simul.group(2))
                    collecting_instruments = []
                    continue

                # Spectrum file line with instrument in parens
                m_spec = re_specfile.search(line)
                if m_spec and current_srcid is not None:
                    collecting_instruments.append(m_spec.group(1))
                    continue

                # Single-spectrum fit (no "Spectrum file" lines
                # follow with instrument in parens — derive from
                # the spectrum filename in the log)
                m_single = re_single.search(line)
                if (m_single and current_srcid is not None
                        and int(m_single.group(1)) == 1):
                    # Single spectrum — instrument not always
                    # logged with parens. Look for it in recent
                    # context or mark as single.
                    if not collecting_instruments:
                        expected_nspec = 1

                # Fit completed — save what we have
                if ("Fit completed successfully" in line
                        and current_srcid is not None):
                    if collecting_instruments:
                        _save_instruments(
                            srcid_info, current_srcid,
                            collecting_instruments)
                    elif expected_nspec == 1:
                        # Single spectrum: no IIN, record as
                        # single with instrument unknown from
                        # this log line. We'll try to infer it
                        # from the spectrum filename logged
                        # earlier.
                        pass
                    current_srcid = None
                    collecting_instruments = []

        # Handle last source in the log file
        if current_srcid is not None and collecting_instruments:
            _save_instruments(
                srcid_info, current_srcid,
                collecting_instruments)

    # Second pass: for single-spectrum sources without instrument
    # info, try to extract from "Loading spectrum:" or similar log
    # lines. For now, mark them with n_spectra=1.
    _backfill_single_spectrum_from_logs(srcid_info, log_files)

    print(f"  Found instrument info for {len(srcid_info)} sources")
    return srcid_info


def _save_instruments(srcid_info, srcid, instruments):
    """Save the instrument list for a SRCID."""
    srcid_info[srcid] = {
        "instruments": ",".join(instruments),
        "n_spectra": len(instruments),
    }


def _backfill_single_spectrum_from_logs(srcid_info, log_files):
    """For sources with n_spectra=1 and no instrument info,
    scan logs for the spectrum filename to infer the instrument.

    Spectrum filenames encode the instrument:
        P<10-digit obsid>PNS... → pn
        P<10-digit obsid>M1S... → MOS1
        P<10-digit obsid>M2S... → MOS2
    """
    re_processing = re.compile(r"\] Processing SRCID (\d+)")
    re_bxa_start = re.compile(
        r"Starting BXA fit on 1 spectrum file")
    re_specname = re.compile(
        r"(P\d{10})(PN|M1|M2)\S*\.(FTZ|fits)", re.IGNORECASE)
    inst_map = {"PN": "pn", "M1": "MOS1", "M2": "MOS2"}

    for log_file in log_files:
        current_srcid = None
        with open(log_file, "r", errors="replace") as fh:
            for line in fh:
                m = re_processing.search(line)
                if m:
                    current_srcid = int(m.group(1))

                if (current_srcid is not None
                        and current_srcid not in srcid_info):
                    # Try to find the instrument from the
                    # spectrum filename in this line
                    m_spec = re_specname.search(line)
                    if m_spec:
                        inst = inst_map.get(
                            m_spec.group(2).upper(), "unknown")
                        srcid_info[current_srcid] = {
                            "instruments": inst,
                            "n_spectra": 1,
                        }


def backfill_fits_file(fits_path, srcid_info, dry_run=False):
    """Add instruments and n_spectra columns to a FITS file.

    Parameters
    ----------
    fits_path : str
        Path to the FITS result file.
    srcid_info : dict
        {srcid: {"instruments": str, "n_spectra": int}}
    dry_run : bool
        If True, print what would be done but don't modify files.

    Returns
    -------
    int
        Number of rows updated.
    """
    if not os.path.isfile(fits_path):
        print(f"  SKIP: {fits_path} not found")
        return 0

    t = Table.read(fits_path)
    n_rows = len(t)

    if n_rows == 0:
        print(f"  SKIP: {fits_path} is empty")
        return 0

    # Check if columns already exist and are populated
    if "instruments" in t.colnames:
        n_filled = sum(1 for v in t["instruments"]
                       if v and str(v).strip())
        if n_filled == n_rows:
            print(f"  SKIP: {os.path.basename(fits_path)} "
                  f"already has instruments column "
                  f"({n_filled}/{n_rows} filled)")
            return 0
        else:
            print(f"  {os.path.basename(fits_path)}: "
                  f"instruments column exists but only "
                  f"{n_filled}/{n_rows} filled — backfilling")

    # Build new columns
    instruments_col = []
    n_spectra_col = []
    n_matched = 0

    for row in t:
        srcid = int(row["SRCID"])
        info = srcid_info.get(srcid)
        if info is not None:
            instruments_col.append(info["instruments"])
            n_spectra_col.append(info["n_spectra"])
            n_matched += 1
        else:
            instruments_col.append("")
            n_spectra_col.append(0)

    basename = os.path.basename(fits_path)
    pct = 100 * n_matched / n_rows if n_rows > 0 else 0
    print(f"  {basename}: {n_rows} rows, "
          f"{n_matched} matched ({pct:.0f}%)")

    if dry_run:
        print(f"    DRY RUN — would write to {fits_path}")
        return n_matched

    # Add or replace columns
    t["instruments"] = instruments_col
    t["n_spectra"] = n_spectra_col

    # Atomic write
    tmp_path = fits_path + ".backfill_tmp"
    t.write(tmp_path, format="fits", overwrite=True)
    os.replace(tmp_path, fits_path)

    return n_matched


def main():
    parser = argparse.ArgumentParser(
        description="Backfill instruments/n_spectra columns "
        "in pipeline FITS results using chunk log files.")
    parser.add_argument(
        "output_dir",
        help="Pipeline output directory "
             "(contains chunk_logs/ and fit_results_chunk_*.fits)")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show what would be done without modifying files")
    parser.add_argument(
        "--include-merged", action="store_true",
        help="Also backfill fit_results_all.fits (the merged file)")
    args = parser.parse_args()

    output_dir = args.output_dir
    log_dir = os.path.join(output_dir, "chunk_logs")

    if not os.path.isdir(log_dir):
        print(f"ERROR: Log directory not found: {log_dir}")
        sys.exit(1)

    # Step 1: Parse all logs
    print("Step 1: Parsing chunk logs for instrument info...")
    srcid_info = parse_logs_for_instruments(log_dir)

    if not srcid_info:
        print("ERROR: No instrument info found in logs.")
        sys.exit(1)

    # Step 2: Backfill each chunk FITS file
    print(f"\nStep 2: Backfilling FITS files"
          f"{' (dry run)' if args.dry_run else ''}...")
    chunk_files = sorted([
        os.path.join(output_dir, f) for f in os.listdir(output_dir)
        if f.startswith("fit_results_chunk_") and f.endswith(".fits")
    ])

    total_matched = 0
    total_rows = 0
    for fits_path in chunk_files:
        n = backfill_fits_file(fits_path, srcid_info, args.dry_run)
        total_matched += n
        try:
            total_rows += len(Table.read(fits_path))
        except Exception:
            pass

    print(f"\n  Total: {total_matched} rows backfilled "
          f"across {len(chunk_files)} files")

    # Step 3: Optionally backfill the merged file
    if args.include_merged:
        merged_path = os.path.join(output_dir, "fit_results_all.fits")
        if os.path.isfile(merged_path):
            print(f"\nStep 3: Backfilling merged file...")
            backfill_fits_file(merged_path, srcid_info, args.dry_run)
        else:
            print(f"\nStep 3: Merged file not found "
                  f"(will be created after pipeline finishes)")

    print("\nDone.")


if __name__ == "__main__":
    main()
