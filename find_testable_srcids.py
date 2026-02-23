#!/usr/bin/env python3
"""
find_testable_srcids.py — Find SRCIDs that actually have spectra on disk.

Usage (on bell):
    python3 find_testable_srcids.py \
        /mnt/xmmcat/5XMM_data/Spectra \
        /home/mcoriat/Work/XMM/5XMM/automated_fits/5xmm_matched_for_pipeline.fits \
        --subdir product \
        --n_test 20 \
        --output test_srcids.txt

This script:
  1. Scans the data directory for OBS_IDs that contain SRSPEC files.
  2. Extracts the hex SRC_NUM from the spectrum filenames.
  3. Cross-references with the matched catalog to find the
     corresponding SRCIDs.
  4. Writes a file of N SRCIDs known to have spectra on disk.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np
from astropy.io import fits


def find_obsids_with_spectra(data_dir, subdir='product',
                             max_obsids=500):
    """
    Scan the data directory for OBS_IDs that contain at least one
    SRSPEC file.  Returns a dict:
      {obsid: [(hex_srcnum, full_path), ...]}
    """
    obsid_spectra = {}
    n_scanned = 0

    # List all subdirectories in data_dir (each is an OBS_ID)
    try:
        all_dirs = sorted(os.listdir(data_dir))
    except FileNotFoundError:
        print(f"ERROR: Data directory not found: {data_dir}")
        sys.exit(1)

    for obsid_dir in all_dirs:
        spec_path = os.path.join(data_dir, obsid_dir, subdir)
        if not os.path.isdir(spec_path):
            continue

        n_scanned += 1
        try:
            files = os.listdir(spec_path)
        except PermissionError:
            continue

        srspec_files = [f for f in files if 'SRSPEC' in f]
        if srspec_files:
            entries = []
            for f in srspec_files:
                # Extract hex SRC_NUM from filename pattern
                # e.g. P0677180601M1S002SRSPEC0026.FTZ
                # The 4 hex chars after SRSPEC are the SRC_NUM
                m = re.search(r'SRSPEC([0-9A-Fa-f]{4})', f)
                if m:
                    hex_srcnum = m.group(1)
                    entries.append(
                        (hex_srcnum, os.path.join(spec_path, f)))
            if entries:
                obsid_spectra[obsid_dir] = entries

        if len(obsid_spectra) >= max_obsids:
            break

    print(f"  Scanned {n_scanned} OBS_ID directories")
    print(f"  Found {len(obsid_spectra)} OBS_IDs with SRSPEC files")
    return obsid_spectra


def match_with_catalog(obsid_spectra, catalog_file):
    """
    Cross-reference the spectra found on disk with the matched
    catalog to find corresponding SRCIDs.

    Returns a list of (srcid, obsid, srcnum, n_spectra_on_disk).
    """
    print(f"\n  Loading catalog: {catalog_file}")
    with fits.open(catalog_file) as hdul:
        data = hdul[1].data
        cat_srcids = data['SRCID']
        cat_obsids = data['OBS_ID']
        cat_srcnums = data['SRC_NUM']

    print(f"  Catalog has {len(cat_srcids)} rows")

    # Build a lookup: (obsid, srcnum_int) -> srcid
    cat_lookup = {}
    for i in range(len(cat_srcids)):
        obsid = str(cat_obsids[i]).strip()
        srcnum = int(cat_srcnums[i])
        srcid = int(cat_srcids[i])
        key = (obsid, srcnum)
        # Keep the first match (catalog may have duplicates for
        # different GroupIDs)
        if key not in cat_lookup:
            cat_lookup[key] = srcid

    print(f"  Built lookup with {len(cat_lookup)} "
          f"(OBS_ID, SRC_NUM) entries")

    # Cross-reference
    matched = []
    srcid_set = set()
    for obsid, entries in obsid_spectra.items():
        for hex_srcnum, filepath in entries:
            srcnum_int = int(hex_srcnum, 16)
            key = (obsid, srcnum_int)
            if key in cat_lookup:
                srcid = cat_lookup[key]
                if srcid not in srcid_set:
                    srcid_set.add(srcid)
                    matched.append(
                        (srcid, obsid, srcnum_int,
                         len(entries)))

    print(f"  Matched {len(matched)} unique SRCIDs with "
          f"spectra on disk")
    return matched


def main():
    parser = argparse.ArgumentParser(
        description="Find SRCIDs with spectra on disk")
    parser.add_argument('data_dir',
                        help="Top data directory with OBS_IDs")
    parser.add_argument('catalog',
                        help="Matched catalog FITS file")
    parser.add_argument('--subdir', default='product',
                        help="Subdirectory name (default: product)")
    parser.add_argument('--n_test', type=int, default=20,
                        help="Number of SRCIDs to output")
    parser.add_argument('--output', default='test_srcids.txt',
                        help="Output file for SRCIDs")
    parser.add_argument('--max_obsids', type=int, default=500,
                        help="Max OBS_IDs to scan (for speed)")
    args = parser.parse_args()

    print("=" * 50)
    print(" Finding SRCIDs with spectra on disk")
    print("=" * 50)

    # Step 1: Find OBS_IDs with spectra
    print(f"\n Step 1: Scanning {args.data_dir} "
          f"(subdir={args.subdir})...")
    obsid_spectra = find_obsids_with_spectra(
        args.data_dir, args.subdir, args.max_obsids)

    if not obsid_spectra:
        print("\n  ERROR: No SRSPEC files found. Check:")
        print(f"    - Data dir: {args.data_dir}")
        print(f"    - Subdir: {args.subdir}")
        print("    - Are spectra extracted yet?")
        sys.exit(1)

    # Show a few examples
    print("\n  Example spectra found:")
    for obsid in list(obsid_spectra.keys())[:3]:
        entries = obsid_spectra[obsid]
        print(f"    OBS_ID={obsid}: "
              f"{len(entries)} SRSPEC files")
        for hex_s, path in entries[:2]:
            print(f"      SRC_NUM=0x{hex_s} "
                  f"(={int(hex_s, 16)}): "
                  f"{os.path.basename(path)}")

    # Step 2: Cross-reference with catalog
    print(f"\n Step 2: Cross-referencing with catalog...")
    matched = match_with_catalog(
        obsid_spectra, args.catalog)

    if not matched:
        print("\n  ERROR: No matches found between disk spectra "
              "and catalog.")
        print("  This might mean:")
        print("    - Catalog OBS_IDs don't match disk OBS_IDs")
        print("    - SRC_NUM encoding mismatch")
        sys.exit(1)

    # Step 3: Select N sources for testing
    # Prefer sources with multiple spectra for more interesting
    # fits
    matched.sort(key=lambda x: -x[3])  # most spectra first
    selected = matched[:args.n_test]

    print(f"\n Step 3: Selected {len(selected)} SRCIDs for "
          f"testing:")
    for srcid, obsid, srcnum, n_spec in selected:
        print(f"    SRCID={srcid}  "
              f"(OBS_ID={obsid}, SRC_NUM={srcnum}, "
              f"{n_spec} spectra in that OBS_ID)")

    # Step 4: Write output file
    with open(args.output, 'w') as f:
        for srcid, _, _, _ in selected:
            f.write(f"{srcid}\n")

    print(f"\n  Written {len(selected)} SRCIDs to {args.output}")
    print(f"  Use with: bash run_pipeline.sh --test "
          f"--srcid_file {args.output}")
    print("=" * 50)


if __name__ == "__main__":
    main()
