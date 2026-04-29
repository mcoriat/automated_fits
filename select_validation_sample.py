#!/usr/bin/env python
"""
select_validation_sample.py
============================

Select a stratified random sample of SRCIDs from the existing
prefit=True production catalogue for the validation experiment
(prefit vs no_prefit comparison).

The sampling is stratified by signal-to-noise ratio (SNR) so that
the validation covers the regime where LM is most likely to fail
(low SNR sources), which is where Francisco's nH histogram peak
at 5e20 cm^-2 is most likely to come from.

Usage
-----
python select_validation_sample.py \\
    --catalog fit_results_all.fits \\
    --n_per_bin 200 \\
    --output validation_srcids.txt

This writes one SRCID per line, ready to feed `automated_fits.py
--srcid_file`.

The script also writes a metadata file alongside the SRCID list,
containing the SNR bin assignment and original prefit-run
parameter values for every selected source — used later by
compare_prefit_validation.py.
"""

import argparse
import os
import sys
import numpy as np
from astropy.table import Table


def main():
    parser = argparse.ArgumentParser(
        description="Select a stratified validation sample")
    parser.add_argument(
        "--catalog", required=True,
        help="FITS catalogue from production run (prefit=True). "
             "Must contain SRCID and a SNR-like column "
             "(default name: SNR or NET_SNR).")
    parser.add_argument(
        "--snr_col", default=None,
        help="Name of the SNR column in the catalogue. If "
             "omitted, the script will try common defaults: "
             "SNR, NET_SNR, S_N, snr.")
    parser.add_argument(
        "--n_per_bin", type=int, default=200,
        help="Number of sources per SNR bin (default 200). "
             "Total sample = 3 * n_per_bin (low/mid/high).")
    parser.add_argument(
        "--snr_bins", nargs=2, type=float,
        default=[3.0, 10.0],
        help="SNR cuts separating low/mid/high bins "
             "(default 3.0 10.0)")
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility")
    parser.add_argument(
        "--output", required=True,
        help="Output text file with one SRCID per line")
    parser.add_argument(
        "--metadata", default=None,
        help="Optional FITS file to write the per-source "
             "metadata (SNR bin, original parameter values). "
             "Default: <output>.meta.fits")
    args = parser.parse_args()

    # ---- Load catalogue ----
    print(f"Reading catalogue: {args.catalog}")
    tbl = Table.read(args.catalog)
    print(f"  {len(tbl)} rows, columns: "
          f"{tbl.colnames[:10]}{'...' if len(tbl.colnames)>10 else ''}")

    # ---- Find stratification column ----
    # Priority: user-specified > SNR-like > nH_median > random
    if args.snr_col is None:
        for cand in ("SNR", "NET_SNR", "S_N", "snr", "net_snr"):
            if cand in tbl.colnames:
                args.snr_col = cand
                break
    if args.snr_col is None:
        # Fallback: stratify by nH_median — ideal for the prefit
        # validation because it directly targets the 5e20 pile-up
        for cand in ("nH_median", "NH_MEDIAN", "SPEC_NH_PL",
                      "nH", "NH"):
            if cand in tbl.colnames:
                args.snr_col = cand
                if args.snr_bins == [3.0, 10.0]:
                    # Override default SNR bins with nH-appropriate
                    # bins: [0.03, 0.5] in 10^22 cm^-2 units
                    # gives low (<3e20), mid (3e20-5e21),
                    # high (>5e21)
                    args.snr_bins = [0.03, 0.5]
                    print("  Auto-selected nH stratification "
                          f"with bins {args.snr_bins}")
                break
    if args.snr_col is None or args.snr_col not in tbl.colnames:
        # Last resort: no stratification, uniform random sample
        print("  WARNING: No suitable column for stratification. "
              "Falling back to uniform random sampling.")
        args.snr_col = None
    if args.snr_col is not None:
        print(f"  Stratifying on: {args.snr_col}, "
              f"bins: {args.snr_bins}")

    if "SRCID" not in tbl.colnames:
        print("ERROR: catalogue has no SRCID column")
        sys.exit(1)

    # ---- Stratified or random selection ----
    rng = np.random.default_rng(args.seed)
    total_n = 3 * args.n_per_bin  # target sample size

    if args.snr_col is not None:
        # Drop rows with missing/invalid values
        strat_vals = np.asarray(tbl[args.snr_col], dtype=float)
        valid = np.isfinite(strat_vals) & (strat_vals > 0)
        print(f"  Rows with valid {args.snr_col}: "
              f"{valid.sum()} / {len(tbl)}")
        tbl = tbl[valid]
        strat_vals = strat_vals[valid]

        # Stratify
        snr_lo, snr_hi = args.snr_bins
        bins = {
            "low":  (strat_vals < snr_lo),
            "mid":  (strat_vals >= snr_lo)
                    & (strat_vals < snr_hi),
            "high": (strat_vals >= snr_hi),
        }
        print(f"\n{args.snr_col} bins "
              f"(cuts {snr_lo} / {snr_hi}):")
        for name, mask in bins.items():
            print(f"  {name:>5}: {mask.sum():>8} sources")

        selected_idx = []
        bin_label = []
        for name, mask in bins.items():
            idx = np.where(mask)[0]
            if len(idx) == 0:
                print(f"  WARNING: bin {name} is empty, "
                      "skipping")
                continue
            n_sel = min(args.n_per_bin, len(idx))
            chosen = rng.choice(idx, size=n_sel, replace=False)
            selected_idx.extend(chosen.tolist())
            bin_label.extend([name] * n_sel)
            print(f"  selected {n_sel} from {name}")

        selected_idx = np.array(selected_idx)
        bin_label = np.array(bin_label)
    else:
        # Uniform random sampling (no stratification)
        n_sel = min(total_n, len(tbl))
        selected_idx = rng.choice(
            len(tbl), size=n_sel, replace=False)
        bin_label = np.array(["random"] * n_sel)
        print(f"\nRandom sample: {n_sel} sources "
              f"(from {len(tbl)})")

    sub = tbl[selected_idx]

    # ---- Write SRCID list ----
    with open(args.output, 'w') as f:
        for sid in sub["SRCID"]:
            f.write(f"{int(sid)}\n")
    print(f"\nWrote {len(sub)} SRCIDs to {args.output}")

    # ---- Write metadata FITS for later comparison ----
    meta_out = args.metadata or (args.output + ".meta.fits")
    sub_meta = sub.copy()
    sub_meta["VALIDATION_BIN"] = bin_label
    if args.snr_col is not None:
        sub_meta["VALIDATION_STRAT"] = strat_vals[selected_idx]
    else:
        sub_meta["VALIDATION_STRAT"] = np.zeros(len(sub))
    sub_meta.write(meta_out, overwrite=True)
    print(f"Wrote metadata to {meta_out}")
    print("\nNext: run on bell with `--no_prefit` and the same "
          "model, then compare with compare_prefit_validation.py")


if __name__ == "__main__":
    main()
