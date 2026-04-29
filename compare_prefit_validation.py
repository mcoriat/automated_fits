#!/usr/bin/env python
"""
compare_prefit_validation.py
============================

Compare two BXA fit catalogues for the same set of SRCIDs:

  catalog A: production run with prefit=True   (current default)
  catalog B: validation run with --no_prefit   (wide priors)

Produces:
  - histograms of nH, PhoIndex, lg10Flux for A vs B,
    side-by-side and stratified by validation bin (from the metadata
    file written by select_validation_sample.py)
  - per-source scatter plots (A median vs B median) with the
    1:1 line and ±1 sigma error bars
  - summary statistics: median offset, KS-test p-value,
    Mann-Whitney U-test p-value
  - a quantitative answer to: "does the 5e20 cm^-2 peak in A
    disappear in B?" (counts of sources with nH within ±0.05
    dex of the LM starting value 0.05x10^22 cm^-2)
  - a JSON summary file with the headline numbers, easy to
    paste into an email reply to Francisco/Giorgos

Usage
-----
python compare_prefit_validation.py \\
    --prefit fit_results_prefit.fits \\
    --no_prefit fit_results_no_prefit.fits \\
    --metadata validation_srcids.txt.meta.fits \\
    --output_dir validation_plots/

Required catalogue columns (or close equivalents — see
PARAM_CANDIDATES below):
  SRCID, NH_MEDIAN, PHOINDEX_MEDIAN, LG10FLUX_MEDIAN
  (and optionally _LO/_HI for percentiles)
"""

import argparse
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, join
from scipy.stats import ks_2samp, mannwhitneyu


# Map "logical" parameter names to candidate column names in
# the production FITS tables. The first match in each list is
# used.
PARAM_CANDIDATES = {
    "NH":       ["NH_MEDIAN", "SPEC_NH_PL", "nH_median",
                 "NH", "nH"],
    "PHOINDEX": ["PHOINDEX_MEDIAN", "SPEC_PHOINDEX_PL",
                 "PhoIndex_median", "PHOINDEX", "PhoIndex"],
    "LG10FLUX": ["LG10FLUX_MEDIAN", "SPEC_LG10FLUX_PL",
                 "lg10Flux_median", "LG10FLUX", "lg10Flux",
                 "log10Flux"],
}

# LM starting values (from spectral_fitting_bxa_adapted.py).
# A bias-detection metric: count sources with NH within
# ±NH_STARTING_VALUE_TOL dex of the starting value 5e20 cm^-2
# (= 0.05 in 10^22 units).
NH_STARTING_VALUE = 0.05      # 10^22 cm^-2 = 5e20 cm^-2
NH_STARTING_VALUE_TOL = 0.05  # dex window around 5e20


def find_col(tbl, candidates, label):
    """Return the first candidate column present in `tbl`, or
    raise ValueError."""
    for c in candidates:
        if c in tbl.colnames:
            return c
    raise ValueError(
        f"None of the candidate columns for {label} were "
        f"found. Tried: {candidates}. Available: "
        f"{tbl.colnames}")


def joint_table(prefit_tbl, no_prefit_tbl):
    """Inner-join the two catalogues on SRCID, suffixing
    columns with _PRE / _NOP."""
    pre = prefit_tbl.copy()
    nop = no_prefit_tbl.copy()
    for c in pre.colnames:
        if c != "SRCID":
            pre.rename_column(c, c + "_PRE")
    for c in nop.colnames:
        if c != "SRCID":
            nop.rename_column(c, c + "_NOP")
    return join(pre, nop, keys="SRCID", join_type="inner")


def histograms(joint, param_pre, param_nop, label, out_dir,
               bin_label_arr=None, bins=40):
    """Side-by-side histogram of A vs B, plus per-bin overlays
    if bin_label_arr is provided."""
    a = np.asarray(joint[param_pre], dtype=float)
    b = np.asarray(joint[param_nop], dtype=float)
    valid = np.isfinite(a) & np.isfinite(b)
    a, b = a[valid], b[valid]

    # Use log scale for NH and Flux is implicit (lg10Flux already log)
    if label == "NH":
        a_plot = np.log10(np.clip(a, 1e-5, None))
        b_plot = np.log10(np.clip(b, 1e-5, None))
        xlabel = r"$\log_{10}(N_H / 10^{22}\ \mathrm{cm}^{-2})$"
    else:
        a_plot, b_plot = a, b
        xlabel = label

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5),
                             sharey=True)
    rng = (min(a_plot.min(), b_plot.min()),
           max(a_plot.max(), b_plot.max()))

    axes[0].hist(a_plot, bins=bins, range=rng, alpha=0.65,
                 label=f"prefit=True (N={len(a_plot)})",
                 color="C0")
    axes[0].hist(b_plot, bins=bins, range=rng, alpha=0.65,
                 label=f"no_prefit (N={len(b_plot)})",
                 color="C1")
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel("count")
    axes[0].set_title(f"{label} — pooled")
    axes[0].legend()
    if label == "NH":
        # Mark the LM starting value
        axes[0].axvline(np.log10(NH_STARTING_VALUE),
                        ls="--", color="k", alpha=0.6,
                        label="LM starting NH (5e20)")
        axes[0].legend()

    # Per-bin overlay if stratification bin labels are provided
    if bin_label_arr is not None:
        bin_label_arr = bin_label_arr[valid]
        for bn, color in zip(("low", "mid", "high"),
                             ("C2", "C3", "C4")):
            mask = bin_label_arr == bn
            if mask.sum() == 0:
                continue
            axes[1].hist(a_plot[mask], bins=bins, range=rng,
                         alpha=0.5, label=f"pre {bn}",
                         color=color, histtype="step",
                         linewidth=2, ls="-")
            axes[1].hist(b_plot[mask], bins=bins, range=rng,
                         alpha=0.5, label=f"nop {bn}",
                         color=color, histtype="step",
                         linewidth=2, ls=":")
        axes[1].set_xlabel(xlabel)
        axes[1].set_title(f"{label} — by stratification bin")
        axes[1].legend(fontsize=8, ncol=2)
        if label == "NH":
            axes[1].axvline(np.log10(NH_STARTING_VALUE),
                            ls="--", color="k", alpha=0.6)
    else:
        axes[1].set_visible(False)

    out = os.path.join(out_dir, f"hist_{label}.png")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  wrote {out}")


def scatter(joint, param_pre, param_nop, label, out_dir):
    """A vs B scatter with 1:1 line."""
    a = np.asarray(joint[param_pre], dtype=float)
    b = np.asarray(joint[param_nop], dtype=float)
    valid = np.isfinite(a) & np.isfinite(b)
    a, b = a[valid], b[valid]
    if label == "NH":
        a, b = np.log10(np.clip(a, 1e-5, None)), \
               np.log10(np.clip(b, 1e-5, None))
        unit = r"$\log_{10} N_H$"
    else:
        unit = label

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(a, b, s=4, alpha=0.4)
    lo = min(a.min(), b.min())
    hi = max(a.max(), b.max())
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.5)
    ax.set_xlabel(f"{unit} (prefit=True)")
    ax.set_ylabel(f"{unit} (no_prefit)")
    ax.set_title(f"{label}: prefit vs no_prefit")
    ax.set_aspect("equal")
    out = os.path.join(out_dir, f"scatter_{label}.png")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  wrote {out}")


def stats_summary(joint, param_pre, param_nop, label):
    """Compute KS, Mann-Whitney, median offset, and the
    starting-value-pile-up metric for NH."""
    a = np.asarray(joint[param_pre], dtype=float)
    b = np.asarray(joint[param_nop], dtype=float)
    valid = np.isfinite(a) & np.isfinite(b)
    a, b = a[valid], b[valid]
    if len(a) < 5:
        return {"label": label, "n": int(len(a)),
                "error": "too few points"}

    if label == "NH":
        a_log = np.log10(np.clip(a, 1e-5, None))
        b_log = np.log10(np.clip(b, 1e-5, None))
        # Counts within ±tol dex of LM starting value
        log_start = np.log10(NH_STARTING_VALUE)
        a_pile = np.sum(
            np.abs(a_log - log_start) < NH_STARTING_VALUE_TOL)
        b_pile = np.sum(
            np.abs(b_log - log_start) < NH_STARTING_VALUE_TOL)
        ks = ks_2samp(a_log, b_log)
        mw = mannwhitneyu(a_log, b_log, alternative="two-sided")
        return {
            "label": label,
            "n": int(len(a)),
            "median_prefit": float(np.median(a)),
            "median_no_prefit": float(np.median(b)),
            "median_diff_log": float(np.median(b_log)
                                     - np.median(a_log)),
            "ks_stat": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
            "mw_pvalue": float(mw.pvalue),
            "starting_value_pileup_prefit": int(a_pile),
            "starting_value_pileup_no_prefit": int(b_pile),
            "starting_value_window_dex":
                float(NH_STARTING_VALUE_TOL),
            "starting_value": float(NH_STARTING_VALUE),
        }
    else:
        ks = ks_2samp(a, b)
        mw = mannwhitneyu(a, b, alternative="two-sided")
        return {
            "label": label,
            "n": int(len(a)),
            "median_prefit": float(np.median(a)),
            "median_no_prefit": float(np.median(b)),
            "median_diff": float(np.median(b)
                                 - np.median(a)),
            "ks_stat": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
            "mw_pvalue": float(mw.pvalue),
        }


def main():
    parser = argparse.ArgumentParser(
        description="Compare prefit=True vs no_prefit catalogues")
    parser.add_argument("--prefit", required=True,
        help="FITS catalogue from production (prefit=True)")
    parser.add_argument("--no_prefit", required=True,
        help="FITS catalogue from validation (--no_prefit)")
    parser.add_argument("--metadata", default=None,
        help="Optional metadata FITS from "
             "select_validation_sample.py — used to overlay "
             "stratification-bin histograms")
    parser.add_argument("--output_dir", required=True,
        help="Directory for plots and JSON summary")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    pre = Table.read(args.prefit)
    nop = Table.read(args.no_prefit)
    print(f"prefit catalogue:    {len(pre)} rows")
    print(f"no_prefit catalogue: {len(nop)} rows")

    if "SRCID" not in pre.colnames or "SRCID" not in nop.colnames:
        print("ERROR: both catalogues must have SRCID column")
        sys.exit(1)

    # Resolve parameter column names per catalogue
    cols_pre = {p: find_col(pre, c, p)
                for p, c in PARAM_CANDIDATES.items()}
    cols_nop = {p: find_col(nop, c, p)
                for p, c in PARAM_CANDIDATES.items()}
    print("\nResolved columns:")
    for p in PARAM_CANDIDATES:
        print(f"  {p:>10}: pre={cols_pre[p]}  "
              f"nop={cols_nop[p]}")

    # Subset to just the needed columns (keeps join cheap)
    pre_sub = pre["SRCID",
                  cols_pre["NH"],
                  cols_pre["PHOINDEX"],
                  cols_pre["LG10FLUX"]]
    nop_sub = nop["SRCID",
                  cols_nop["NH"],
                  cols_nop["PHOINDEX"],
                  cols_nop["LG10FLUX"]]
    joint = joint_table(pre_sub, nop_sub)
    print(f"\nJoined on SRCID: {len(joint)} matching sources")

    # Optional bin labels from metadata
    bin_labels = None
    if args.metadata is not None and os.path.exists(args.metadata):
        meta = Table.read(args.metadata)
        if ("SRCID" in meta.colnames
                and "VALIDATION_BIN" in meta.colnames):
            joint_with_bin = join(
                joint, meta["SRCID", "VALIDATION_BIN"],
                keys="SRCID", join_type="left")
            bin_labels = np.asarray(
                joint_with_bin["VALIDATION_BIN"])
            joint = joint_with_bin
            print(f"  bin labels merged "
                  f"({np.unique(bin_labels)})")

    # Compute and write plots + stats
    summary = {}
    print("\nGenerating plots...")
    for p in PARAM_CANDIDATES:
        col_pre = cols_pre[p] + "_PRE"
        col_nop = cols_nop[p] + "_NOP"
        histograms(joint, col_pre, col_nop, p,
                   args.output_dir, bin_labels)
        scatter(joint, col_pre, col_nop, p, args.output_dir)
        s = stats_summary(joint, col_pre, col_nop, p)
        summary[p] = s

    # Write JSON summary
    summary_path = os.path.join(args.output_dir,
                                "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWrote {summary_path}")

    # Print headline numbers
    print("\n" + "=" * 60)
    print("HEADLINE NUMBERS (pasteable to email)")
    print("=" * 60)
    for p, s in summary.items():
        if "error" in s:
            print(f"\n{p}: {s['error']}")
            continue
        print(f"\n{p} (N={s['n']}):")
        print(f"  median(prefit)    = {s['median_prefit']:.4g}")
        print(f"  median(no_prefit) = {s['median_no_prefit']:.4g}")
        diff_key = ("median_diff_log" if "median_diff_log" in s
                    else "median_diff")
        print(f"  Δmedian           = {s[diff_key]:.4g}")
        print(f"  KS p-value        = {s['ks_pvalue']:.3g}")
        print(f"  Mann-Whitney p    = {s['mw_pvalue']:.3g}")
        if p == "NH":
            print(f"\n  Starting-value pile-up "
                  f"(±{s['starting_value_window_dex']} dex of "
                  f"5e20 cm^-2):")
            print(f"    prefit=True : "
                  f"{s['starting_value_pileup_prefit']} sources")
            print(f"    no_prefit   : "
                  f"{s['starting_value_pileup_no_prefit']} sources")
            ratio = (s['starting_value_pileup_prefit']
                     / max(1, s['starting_value_pileup_no_prefit']))
            print(f"    ratio       : {ratio:.2f}x")
            if ratio > 2.0:
                print("    >>> SIGNIFICANT pile-up in prefit "
                      "run: LM-non-convergence artifact "
                      "confirmed")
            else:
                print("    >>> No significant pile-up: prefit "
                      "is innocent for this peak")


if __name__ == "__main__":
    main()
