#!/usr/bin/env python3
"""
Recompute the KS goodness-of-fit for already-fitted 5XMM sources.

Why this exists
---------------
`_compute_goodness_of_fit()` in spectral_fitting_bxa_adapted.py feeds the KS
test with `Plot.y()` and `Plot.model()`.  With a background attached XSPEC
returns *background-subtracted* values from `Plot.y()`, which can be negative,
while `Plot.model()` is the source model alone.  Negative entries then reach
`_ks_2samp()`, where `cumsum()/sum()` stops being a monotonic CDF.

The original (Sherpa) implementation this was ported from did it correctly:
`utils.get_data_model_counts()` returns the raw observed counts together with
the *full* model including the fitted background component.  This script
restores those semantics for the XSPEC path:

    Plot.background = True
    obs = Plot.y(i)     + Plot.backgroundVals(i)     # total observed counts
    mod = Plot.model(i) + Plot.backgroundVals(i)     # total predicted counts

Only the goodness-of-fit is affected: the BXA fit itself used cstat with an
attached background (the W-statistic), which never subtracts and is correct.

No re-fitting is needed.  The spectra are reloaded and the model parameters are
set to the stored posterior medians, so this is ~1 s per source.

Usage
-----
    python3 recompute_gof.py --results fit_results_all.fits \\
        --output_dir /storage/.../pipeline_output \\
        --data_dir /storage/xmmcat/5XMM_data/Spectra --subdir pps \\
        --out gof_chunk_000.fits --chunk 0 --nchunks 100

Outputs one FITS table per chunk with the new and old GoF side by side.
"""

import argparse
import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits
from astropy.table import Table

# Thresholds — must match spectral_fitting_bxa_adapted.py
GOF_PVALUE_THRESHOLD = 0.01
PHOINDEX_PEG_MARGIN = 0.05
NH_PEG_THRESHOLD = 100.0          # 1e24 cm^-2, Compton-thick
FIT_BAND = "**-0.3 10.0-**"
CFLUX_BAND = (0.5, 10.0)
NITER = 1000


# --------------------------------------------------------------------------
# statistics (same definitions as the pipeline / upstream stats.py)
# --------------------------------------------------------------------------
def ks_2samp(a, b):
    """KS statistic between two count arrays (CDF comparison)."""
    return np.abs(a.cumsum() / a.sum() - b.cumsum() / b.sum()).max()


def ks_permutation(data, model, niter=NITER, seed=None):
    """Fixed-label permutation test: p = fraction of permutations with KS >= observed.

    Seeded so results are reproducible (the pipeline version was unseeded).
    """
    ks_obs = ks_2samp(data, model)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(niter):
        mask = rng.choice([False, True], size=len(data))
        if ks_2samp(np.where(mask, data, model),
                    np.where(mask, model, data)) >= ks_obs:
            count += 1
    return float(ks_obs), count / niter


def quality_flag(ks_p, phoindex, nh):
    """flag = 4 + (gof_bad + 2*pho_pegged + 4*nh_pegged), 0 if no issue."""
    gof_bad = bool(np.isfinite(ks_p) and ks_p < GOF_PVALUE_THRESHOLD)
    pho_peg = bool(phoindex is not None and np.isfinite(phoindex) and
                   (phoindex <= 1.0 + PHOINDEX_PEG_MARGIN or
                    phoindex >= 3.0 - PHOINDEX_PEG_MARGIN))
    nh_peg = bool(nh is not None and np.isfinite(nh) and nh >= NH_PEG_THRESHOLD)
    combined = int(gof_bad) + 2 * int(pho_peg) + 4 * int(nh_peg)
    return (4 + combined) if combined > 0 else 0


# --------------------------------------------------------------------------
# per-source spectral setup
# --------------------------------------------------------------------------
def _obsid_from_pps_name(name):
    """'P0672050201PNS003BGSPEC0002.FTZ' -> '0672050201'."""
    base = os.path.basename(name)
    return base[1:11] if base.startswith("P") and len(base) > 11 else None


def resolve_files(srcid, inst, src_dir, data_dir, subdir):
    """Return (grp, rmf, arf, bkg) absolute paths for one instrument.

    The .grp and .pha are real files and the RMF symlink in the source
    directory is still valid, so the RMF is resolved through it (this handles
    the PN '_v22.0' suffix and the MOS '<N>eV/' subdirectory automatically).
    The BGSPEC/SRCARF symlinks point at the old '<obsid>/product/' layout and
    are dead, so those are rebuilt from the header names against `subdir`.
    """
    grp = os.path.join(src_dir, f"{srcid}_{inst}.grp")
    if not os.path.exists(grp):
        raise FileNotFoundError(f"missing {grp}")

    with fits.open(grp) as h:
        hdr = h[1].header
        backfile = hdr.get("BACKFILE")
        ancrfile = hdr.get("ANCRFILE")
        respfile = hdr.get("RESPFILE")

    if not (backfile and ancrfile and respfile):
        raise ValueError(f"{grp}: missing BACKFILE/ANCRFILE/RESPFILE keyword")

    # RMF via the (valid) symlink next to the spectrum
    rmf = os.path.realpath(os.path.join(src_dir, os.path.basename(respfile)))
    if not os.path.exists(rmf):
        raise FileNotFoundError(f"missing RMF {rmf}")

    obsid = _obsid_from_pps_name(backfile)
    if obsid is None:
        raise ValueError(f"cannot parse OBS_ID from {backfile}")
    prod = os.path.join(data_dir, obsid, subdir)
    arf = os.path.join(prod, os.path.basename(ancrfile))
    bkg = os.path.join(prod, os.path.basename(backfile))
    for f in (arf, bkg):
        if not os.path.exists(f):
            raise FileNotFoundError(f"missing {f}")
    return grp, rmf, arf, bkg


def load_source(srcid, instruments, src_dir, data_dir, subdir, scratch):
    """Load all instruments of one source into XSPEC data groups.

    Returns the list of instruments actually loaded, in group order.
    """
    from xspec import AllData, AllModels

    AllData.clear()
    AllModels.clear()

    specs = []
    for inst in instruments:
        grp, rmf, arf, bkg = resolve_files(srcid, inst, src_dir,
                                           data_dir, subdir)
        # XSPEC resolves BACKFILE/ANCRFILE/RESPFILE at load time and would hit
        # the dead symlinks, so work on a copy with those keywords cleared.
        local = os.path.join(scratch, os.path.basename(grp))
        shutil.copy2(grp, local)
        with fits.open(local, mode="update") as h:
            for k in ("BACKFILE", "ANCRFILE", "RESPFILE"):
                if k in h[1].header:
                    h[1].header[k] = "none"
        specs.append((inst, local, rmf, arf, bkg))

    expr = " ".join(f"{i+1}:{i+1} {s[1]}" for i, s in enumerate(specs))
    AllData(expr)
    for i, (inst, local, rmf, arf, bkg) in enumerate(specs, start=1):
        AllData(i).response = rmf
        AllData(i).response.arf = arf
        AllData(i).background = bkg
    AllData.ignore(FIT_BAND)
    return [s[0] for s in specs]


def build_model(loaded, nh, phoindex, lg10flux, factor):
    """Rebuild the fitted model and set parameters to the posterior medians."""
    from xspec import AllModels, Fit, Model

    Fit.statMethod = "cstat"
    use_iin = len(loaded) > 1
    expr = ("constant*" if use_iin else "") + "cflux*phabs*zpowerlw"
    m = Model(expr)
    m.zpowerlw.Redshift = 0.0
    m.zpowerlw.Redshift.frozen = True
    # cflux carries the normalisation, so zpowerlw.norm is frozen (as the
    # pipeline does).  Without this the dof is one too low.
    m.zpowerlw.norm.frozen = True
    m.cflux.Emin, m.cflux.Emax = CFLUX_BAND
    m.phabs.nH.values = [float(nh)]
    m.zpowerlw.PhoIndex.values = [float(phoindex)]
    m.cflux.lg10Flux.values = [float(lg10flux)]

    if use_iin:
        # pn is the reference (constant frozen at 1); others carry the IIN
        try:
            ref = loaded.index("pn") + 1
        except ValueError:
            ref = 1
        m.constant.factor.values = "1.0"
        m.constant.factor.frozen = True
        for gi in range(2, len(loaded) + 1):
            g = AllModels(gi)
            g.constant.factor.untie()
            if gi == ref:
                g.constant.factor.values = "1.0"
                g.constant.factor.frozen = True
            else:
                g.constant.factor.frozen = False
                g.constant.factor.values = [float(factor)]
        if ref != 1:
            m.constant.factor.frozen = False
            m.constant.factor.values = [float(factor)]
    return m


def gof_arrays(nspec, corrected):
    """Concatenated (data, model) arrays.

    corrected=True  -> raw total counts on both sides (background re-added)
    corrected=False -> the original, background-subtracted behaviour
    """
    from xspec import Plot

    Plot.device = "/null"
    Plot.xAxis = "channel"
    Plot.background = corrected
    Plot("counts")

    data, model = [], []
    for si in range(1, nspec + 1):
        y = np.array(Plot.y(si))
        mo = np.array(Plot.model(si))
        if corrected:
            bg = np.array(Plot.backgroundVals(si))
            y = y + bg
            mo = mo + bg
        data.append(y)
        model.append(mo)
    d = np.concatenate(data)
    mo = np.concatenate(model)

    mask = mo > 0
    if corrected:
        mask &= d > 0          # also drop non-positive DATA (the original bug)
    return d[mask], mo[mask]


def process_one(row, output_dir, data_dir, subdir, scratch, also_old=False):
    """Recompute the GoF for one source. Returns a result dict."""
    from xspec import AllData, AllModels, Fit

    srcid = int(row["SRCID"])
    out = {
        "SRCID": srcid, "ks_stat_new": np.nan, "ks_pvalue_new": np.nan,
        "ks_stat_old_repro": np.nan, "ks_pvalue_old_repro": np.nan,
        "cstat_new": np.nan, "dof_new": -1, "flag_new": -1,
        "n_bins": -1, "n_negative_old": -1, "status": "ok",
    }

    inst_str = str(row["instruments"]) if row["instruments"] else ""
    instruments = [i.strip() for i in inst_str.split(",") if i.strip()]
    if not instruments:
        instruments = ["pn", "MOS"]
    # normalise MOS1/MOS2 -> the merged 'MOS' spectrum written by merge_spectra
    instruments = ["MOS" if i.upper().startswith("M") else "pn"
                   for i in instruments]
    seen, ordered = set(), []
    for i in instruments:
        if i not in seen:
            seen.add(i)
            ordered.append(i)

    src_dir = os.path.join(output_dir, str(srcid))
    try:
        loaded = load_source(srcid, ordered, src_dir, data_dir,
                             subdir, scratch)
        build_model(loaded, row["nH_median"], row["PhoIndex_median"],
                    row["lg10Flux_median"], row.get("factor_median", 1.0))

        out["cstat_new"] = float(Fit.statistic)
        out["dof_new"] = int(Fit.dof)

        d, mo = gof_arrays(len(loaded), corrected=True)
        if d.size == 0:
            out["status"] = "no_bins"
            return out
        k, p = ks_permutation(d, mo, seed=srcid)   # seeded on SRCID
        out["ks_stat_new"], out["ks_pvalue_new"] = k, p
        out["n_bins"] = int(d.size)
        out["flag_new"] = quality_flag(p, row["PhoIndex_median"],
                                       row["nH_median"])

        if also_old:
            d0, mo0 = gof_arrays(len(loaded), corrected=False)
            out["n_negative_old"] = int((d0 < 0).sum())
            k0, p0 = ks_permutation(d0, mo0, seed=srcid)
            out["ks_stat_old_repro"], out["ks_pvalue_old_repro"] = k0, p0
    except Exception as exc:                       # noqa: BLE001
        out["status"] = f"error: {type(exc).__name__}: {exc}"[:120]
    finally:
        try:
            AllData.clear()
            AllModels.clear()
        except Exception:                          # noqa: BLE001
            pass
        for f in os.listdir(scratch):
            try:
                os.remove(os.path.join(scratch, f))
            except OSError:
                pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True,
                    help="fit_results_all.fits from the production run")
    ap.add_argument("--output_dir", required=True,
                    help="pipeline_output directory (per-source subdirs)")
    ap.add_argument("--data_dir", required=True,
                    help="Spectra tree root")
    ap.add_argument("--subdir", default="pps",
                    help="subdirectory under each OBS_ID (default: pps)")
    ap.add_argument("--out", required=True, help="output FITS table")
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0,
                    help="process at most N sources (for testing)")
    ap.add_argument("--also_old", action="store_true",
                    help="also reproduce the old (subtracted) GoF, for comparison")
    ap.add_argument("--flush_every", type=int, default=200)
    args = ap.parse_args()

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    t = Table.read(args.results)
    # only rows with an actual fit
    fitted = np.isin(np.array(t["flag"]), [0, 5, 6, 7, 8, 9, 10, 11])
    t = t[fitted]
    if args.nchunks > 1:
        t = t[args.chunk::args.nchunks]
    if args.limit:
        t = t[:args.limit]
    print(f"[chunk {args.chunk}] {len(t)} sources to process", flush=True)

    scratch = tempfile.mkdtemp(prefix=f"gof_{args.chunk}_")
    rows = []
    try:
        for n, row in enumerate(t, 1):
            rows.append(process_one(dict(zip(row.colnames, row)),
                                    args.output_dir, args.data_dir,
                                    args.subdir, scratch,
                                    also_old=args.also_old))
            if n % args.flush_every == 0:
                Table(rows=rows).write(args.out, overwrite=True)
                print(f"[chunk {args.chunk}] {n}/{len(t)}", flush=True)
    finally:
        if rows:
            Table(rows=rows).write(args.out, overwrite=True)
        shutil.rmtree(scratch, ignore_errors=True)

    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"[chunk {args.chunk}] done: {ok}/{len(rows)} ok -> {args.out}",
          flush=True)


if __name__ == "__main__":
    main()
