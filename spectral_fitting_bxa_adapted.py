import os
import gc
import signal
import datetime

# Prevent numpy/MKL/OpenMP from spawning threads that
# conflict with process-level parallelism (xargs -P / MPI)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

# Use offscreen Qt backend to prevent ICE/PyQt5 segfault
# at exit (no X11 session manager needed for headless fits)
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from xspec import *
Fit.query = "no"  # Disable interactive prompts during fits
import bxa.xspec as bxa
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.table import Table, vstack
import glob
import logging

# python3 automated_fits.py 3067718060100029 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits dummy_output.txt --use_bxa --model_name=powerlaw --redshift=1.0 --overwrite=1 --export_results_fits --export_filename=fit_results.fits --bxa_output_dir=bxa_fit_results

#logger = logging.getLogger(__name__)
logger = logging.getLogger()   # root logger


class BXATimeout(Exception):
    """Raised when a BXA fit exceeds its wallclock time limit."""
    pass


def _bxa_timeout_handler(signum, frame):
    raise BXATimeout("BXA fit exceeded time limit")


def _get_open_fds():
    """Return set of open fd numbers (Linux /proc/self/fd).
    Used as a safety net to detect and close leaked file
    descriptors after each BXA fit, regardless of source."""
    try:
        return set(int(fd) for fd in os.listdir('/proc/self/fd'))
    except OSError:
        return set()


# Extensions / prefixes of fds that must NOT be closed by the
# fd fence — they are long-lived library caches, not leaks.
_FD_KEEP_SUFFIXES = (
    '.ttf', '.otf', '.woff',   # matplotlib font cache
    '.so', '.dylib',            # shared libraries
    '.pyc', '.pyo',             # compiled Python
)
_FD_KEEP_PREFIXES = (
    '/dev/', '/proc/', '/sys/', # OS virtual filesystems
    'pipe:', 'socket:', 'anon_inode:',  # non-file fds
)


def _is_safe_to_close(fd):
    """Check if a leaked fd can be safely closed.
    Returns True for fit-related files (HDF5, logs, FITS),
    False for library caches (fonts, .so) and OS resources."""
    try:
        target = os.readlink(f'/proc/self/fd/{fd}')
    except OSError:
        return False  # can't read → don't touch
    for pfx in _FD_KEEP_PREFIXES:
        if target.startswith(pfx):
            return False
    for sfx in _FD_KEEP_SUFFIXES:
        if target.endswith(sfx):
            return False
    return True


def _create_logflat_prior_for(model, param):
    """Log-uniform (Jeffreys) prior for an XSPEC parameter.

    Samples uniformly in log-space — symmetric around 1.0 in
    ratio — but stores the LINEAR value in the chain.  This
    avoids bxa.create_loguniform_prior_for's post-run bug where
    BXA tries to set the log10 value back to the XSPEC parameter
    without exponentiating.

    Works by creating a standard BXA uniform prior (to discover
    BXA's expected return format), then replacing the transform
    with one that maps [0,1] → 10^(loglo + u*(loghi − loglo)).
    The chain column keeps the original XSPEC name ("factor",
    not "log(factor)").

    Handles both BXA prior formats:
      - newer BXA: returns a dict  {'paramname':…, 'transform':…, …}
      - older BXA: returns a callable with .paramname attribute
    """
    lo = float(param.values[2])   # XSPEC hard min
    hi = float(param.values[5])   # XSPEC hard max
    log_lo = np.log10(lo)
    log_hi = np.log10(hi)

    def logflat_transform(cube):
        return 10 ** (log_lo + cube * (log_hi - log_lo))

    # Build a standard BXA prior to discover expected return format
    base = bxa.create_uniform_prior_for(model, param)

    if isinstance(base, dict):
        # Newer BXA: priors are dicts — replace the transform in-place
        base['transform'] = logflat_transform
        return base
    else:
        # Classic BXA: priors are callables with .paramname attribute
        logflat_transform.paramname = base.paramname
        if hasattr(base, 'param'):
            logflat_transform.param = base.param
        return logflat_transform


def get_model_and_priors(model_name, redshift=0.0,
                         flux_band=(0.5, 10.0),
                         prefit=True, n_groups=1,
                         instruments=None,
                         use_iin=True,
                         iin_bounds=(0.5, 2.0),
                         observed_flux=True):
    """
    Construct an XSPEC model wrapped with cflux, optionally
    run a quick pre-fit to tighten priors, then build BXA
    priors.

    Parameters
    ----------
    model_name : str
        Name of the physical model
        ("powerlaw", "apec_single", "blackbody", "bremss").
    redshift : float
        Redshift for models that need it.
        Defaults to 0.0 (no redshift correction).
    flux_band : tuple
        Energy band (Emin, Emax) in keV for cflux.
    prefit : bool
        If True, run a quick Fit.perform() and tighten the
        parameter ranges around the best-fit before creating
        priors. This dramatically speeds up ultranest by
        reducing the prior volume.
    n_groups : int
        Number of XSPEC data groups (one per spectrum).
        Used to configure per-group IIN constants.
    instruments : list of str or None
        Instrument identifier per data group (e.g.
        ["pn", "MOS", "MOS"]). Used to determine which
        group is the reference (pn → frozen constant=1).
        If None, group 1 is used as reference.
    use_iin : bool
        If True, prepend a multiplicative constant to the
        model for inter-instrument normalisation (IIN).
        The constant is frozen at 1.0 for the reference
        instrument (pn) and free for others (MOS).
    iin_bounds : tuple
        (lo, hi) bounds for the free IIN constant
        (default 0.5–2.0). The prior is log-uniform
        (flat in log-space), so it is symmetric around
        1.0 in ratio (per Carrera, following XMM2ATHENA).
    observed_flux : bool
        If True (default), use cflux*phabs*model so that
        lg10Flux measures the observed (absorbed) flux
        (XMM2ATHENA convention, Viitanen+25).
        If False, use phabs*cflux*model so that lg10Flux
        measures the intrinsic (unabsorbed) flux
        (5XMM legacy convention).
    """
    # Defensive: ensure redshift is a valid float
    if redshift is None:
        redshift = 0.0

    # Build the model expression string
    # Absorption + emission core
    if observed_flux:
        # cflux wraps phabs*emission → observed flux
        flux_prefix = "cflux*phabs*"
    else:
        # phabs wraps cflux*emission → intrinsic flux
        flux_prefix = "phabs*cflux*"

    # Prepend IIN constant if requested
    iin_prefix = "constant*" if use_iin else ""

    if model_name == "powerlaw":
        expr = f"{iin_prefix}{flux_prefix}zpowerlw"
        model = Model(expr)
        model.zpowerlw.Redshift = redshift
        model.zpowerlw.Redshift.frozen = True
        model.phabs.nH.values = \
            "0.05,,0.001,0.001,1000.0,1000.0"
        model.zpowerlw.PhoIndex.values = \
            "2.0,,1.0,1.0,3.0,3.0"
        model.zpowerlw.norm.frozen = True

    elif model_name == "apec_single":
        expr = f"{iin_prefix}{flux_prefix}apec"
        model = Model(expr)
        model.phabs.nH.values = \
            "0.05,,0.001,0.001,1000.0,1000.0"
        model.apec.kT.values = \
            "1.0,,0.1,0.1,10.0,10.0"
        model.apec.norm.frozen = True

    elif model_name == "blackbody":
        expr = f"{iin_prefix}{flux_prefix}bbody"
        model = Model(expr)
        model.phabs.nH.values = \
            "0.05,,0.001,0.001,1000.0,1000.0"
        model.bbody.kT.values = \
            "0.1,,0.01,0.01,2.0,2.0"
        model.bbody.norm.frozen = True

    elif model_name == "bremss":
        expr = f"{iin_prefix}{flux_prefix}bremss"
        model = Model(expr)
        model.phabs.nH.values = \
            "0.05,,0.001,0.001,1000.0,1000.0"
        model.bremss.kT.values = \
            "5.0,,0.1,0.1,20.0,20.0"
        model.bremss.norm.frozen = True

    else:
        raise ValueError(f"Unknown model: {model_name}")

    logger.info(f"   Model expression: {expr}")

    # Configure cflux energy range
    model.cflux.Emin = flux_band[0]
    model.cflux.Emax = flux_band[1]

    # Starting value for lg10Flux
    model.cflux.lg10Flux.values = \
        "-12.0,,-15.0,-15.0,-9.0,-9.0"

    # -------------------------------------------------
    # Configure IIN constants per data group
    # -------------------------------------------------
    if use_iin:
        # Determine which group is the reference (pn)
        ref_group = 1  # default: first group is reference
        if instruments is not None:
            for gi, inst in enumerate(instruments):
                if inst.lower() == "pn":
                    ref_group = gi + 1  # 1-indexed
                    break

        # Group 1 (master model): freeze constant
        model.constant.factor.values = "1.0"
        model.constant.factor.frozen = True

        # Groups 2..N: untie and configure constant
        iin_lo, iin_hi = iin_bounds
        for gi in range(2, n_groups + 1):
            gmod = AllModels(gi)
            gmod.constant.factor.untie()
            if gi == ref_group:
                # This group is pn (reference) — freeze
                gmod.constant.factor.values = "1.0"
                gmod.constant.factor.frozen = True
            else:
                # MOS group — free constant
                gmod.constant.factor.frozen = False
                gmod.constant.factor.values = (
                    f"1.0,,{iin_lo},{iin_lo},"
                    f"{iin_hi},{iin_hi}")

        # If group 1 is NOT the reference (pn is in a later
        # group), swap: free group 1, freeze group ref_group
        if ref_group != 1:
            # Free group 1 (it's MOS)
            model.constant.factor.frozen = False
            model.constant.factor.values = (
                f"1.0,,{iin_lo},{iin_lo},"
                f"{iin_hi},{iin_hi}")
            # Freeze the pn group (already done above)

        logger.info(
            f"   IIN: {n_groups} groups, reference=group "
            f"{ref_group} "
            f"({'pn' if instruments else 'first'}), "
            f"free range [{iin_lo}, {iin_hi}]")

    # -------------------------------------------------
    # Quick pre-fit to tighten parameter ranges
    # -------------------------------------------------
    if prefit:
        _prefit_and_tighten(model, model_name, logger)

    # -------------------------------------------------
    # Build BXA priors from (possibly tightened) ranges
    # -------------------------------------------------
    # nH: log-uniform (Jeffreys) prior — flat in log-space.
    # Matches the D6.2/Viitanen+25 (aa55271-25) convention of
    # a flat prior on log10(nH/cm^-2). Uses our custom
    # _create_logflat_prior_for to store the LINEAR value in
    # the chain, avoiding BXA's post-run conversion bug for
    # XSPEC parameters that don't natively live in log space.
    # The previous linear-uniform prior is the upstream-port
    # regression that produced the 5e20 cm^-2 starting-value
    # pile-up in 5XMM nH histograms.
    priors = [
        _create_logflat_prior_for(
            model, model.phabs.nH),
    ]

    if model_name == "powerlaw":
        priors.append(bxa.create_uniform_prior_for(
            model, model.zpowerlw.PhoIndex))
    elif model_name == "apec_single":
        priors.append(bxa.create_uniform_prior_for(
            model, model.apec.kT))
    elif model_name == "blackbody":
        priors.append(bxa.create_uniform_prior_for(
            model, model.bbody.kT))
    elif model_name == "bremss":
        priors.append(bxa.create_uniform_prior_for(
            model, model.bremss.kT))

    # lg10Flux: uniform in log-space = log-uniform on flux
    priors.append(bxa.create_uniform_prior_for(
        model, model.cflux.lg10Flux))

    # IIN constant priors (one per free constant).
    # Log-uniform (Jeffreys) prior — flat in log-space, so
    # the prior is symmetric around 1.0 in ratio: P(factor=2)
    # = P(factor=0.5).  Per Carrera / XMM2ATHENA convention.
    # Uses our custom _create_logflat_prior_for (not BXA's
    # create_loguniform_prior_for) to store the LINEAR value
    # in the chain and avoid BXA's post-run conversion bug.
    if use_iin:
        for gi in range(1, n_groups + 1):
            gmod = AllModels(gi)
            if not gmod.constant.factor.frozen:
                priors.append(
                    _create_logflat_prior_for(
                        gmod, gmod.constant.factor))

    return model, priors


def _prefit_and_tighten(model, model_name, logger):
    """
    Run a quick XSPEC fit (Levenberg-Marquardt) and tighten
    the hard parameter limits around the best-fit values for
    the well-constrained parameters (lg10Flux, shape).

    This reduces the prior volume for those parameters,
    making ultranest converge faster.

    The tightening strategy:
    - nH:       NEVER tightened (log-uniform prior over 6 dex
                is sampled efficiently by ultranest; tightening
                it creates the 5e20 cm^-2 pile-up artefact
                described in D6.2 and Viitanen+25).
    - lg10Flux: best-fit ± 2 dex, clamped to [-15, -9].
                Tighten only if LM actually moved it.
                (Flux is well-determined by total counts, so
                the tightening is artefact-safe and provides
                a modest BXA speed-up.)
    - PhoIndex/kT: NEVER tightened. Same rationale as nH —
                   the original priors are already narrow
                   enough for ultranest, and tightening
                   creates a starting-value pile-up at the
                   LM result for low-S/N sources.

    Note: the IIN constant (if present) is NOT tightened here.
    It is configured separately in get_model_and_priors() and
    always keeps its physical range [iin_lo, iin_hi]. LM will
    fit it during Fit.perform() but the result is not used to
    modify the constant's prior bounds.

    Per-parameter convergence check:
    For lg10Flux and the shape parameter, the prior is
    tightened ONLY if Levenberg-Marquardt actually moved it
    from its starting value (relative change >
    PARAM_MOVE_THRESHOLD or, for log-space lg10Flux, absolute
    change > FLUX_MOVE_THRESHOLD). If a parameter did not
    move, its prior is left at the original wide range.

    If the pre-fit raises an exception, ALL priors are left
    at the original wide ranges (safe fallback).
    """
    # Thresholds for "did LM actually move this parameter?"
    PARAM_MOVE_THRESHOLD = 0.01   # 1% relative change
    FLUX_MOVE_THRESHOLD = 0.05    # 0.05 dex in log10(flux)

    # ----- Record starting values BEFORE Fit.perform() -----
    nh_start = float(model.phabs.nH.values[0])
    flux_start = float(model.cflux.lg10Flux.values[0])

    if model_name == "powerlaw":
        shape_start = float(model.zpowerlw.PhoIndex.values[0])
        shape_name = "PhoIndex"
    elif model_name == "apec_single":
        shape_start = float(model.apec.kT.values[0])
        shape_name = "kT(apec)"
    elif model_name == "blackbody":
        shape_start = float(model.bbody.kT.values[0])
        shape_name = "kT(bbody)"
    elif model_name == "bremss":
        shape_start = float(model.bremss.kT.values[0])
        shape_name = "kT(bremss)"
    else:
        shape_start = None
        shape_name = None

    try:
        Fit.perform()
        logger.info("   Pre-fit completed: "
                     f"chi2/dof = {Fit.statistic:.1f}"
                     f"/{Fit.dof}")

        # ----- Read post-fit values -----
        nh_val = float(model.phabs.nH.values[0])
        flux_val = float(model.cflux.lg10Flux.values[0])

        # Per-parameter convergence checks
        nh_moved = (abs(nh_val - nh_start)
                    / max(abs(nh_start), 1e-6)
                    > PARAM_MOVE_THRESHOLD)
        flux_moved = (abs(flux_val - flux_start)
                      > FLUX_MOVE_THRESHOLD)

        # Log overall convergence status
        logger.info(
            f"   Pre-fit movement: "
            f"nH {nh_start:.4g}→{nh_val:.4g} "
            f"({'moved' if nh_moved else 'STUCK'}), "
            f"lg10Flux {flux_start:.2f}→{flux_val:.2f} "
            f"({'moved' if flux_moved else 'STUCK'})")

        # --- nH: NEVER tighten ---
        # With a log-uniform prior over 6 dex [10^19, 10^25],
        # ultranest samples efficiently without prefit help.
        # Tightening Nh is precisely what creates the
        # starting-value pile-up at 5e20 (the "attractor"
        # effect): for low-S/N sources, LM lands near 0.05,
        # the prior collapses around it, and the BXA posterior
        # median follows. We log the LM movement for diagnostics
        # only.
        logger.info(
            f"   nH NOT tightened (log-uniform prior kept "
            f"wide); LM {'moved' if nh_moved else 'STUCK'} "
            f"at {nh_val:.4g}")

        # --- lg10Flux: tighten only if LM moved it ---
        if flux_moved:
            flux_lo = max(-15.0, flux_val - 2.0)
            flux_hi = min(-9.0, flux_val + 2.0)
            model.cflux.lg10Flux.values = (
                f"{flux_val},,{flux_lo},{flux_lo},"
                f"{flux_hi},{flux_hi}")
            logger.info(f"   Tightened lg10Flux: range "
                        f"[{flux_lo:.2f}, {flux_hi:.2f}]")
        else:
            logger.info(
                "   lg10Flux NOT tightened (LM did not move "
                "it); keeping original wide prior [-15, -9]")

        # --- Shape parameter: NEVER tighten ---
        # Same rationale as for nH: tightening the prior around
        # the LM result creates an attractor at the starting
        # value for low-S/N sources (PhoIndex pile-up at 2.0,
        # kT pile-up at the kT start). The original priors are
        # already narrow enough for ultranest to sample
        # efficiently. Diagnostic-only logging.
        if model_name == "powerlaw":
            ph_val = float(
                model.zpowerlw.PhoIndex.values[0])
            ph_moved = (abs(ph_val - shape_start)
                        / max(abs(shape_start), 1e-6)
                        > PARAM_MOVE_THRESHOLD)
            logger.info(
                f"   PhoIndex NOT tightened (kept wide); "
                f"LM {'moved' if ph_moved else 'STUCK'} "
                f"at {ph_val:.3f}")

        elif model_name == "apec_single":
            kt_val = float(model.apec.kT.values[0])
            kt_moved = (abs(kt_val - shape_start)
                        / max(abs(shape_start), 1e-6)
                        > PARAM_MOVE_THRESHOLD)
            logger.info(
                f"   kT(apec) NOT tightened (kept wide); "
                f"LM {'moved' if kt_moved else 'STUCK'} "
                f"at {kt_val:.3f}")

        elif model_name == "blackbody":
            kt_val = float(model.bbody.kT.values[0])
            kt_moved = (abs(kt_val - shape_start)
                        / max(abs(shape_start), 1e-6)
                        > PARAM_MOVE_THRESHOLD)
            logger.info(
                f"   kT(bbody) NOT tightened (kept wide); "
                f"LM {'moved' if kt_moved else 'STUCK'} "
                f"at {kt_val:.3f}")

        elif model_name == "bremss":
            kt_val = float(model.bremss.kT.values[0])
            kt_moved = (abs(kt_val - shape_start)
                        / max(abs(shape_start), 1e-6)
                        > PARAM_MOVE_THRESHOLD)
            logger.info(
                f"   kT(bremss) NOT tightened (kept wide); "
                f"LM {'moved' if kt_moved else 'STUCK'} "
                f"at {kt_val:.3f}")

    except Exception as e:
        logger.warning(
            f"   Pre-fit failed ({e}), using wide priors")
    
    
# === NEW: unified background check that always returns flag=3 on any problem ===
def check_background_fit(spectrum_file, background_file, rmf_file, arf_file,
                         model_name, redshift, logger, srcid="unknown"):
    """
    Returns:
      pval (float) if background is OK (>= 0.01),
      None if background cannot be fit or p<0.01 or any exception occurred.
    """
    from xspec import AllData, AllModels, Fit, Spectrum
    import numpy as np
    # Defensive: ensure redshift is a valid float
    if redshift is None:
        redshift = 0.0
    try:
        # Clean local XSPEC state for the quick test
        AllData.clear()
        AllModels.clear()

        # Change to the spectrum directory so XSPEC can resolve
        # filenames referenced in the FITS header
        spec_dir = os.path.dirname(spectrum_file)
        if spec_dir:
            os.chdir(spec_dir)

        s = Spectrum(spectrum_file)
        s.background = background_file
        s.response = rmf_file
        s.response.arf = arf_file

        AllData.ignore("**-0.3 10.0-**")


        # Build the same model you'll use (so the set-up is consistent),
        # but *do not* run BXA here—just a quick XSPEC fit to get p-value.
        # prefit=False: we do our own Fit.perform() below
        try:
            _model, _priors = get_model_and_priors(model_name, redshift, prefit=False)
        except Exception as e:
            logger.warning(f"   Background check: could not build model for {srcid}: {e}")
            return None

        # Perform a quick local fit
        Fit.perform()

        # Convert the fit test statistic to a p-value
        try:
            from scipy.stats import chi2 as chi2_dist
            bg_stat = Fit.testStatistic   # property, not method
            dof = Fit.dof
            pval = 1.0 - chi2_dist.cdf(bg_stat, dof)
        except Exception as e:
            logger.warning(f"   Background check: cannot compute p-value for {srcid}: {e}")
            return None

        if (pval is None) or (not np.isfinite(pval)) or (pval < 0.01):
            logger.warning(f"   Background check: p={pval:.4f} (<0.01) → NOT OK for {srcid}")
            return None

        logger.info(f"   Background check: p={pval:.4f} → OK for {srcid}")
        return pval

    except Exception as e:
        logger.warning(f"   Background check failed for {srcid}: {e}")
        return None
    finally:
        # Always clean up XSPEC state to avoid fd leaks
        AllData.clear()
        AllModels.clear()



def _compute_goodness_of_fit(labels, posterior_median,
                             samples_array, logger):
    """Compute goodness-of-fit statistics while XSPEC is loaded.

    Called after BXA has finished but BEFORE AllData.clear().
    Sets model params to posterior medians, then computes:
      - C-stat and dof from XSPEC
      - KS statistic and bootstrap p-value (data vs model counts)

    Returns a dict with keys: cstat, dof, ks_stat, ks_pvalue.
    """
    gof = {"cstat": np.nan, "dof": -1,
           "ks_stat": np.nan, "ks_pvalue": np.nan}

    try:
        # Set model parameters to posterior medians.
        # Must handle two BXA naming conventions:
        #   - uniform priors:     label = XSPEC name (e.g. "nH")
        #   - log-uniform priors: label = "log(name)" and the
        #     chain stores log10(value), so convert back.
        # Search across ALL data groups (not just group 1)
        # so that IIN constants on groups 2+ get set too.
        n_groups = AllData.nSpectra
        for name, value in zip(labels, posterior_median):
            if name.startswith("log(") and name.endswith(")"):
                xspec_name = name[4:-1]   # "log(factor)" → "factor"
                xspec_value = 10**value   # log10 → linear
            else:
                xspec_name = name
                xspec_value = float(value)

            set_ok = False
            for gi in range(1, n_groups + 1):
                m = AllModels(gi)
                for ci in range(1, m.nParameters + 1):
                    p = m(ci)
                    if p.name == xspec_name and not p.frozen:
                        p.values = [float(xspec_value)]
                        set_ok = True
                        break
                if set_ok:
                    break

        # C-stat and dof — just read, do NOT re-fit.
        # Fit.perform() would re-optimise from here and
        # can step outside tightened hard limits.  After
        # setting parameters, Plot("counts") forces XSPEC
        # to evaluate the folded model, which updates the
        # statistic.
        Plot.device = "/null"
        Plot.xAxis = "channel"
        Plot("counts")
        gof["cstat"] = float(Fit.statistic)
        gof["dof"] = int(Fit.dof)

        # Get observed and model-predicted counts per channel
        # across all loaded spectra (for KS test).
        # Plot("counts") already called above; data is ready.
        data_all = []
        model_all = []
        for si in range(1, AllData.nSpectra + 1):
            # Plot.y(si) = observed counts,
            # Plot.model(si) = folded model counts
            obs = np.array(Plot.y(si))
            mod = np.array(Plot.model(si))
            data_all.append(obs)
            model_all.append(mod)

        data_total = np.concatenate(data_all)
        model_total = np.concatenate(model_all)

        # Remove channels with zero model prediction
        # (ignored channels or zero-exposure bins)
        mask = model_total > 0
        data_total = data_total[mask]
        model_total = model_total[mask]

        if len(data_total) > 0:
            ks, ks_pval = _ks_bootstrap(data_total, model_total)
            gof["ks_stat"] = float(ks)
            gof["ks_pvalue"] = float(ks_pval)

        logger.info(
            f"   Goodness-of-fit: cstat/dof = "
            f"{gof['cstat']:.1f}/{gof['dof']}, "
            f"KS p-value = {gof['ks_pvalue']:.4f}")

    except Exception as e:
        logger.warning(
            f"   Could not compute goodness-of-fit: {e}")

    return gof


def _ks_2samp(sample1, sample2):
    """KS statistic between two count arrays (CDF comparison)."""
    cdf1 = sample1.cumsum() / sample1.sum()
    cdf2 = sample2.cumsum() / sample2.sum()
    return np.abs(cdf1 - cdf2).max()


def _ks_bootstrap(data, model, niter=1000):
    """KS test with permutation-based p-value.

    Standard KS p-values are invalid when the model was fitted
    to the data. Bootstrap permutation corrects for this by
    estimating the null distribution of the KS statistic under
    random relabelling of data/model counts.

    See Babu & Rao (2004) and
    https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
    """
    ks_obs = _ks_2samp(data, model)

    rng = np.random.default_rng()
    count = 0
    for _ in range(niter):
        # Fixed-label permutation (swap data↔model per bin)
        mask = rng.choice([False, True], size=len(data))
        perm1 = np.where(mask, data, model)
        perm2 = np.where(mask, model, data)
        if _ks_2samp(perm1, perm2) >= ks_obs:
            count += 1

    return ks_obs, count / niter


# ============================================================
# Quality flag definitions (flags 0–11)
# ============================================================
# 0  No issues detected
# 1  Zero or negative source counts
# 2  Zero or negative net counts
# 3  Could not create merged spectrum / background fit failed
# 4  BXA solver failed (timeout, crash, no chain)
# 5  Poor goodness-of-fit (KS p-value < 0.01)
# 6  PhoIndex pegged at prior boundary
# 7  Poor GoF + PhoIndex pegged
# 8  nH pegged at upper boundary
# 9  Poor GoF + nH pegged
# 10 PhoIndex pegged + nH pegged
# 11 Poor GoF + PhoIndex pegged + nH pegged
# ============================================================

# Thresholds for quality flag checks
GOF_PVALUE_THRESHOLD = 0.01
PHOINDEX_PEG_MARGIN = 0.05   # within 0.05 of [1.0, 3.0]
NH_PEG_THRESHOLD = 100.0     # = 10^24 cm^-2 (Compton-thick);
                             # prior cap is now 1000.0 (10^25).


def _compute_quality_flag(results, model_name, logger):
    """Compute post-fit quality flag (5–11) from posterior.

    Called after a successful BXA fit (flag=0). Checks for:
    - Poor goodness-of-fit (KS p-value < 0.01)
    - PhoIndex pegged at prior limits (powerlaw only)
    - nH pegged at upper boundary

    These three issues combine into flags 5–11 via:
        combined = gof_bad + 2*pho_pegged + 4*nh_pegged
        flag = 4 + combined  (if combined > 0, else 0)

    Returns the flag value (0 if no issues, 5–11 otherwise).
    """
    labels = results.get("parameter_names", [])
    medians = results.get("posterior_median", [])
    param_dict = dict(zip(labels, medians))

    # 1. Goodness-of-fit check
    ks_p = results.get("ks_pvalue", np.nan)
    gof_bad = (np.isfinite(ks_p) and ks_p < GOF_PVALUE_THRESHOLD)

    # 2. PhoIndex pegged check (powerlaw model only)
    pho_pegged = False
    if model_name == "powerlaw":
        pho_med = param_dict.get("PhoIndex", None)
        if pho_med is not None:
            pho_pegged = (pho_med <= 1.0 + PHOINDEX_PEG_MARGIN
                          or pho_med >= 3.0 - PHOINDEX_PEG_MARGIN)

    # 3. nH pegged check
    nh_pegged = False
    nh_med = param_dict.get("nH", None)
    if nh_med is not None:
        nh_pegged = (nh_med >= NH_PEG_THRESHOLD)

    combined = int(gof_bad) + 2 * int(pho_pegged) + 4 * int(nh_pegged)
    flag = (4 + combined) if combined > 0 else 0

    if flag > 0:
        issues = []
        if gof_bad:
            issues.append(f"poor GoF (KS p={ks_p:.4f})")
        if pho_pegged:
            issues.append(
                f"PhoIndex pegged ({pho_med:.3f})")
        if nh_pegged:
            issues.append(
                f"nH pegged ({nh_med:.3f})")
        logger.warning(
            f"   Quality flag={flag}: "
            + ", ".join(issues))

    return flag


def fit_spectrum_bxa(spectrum_files, background_files, rmf_files, arf_files,
                     redshift=0.0, model_name="powerlaw",
                     output_base="bxa_fit_results", srcid="unknown", log_file="fit_spectrum_bxa.log",
                     prefit=True, instruments=None,
                     use_iin=True, iin_bounds=(0.5, 2.0),
                     observed_flux=True):
    # Defensive: ensure redshift is a valid float
    if redshift is None:
        redshift = 0.0

    logger.info('\n')
    logger.info(f'Starting BXA fit on {len(spectrum_files)} spectrum file(s)')
    dirname = os.path.dirname(spectrum_files[0])
    os.chdir(dirname)
    logger.info(f'   Changing focus to {dirname}')

    # Background check already done by the caller (automated_fits.py)
    AllData.clear()
    AllModels.clear()

    Fit.statMethod = "cstat"
    Plot.device = "/null"

    # Load each spectrum into its own data group so that the
    # IIN constant can be configured per instrument.
    # Syntax: "group:group filename" assigns a spectrum to a
    # specific data group in XSPEC.
    n_spectra = len(spectrum_files)
    if n_spectra == 1:
        # Single spectrum: standard loading (no IIN needed)
        s = Spectrum(spectrum_files[0])
        s.background = background_files[0]
        s.response = rmf_files[0]
        s.response.arf = arf_files[0]
    else:
        # Multiple spectra: load into separate data groups
        # so each gets its own model copy (for IIN constant)
        for i in range(n_spectra):
            gi = i + 1  # 1-indexed group
            AllData(f"{gi}:{gi} {spectrum_files[i]}")
            AllData(gi).background = background_files[i]
            AllData(gi).response = rmf_files[i]
            AllData(gi).response.arf = arf_files[i]

    AllData.ignore("**-0.3 10.0-**")

    # Disable IIN for single-spectrum fits (no cross-calibration)
    effective_iin = use_iin and n_spectra > 1

    # Model + priors. prefit=False yields the original WIDE
    # prior ranges and runs ultranest without LM bootstrapping
    # (used for the validation experiment comparing prefit vs
    # no-prefit posteriors; see EXECUTIVE_SUMMARY.md §validation).
    if not prefit:
        logger.info("   PREFIT DISABLED for this fit "
                    "(--no_prefit). Using original WIDE priors. "
                    "Expect 10-100x slower BXA convergence.")
    model, priors_list = get_model_and_priors(
        model_name, redshift, prefit=prefit,
        n_groups=n_spectra, instruments=instruments,
        use_iin=effective_iin, iin_bounds=iin_bounds,
        observed_flux=observed_flux)

    # Output dir
    timestamp = datetime.datetime.now().strftime("%d%m%Y_%H%M")
    model_dirname = f"{model_name}_{timestamp}"
    output_dir = os.path.abspath(os.path.join(output_base, str(srcid), model_dirname))
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f'   Setting the output directory for the fits {output_dir}')

    # Run BXA
    # Number of free parameters determines optimal settings
    n_free = len(priors_list)
    solver = bxa.BXASolver(
        transformations=priors_list,
        outputfiles_basename=os.path.join(output_dir))

    # Speed-tuned run:
    # - n_live_points: 35 per free parameter (min 100);
    #   sufficient for smooth unimodal posteriors like
    #   absorbed powerlaw (3 params → 105 live points)
    # - evidence_tolerance=1.0: posteriors converge well
    #   before evidence integral; safe for parameter
    #   estimation (use 0.5 for model comparison)
    # - speed="safe": use BXA's default step sampler
    # - Lepsilon=0.1: default likelihood tolerance
    # - SIGALRM timeout: 30 min hard cap to prevent
    #   pathological sources from blocking a worker
    BXA_TIMEOUT = 1800  # seconds (30 minutes)
    n_live = max(35 * n_free, 100)
    logger.info(
        f"   BXA run: {n_free} free params, "
        f"{n_live} live points, dlogz=1.0")

    # Snapshot open fds BEFORE the fit so we can detect leaks
    fds_before = _get_open_fds()

    # Set wallclock alarm so pathological fits don't run forever
    bxa_failed = False
    old_handler = signal.signal(signal.SIGALRM, _bxa_timeout_handler)
    signal.alarm(BXA_TIMEOUT)
    try:
        solver.run(
            resume=False,
            n_live_points=n_live,
            evidence_tolerance=1.0,
            speed="safe",
            Lepsilon=0.1)
    except BXATimeout:
        logger.warning(
            f"   BXA fit timed out after {BXA_TIMEOUT}s "
            f"for source {srcid}")
        bxa_failed = True
    except Exception as e:
        logger.error(
            f"   BXA solver failed for {srcid}: {e}")
        bxa_failed = True
    finally:
        signal.alarm(0)  # cancel the alarm
        signal.signal(signal.SIGALRM, old_handler)

        # --- Thorough fd cleanup (MUST run regardless) ---
        # Save paramnames before releasing the solver
        labels = solver.paramnames

        # 1. Close ultranest's HDF5 pointstore explicitly
        try:
            sampler = getattr(solver, 'solver', None)
            if sampler is not None:
                ps = getattr(sampler, 'pointstore', None)
                if ps is not None and hasattr(ps, 'close'):
                    ps.close()
        except Exception:
            pass

        # 2. Delete the solver object
        del solver

        # 3. Close ultranest's logger FileHandlers (it creates
        #    a new debug.log handler per run, never closes them)
        for lname in ('ultranest', 'ultranest.integrator',
                      'ultranest.mlfriends',
                      'ultranest.stepsampler'):
            ul = logging.getLogger(lname)
            for h in ul.handlers[:]:
                if isinstance(h, logging.FileHandler):
                    ul.removeHandler(h)
                    h.close()

        # 4. FD fence: force-close leaked file descriptors.
        #    Only closes fds that are (a) new since our snapshot
        #    and (b) safe to close (not library caches like
        #    matplotlib fonts or shared objects).
        fds_after = _get_open_fds()
        leaked = fds_after - fds_before
        if leaked:
            n_closed = 0
            for fd in leaked:
                if _is_safe_to_close(fd):
                    try:
                        os.close(fd)
                        n_closed += 1
                    except OSError:
                        pass
            if n_closed:
                logger.warning(
                    f"   Closed {n_closed} leaked fds "
                    f"(of {len(leaked)} new)")

    # If solver.run() failed, clean up XSPEC and return
    if bxa_failed:
        AllData.clear()
        AllModels.clear()
        plt.close('all')
        gc.collect()
        return {"flag": 4}

    # Read chain, make plots, return summaries
    chain_file = os.path.join(output_dir, "chain.fits")
    if os.path.exists(chain_file):
        with fits.open(chain_file) as hdul:
            samples = hdul[1].data
            samples_array = np.column_stack([samples[name] for name in samples.names])
            stds = np.std(samples_array, axis=0)
            valid_cols = stds > 0
            filtered_samples = samples_array[:, valid_cols]
            filtered_labels = [label for i, label in enumerate(labels) if valid_cols[i]]

            if filtered_samples.shape[1] > 0:
                fig = corner.corner(filtered_samples, labels=filtered_labels, show_titles=True, title_fmt=".3e")
                corner_path = os.path.join(output_dir, "corner.png")
                fig.savefig(corner_path)
                plt.close(fig)
                logger.info(f'   Saved corner plot to file {corner_path} ')

        posterior_median = np.median(samples_array, axis=0)
        posterior_p16    = np.percentile(samples_array, 16, axis=0)
        posterior_p84    = np.percentile(samples_array, 84, axis=0)

        # ----- Goodness-of-fit statistics -----
        # Compute while XSPEC model + data are still loaded.
        # Set model parameters to posterior medians so the
        # folded model reflects the Bayesian best estimate.
        gof = _compute_goodness_of_fit(
            labels, posterior_median, samples_array, logger)

        AllData.clear()
        AllModels.clear()
        plt.close('all')
        gc.collect()

        result = {
            "parameter_names": labels,
            "posterior_median": posterior_median,
            "posterior_p16": posterior_p16,
            "posterior_p84": posterior_p84,
            "output_dir": output_dir,
            "flag": 0,
            "n_spectra": len(spectrum_files),
            "instruments": instruments,
        }
        result.update(gof)

        # ----- Post-fit quality flag (5–11) -----
        quality_flag = _compute_quality_flag(
            result, model_name, logger)
        if quality_flag > 0:
            result["flag"] = quality_flag

        return result

    else:
        logger.error(f'   Chain file {chain_file} not found after BXA run ')
        AllData.clear()
        AllModels.clear()
        plt.close('all')
        gc.collect()
        return {"flag": 4}


def export_bxa_results_to_fits(srcid, output_base="bxa_fit_results", fits_filename="fit_results.fits", log_file="fit_spectrum_bxa.log", global_results=False):
    # directory containing the fit results
    src_dir = os.path.join(output_base, str(srcid))
    #
    # generating filename to save the summary fit results to
    if (global_results):
        # adding them to a file including all previous results
        fits_path = os.path.join(output_base, fits_filename)
    else:
        # adding/writing them to a file in the directory with the fit results
        fits_path = os.path.join(src_dir, fits_filename)
    #
    # print('\n\n Inside export...')
    # print(f'    output_base=({output_base})')
    # print(f'    fits_path=({fits_path})')
    # print(f'    src_dir=({src_dir})')
    
    logger.info('\n')
    logger.info(f'Exporting BXA fit results to FITS file {fits_path}')

    os.makedirs(output_base, exist_ok=True)

    # Get all models in SRCID directory
    #    sorted so that, for each model, the last in the list is the latest fit
    model_dirs = sorted( [d for d in os.listdir(src_dir) if os.path.isdir(os.path.join(src_dir, d))] )

    short_map = {"powerlaw": "PL", "blackbody": "BB", "bremss": "BR", "apec_single": "AP"}

    model_data = {}
    for mdir in model_dirs:
        for long_name, short_name in short_map.items():
            if mdir.startswith(long_name):
                model_chain = os.path.join(src_dir, mdir, "chain.fits")
                if os.path.exists(model_chain):
                    with fits.open(model_chain) as hdul:
                        samples = hdul[1].data
                        param_names = hdul[1].columns.names
                        medians = [np.percentile(samples[p], 50) for p in param_names]
                        p16 = [np.percentile(samples[p], 16) for p in param_names]
                        p84 = [np.percentile(samples[p], 84) for p in param_names]
                        model_data[short_name] = {
                            "names": param_names,
                            "medians": medians,
                            "p16": p16,
                            "p84": p84
                        }

    # Create or update table without Models column
    if os.path.exists(fits_path):
        table = Table.read(fits_path)
        src_mask = table["SRCID"] == srcid
        if np.any(src_mask):
            # table exists and data for srcid are already there, updating them
            idx = np.where(src_mask)[0][0]
            for model_short, pdata in model_data.items():
                for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                    col_med = f"{model_short}_{pname}_median"
                    col_p16 = f"{model_short}_{pname}_p16"
                    col_p84 = f"{model_short}_{pname}_p84"
                    for col, val in zip([col_med, col_p16, col_p84], [median, p16, p84]):
                        if col not in table.colnames:
                            table[col] = np.full(len(table), np.nan)
                        table[col][idx] = val
        else:
            # table exists, but data for srcid are not there yet, adding them
            new_row = {"SRCID": srcid}
            for model_short, pdata in model_data.items():
                for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                    new_row[f"{model_short}_{pname}_median"] = median
                    new_row[f"{model_short}_{pname}_p16"] = p16
                    new_row[f"{model_short}_{pname}_p84"] = p84
            for col in new_row.keys():
                if col not in table.colnames:
                    table[col] = np.full(len(table), np.nan)
            table.add_row(new_row)
        table.write(fits_path, overwrite=True)
    else:
        # table does not exist, creating it with the data for srcid
        row_data = {"SRCID": [srcid]}
        for model_short, pdata in model_data.items():
            for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                row_data[f"{model_short}_{pname}_median"] = [median]
                row_data[f"{model_short}_{pname}_p16"] = [p16]
                row_data[f"{model_short}_{pname}_p84"] = [p84]
        Table(row_data).write(fits_path, overwrite=True)


def export_bxa_results_to_fits_bulk(
        accumulated_results, output_base="bxa_fit_results",
        fits_filename="fit_results.fits",
        partial_flush_every=10000):
    """
    Write all BXA fit results in one pass — O(N) not O(N²).

    Instead of reading-modifying-writing the FITS table for
    each source (which becomes catastrophically slow at 818K
    sources), this function builds the entire table in memory
    and writes it once.

    For crash safety on very long runs, partial files are
    flushed every `partial_flush_every` sources.

    Parameters:
    - accumulated_results: list of (srcid, results_dict)
        where results_dict has keys 'parameter_names',
        'posterior_median', 'posterior_p16', 'posterior_p84'.
    - output_base (str): Top output directory.
    - fits_filename (str): Output FITS filename.
    - partial_flush_every (int): Flush partial file every
        N sources for crash safety.
    """
    logger.info(
        f'\nBulk exporting {len(accumulated_results)} '
        f'BXA results to {fits_filename}')

    os.makedirs(output_base, exist_ok=True)
    fits_path = os.path.join(output_base, fits_filename)

    # Collect all unique column names first
    all_columns = set()
    for srcid, results in accumulated_results:
        if results is None:
            continue
        for pname in results.get("parameter_names", []):
            all_columns.add(f"{pname}_median")
            all_columns.add(f"{pname}_p16")
            all_columns.add(f"{pname}_p84")

    # Build rows
    rows = []
    for i, (srcid, results) in enumerate(
            accumulated_results):
        if results is None:
            continue
        row = {"SRCID": srcid}
        pnames = results.get("parameter_names", [])
        medians = results.get("posterior_median", [])
        p16s = results.get("posterior_p16", [])
        p84s = results.get("posterior_p84", [])

        for pname, med, lo, hi in zip(
                pnames, medians, p16s, p84s):
            row[f"{pname}_median"] = float(med)
            row[f"{pname}_p16"] = float(lo)
            row[f"{pname}_p84"] = float(hi)

        # Goodness-of-fit columns (present if pipeline
        # version includes _compute_goodness_of_fit)
        for gof_key in ("cstat", "dof", "ks_stat",
                         "ks_pvalue"):
            if gof_key in results:
                row[gof_key] = results[gof_key]

        # Quality flag (0 = clean, 5–11 = quality warnings)
        if "flag" in results:
            row["flag"] = int(results["flag"])

        # Number of spectra fitted simultaneously
        if "n_spectra" in results:
            row["n_spectra"] = int(results["n_spectra"])

        # Instrument identifiers per data group, e.g. "pn,MOS1,MOS2"
        # Allows mapping factor columns to specific instruments.
        inst_list = results.get("instruments")
        if inst_list is not None:
            row["instruments"] = ",".join(inst_list)

        # Fill missing columns with NaN
        for col in all_columns:
            if col not in row:
                row[col] = np.nan

        rows.append(row)

        # Partial flush for crash safety
        if (partial_flush_every > 0 and
                len(rows) % partial_flush_every == 0):
            part_num = len(rows) // partial_flush_every
            partial_path = os.path.join(
                output_base,
                f"{os.path.splitext(fits_filename)[0]}"
                f"_partial_{part_num:04d}.fits")
            _write_rows_to_fits(rows, partial_path)
            logger.info(
                f"   Partial flush: {len(rows)} rows "
                f"→ {partial_path}")

    if len(rows) == 0:
        logger.warning(
            "   No valid results to export in bulk mode")
        return

    # Final write
    _write_rows_to_fits(rows, fits_path)
    logger.info(
        f'   Bulk export complete: {len(rows)} rows '
        f'→ {fits_path}')


def _write_rows_to_fits(rows, fits_path):
    """Helper to write a list of row dicts to a FITS table.

    Uses atomic write (temp file + rename) so that a crash
    during the write cannot destroy previously flushed data.
    """
    if len(rows) == 0:
        return

    # Collect all column names
    all_cols = set()
    for row in rows:
        all_cols.update(row.keys())

    # Build column arrays
    # Detect string columns (any row has a str value for that key)
    string_cols = set()
    for row in rows:
        for k, v in row.items():
            if isinstance(v, str):
                string_cols.add(k)

    col_data = {}
    for col in sorted(all_cols):
        if col == "SRCID":
            col_data[col] = [row.get(col, 0) for row in rows]
        elif col in string_cols:
            col_data[col] = [row.get(col, "") for row in rows]
        else:
            col_data[col] = [
                row.get(col, np.nan) for row in rows]

    table = Table(col_data)

    # Atomic write: write to a temp file first, then rename.
    # If the write crashes, the original file is untouched.
    tmp_path = fits_path + ".tmp"
    table.write(tmp_path, format='fits', overwrite=True)
    os.replace(tmp_path, fits_path)

