'''

Wrapping script to accept a SRCID for a 5XMM source in a given catalogue, find the associated spectra, fit them with the chosen model
  and write the fit results to the specified file.

All output written to{output_dir}/{srcid}/ , which is created if it does not exit

Logger file is {output_dir}/{srcid}/{srcid}_process_log_{model_name}.txt , which is overwritten if it exists

Meaning of the flags:
    -2 : cannot open spectral file
    -1 : cannot open background file
     0 : no issues detected
     1 : zero or negative source counts
     2 : zero or negative source counts (which also implies <=0 net counts)
     3 : could not create merged spectrum / background fit failed
     4 : BXA solver failed (timeout, crash, no chain)
     5 : poor goodness-of-fit (KS p-value < 0.01)
     6 : PhoIndex pegged at prior boundary
     7 : poor GoF + PhoIndex pegged
     8 : nH pegged at upper boundary
     9 : poor GoF + nH pegged
    10 : PhoIndex pegged + nH pegged
    11 : poor GoF + PhoIndex pegged + nH pegged
    Note: flags 5-11 indicate quality warnings (fit completed,
    posteriors are valid but should be reviewed)

Output error codes
    1 : SRCID not present in the catalogue
    2 : SRCID present in the catalogue, but there are no associated OBS_ID,SRC_NUM
    3 : SRCID present in the catalogue, but the associated OBS_ID,SRC_NUM have no associated extracted spectra
    4 : SRCID present in the catalogue, the associated OBS_ID,SRC_NUM have associated extracted spectra,
          but they are not valid (cannot be opened, background file not present or cannot be opened, arf/rmf files not present)
    5 : SRCID present in the catalogue, the associated OBS_ID,SRC_NUM have associated extracted spectra,
          some individual spectra are valid, but no corresponding merged spectra could be created
    6 : Got all the way to the fitting stage, but fit failed to produce a chain file


'''
import argparse
import os
import sys
import gc
import shutil
import logging

import numpy as np


def _safe_print(*args, **kwargs):
    """print() wrapper that survives broken stdout (e.g. NFS glitch
    killing the log file descriptor).  Falls back to logger."""
    try:
        print(*args, **kwargs)
    except OSError:
        # stdout fd is dead — log instead so we don't lose
        # the message entirely
        msg = " ".join(str(a) for a in args)
        logging.getLogger().info(f"(stdout dead) {msg}")


# ── File-descriptor fence helpers ──────────────────────────
# These catch fd leaks from ANY library (XSPEC/cfitsio, BXA,
# ultranest, matplotlib, h5py) by comparing the set of open
# fds before and after processing a source.

def _get_open_fds():
    """Return set of open fd numbers (Linux /proc/self/fd)."""
    try:
        return set(int(fd) for fd in os.listdir('/proc/self/fd'))
    except OSError:
        return set()


# Extensions / prefixes of fds that must NOT be closed —
# they are long-lived library caches, not per-source leaks.
_FD_KEEP_SUFFIXES = (
    '.ttf', '.otf', '.woff', '.woff2',   # matplotlib fonts
    '.so', '.dylib',                       # shared libraries
    '.pyc', '.pyo',                        # compiled Python
)
_FD_KEEP_PREFIXES = (
    '/dev/', '/proc/', '/sys/',            # OS virtual fs
    'pipe:', 'socket:', 'anon_inode:',     # non-file fds
)


def _is_safe_to_close(fd):
    """Check if a leaked fd can be safely closed.
    Returns True for fit-related files (HDF5, FITS, logs),
    False for library caches (fonts, .so) and OS resources."""
    try:
        target = os.readlink(f'/proc/self/fd/{fd}')
    except OSError:
        return False
    for pfx in _FD_KEEP_PREFIXES:
        if target.startswith(pfx):
            return False
    for sfx in _FD_KEEP_SUFFIXES:
        if target.endswith(sfx):
            return False
    return True


from read_stacked_catalog import (read_stacked_catalog,
                                  read_stacked_catalog_batch)
from list_spectra import (list_spectra,
                          build_dir_listing_cache)
from check_spectra import check_spectra
from merge_spectra import merge_spectra
from spectral_fitting import perform_spectrum_fitting
from spectral_fitting_bxa_adapted import (
    fit_spectrum_bxa,
    export_bxa_results_to_fits,
    export_bxa_results_to_fits_bulk,
    check_background_fit
)


logger = logging.getLogger(__name__)


def _load_completed_srcids(output_dir):
    """Scan output_dir for fit_results*.fits files and return
    the set of SRCIDs that already have results.

    Used by --skip_completed when chain files have been cleaned
    up (--cleanup_chains) and are no longer on disk.
    """
    from astropy.io import fits as pyfits
    import glob as globmod
    completed = set()
    pattern = os.path.join(output_dir, "fit_results*.fits")
    for fpath in globmod.glob(pattern):
        try:
            with pyfits.open(fpath) as hdul:
                if len(hdul) > 1 and 'SRCID' in hdul[1].columns.names:
                    for sid in hdul[1].data['SRCID']:
                        completed.add(int(sid))
        except Exception:
            pass  # skip corrupt / incomplete files
    return completed


def _load_previous_results(output_dir):
    """Load existing fit results from fit_results*.fits files
    as (srcid, results_dict) tuples suitable for accumulated_results.

    Used on --resume to pre-seed accumulated_results so that the
    final export preserves results from a previous (interrupted) run.
    Without this, the final export overwrites the FITS file with
    only the new results, losing the earlier ones.
    """
    from astropy.table import Table
    import glob as globmod
    previous = []
    pattern = os.path.join(output_dir, "fit_results*.fits")
    for fpath in sorted(globmod.glob(pattern)):
        try:
            t = Table.read(fpath)
        except Exception:
            continue
        colnames = t.colnames
        # Detect parameter names from *_median columns
        param_names = [c.replace("_median", "")
                       for c in colnames
                       if c.endswith("_median") and c != "SRCID"]
        for row in t:
            srcid = int(row["SRCID"])
            results = {
                "parameter_names": param_names,
                "posterior_median": np.array(
                    [float(row[f"{p}_median"]) for p in param_names]),
                "posterior_p16": np.array(
                    [float(row[f"{p}_p16"]) for p in param_names]),
                "posterior_p84": np.array(
                    [float(row[f"{p}_p84"]) for p in param_names]),
                "flag": int(row["flag"]) if "flag" in colnames else 0
            }
            # Restore GoF columns if present
            for key in ("cstat", "dof", "ks_stat", "ks_pvalue"):
                if key in colnames:
                    val = row[key]
                    results[key] = float(val) if key != "dof" else int(val)
            previous.append((srcid, results))
    return previous


def _cleanup_chain_dir(results, srcid, output_dir):
    """Remove the BXA output directory (chain.fits, corner.png,
    ultranest run files) after summary statistics have been
    extracted into memory.  Keeps the per-source directory and
    the log file.

    Parameters:
    - results: dict returned by fit_spectrum_bxa, must contain
        'output_dir' key pointing to the model subdirectory.
    - srcid: Source ID (for logging).
    - output_dir: Top-level output directory (unused here but
        kept for consistency).
    """
    chain_dir = results.get("output_dir")
    if chain_dir and os.path.isdir(chain_dir):
        try:
            shutil.rmtree(chain_dir)
            logger.info(
                f"   Cleaned up chain directory: {chain_dir}")
            print(
                f"   Cleaned up chain directory: {chain_dir}")
        except OSError as e:
            logger.warning(
                f"   Could not clean up {chain_dir}: {e}")


def process_one_source(srcid, args, output_dir,
                       srcid_obsid_mapping=None,
                       dir_listing_cache=None,
                       completed_srcids=None):
    """
    Process a single SRCID through the full pipeline:
    catalog lookup → list spectra → validate → merge →
    rebin → fit → (optionally export).

    Parameters:
    - srcid (int): The SRCID to process.
    - args: Parsed argparse namespace with all CLI options.
    - output_dir (str): Absolute path to the top output dir.
    - srcid_obsid_mapping (dict or None): Pre-loaded catalog
        mapping from read_stacked_catalog_batch(). If None,
        the catalog is read for this single source.
    - dir_listing_cache (dict or None): Pre-built directory
        listing cache from build_dir_listing_cache(). If None,
        directories are listed on the fly.
    - completed_srcids (set or None): Pre-loaded set of SRCIDs
        already present in fit_results*.fits files. Used by
        --skip_completed when chains have been cleaned up.

    Returns:
    - (error_code, results_dict_or_None)
      error_code: 0 on success, 1-6 on various failures,
                  3 for background failure (flag=3).
      results_dict_or_None: BXA fit results dict when
                  successful, None otherwise.
    """
    # Per-source output directory
    src_dir = os.path.join(output_dir, f'{srcid}')
    if not os.path.exists(src_dir):
        os.makedirs(src_dir, exist_ok=True)

    # Skip if already completed: check chain.fits on disk
    # OR presence in fit_results*.fits tables.
    if getattr(args, 'skip_completed', False):
        # Check 1: chain.fits still on disk
        import glob
        chain_pattern = os.path.join(
            src_dir, f"{args.model_name}*", "chain.fits")
        existing = glob.glob(chain_pattern)
        if existing:
            print(f"   Skipping SRCID {srcid} "
                  f"(chain.fits already exists)")
            return (0, None)
        # Check 2: already in results FITS table
        if completed_srcids and srcid in completed_srcids:
            print(f"   Skipping SRCID {srcid} "
                  f"(already in fit_results)")
            return (0, None)

    # Per-source log file
    log_file = os.path.join(
        src_dir,
        f"{srcid}_process_log_{args.model_name}.txt")
    open(log_file, 'w').close()

    # Set up a file handler for this source
    # (remove previous handlers to avoid cross-contamination)
    root_logger = logging.getLogger()
    # Remove all existing file handlers (may already be closed
    # by the outer fd fence in the batch loop)
    for h in root_logger.handlers[:]:
        if isinstance(h, logging.FileHandler):
            root_logger.removeHandler(h)
            try:
                h.close()
            except Exception:
                pass
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.INFO)
    root_logger.addHandler(fh)
    root_logger.setLevel(logging.INFO)

    logger.info(
        f'\n\n The log file for this source is {log_file}')

    # Log arguments
    logger.info(
        '\n\n Values of all the input arguments for this run')
    for arg in vars(args):
        logger.info(f'    {arg} = {getattr(args, arg)} ')

    # Start
    message = f'\n\n\n Working on SRCID {srcid} '
    logger.info(message)
    print(message)

    # Catalogue
    message = (f'\n\n Getting the OBS_ID and SRCNUM associated '
               f'to SRCID={srcid} in file {args.catalog}')
    logger.info(message)
    print(message)

    if srcid_obsid_mapping is None:
        # Single-source mode: read catalog just for this SRCID
        srcid_obsid_mapping = read_stacked_catalog(
            args.catalog, srcid)

    n_mapping = len(srcid_obsid_mapping)
    if n_mapping == 0 or srcid not in srcid_obsid_mapping:
        logger.error(f"   SRCID {srcid} not found \n\n")
        print(f"\n\n   ERROR 1: SRCID={srcid} not found \n\n")
        return (1, None)
    else:
        n_pairs = len(srcid_obsid_mapping[srcid])
        if n_pairs == 0:
            logger.error(
                f'    There are not OBS_ID,SRC_NUM pairs '
                f'for SRCID={srcid} \n\n')
            print(
                f'\n\n    ERROR 2: There are not OBS_ID,'
                f'SRC_NUM pairs for SRCID={srcid} \n\n')
            return (2, None)
        else:
            logger.info(
                f'    {n_pairs} OBS_ID,SRC_NUM pairs found '
                f'for SRCID={srcid} ')
            logger.info(
                f'        {srcid_obsid_mapping[srcid]}')
            print(
                f'    {n_pairs} OBS_ID,SRC_NUM pairs found '
                f'for SRCID={srcid} ')

    # Spectra
    logger.info(
        '\n\n Finding which of those combinations '
        'correspond to existing spectra on disk')
    srcid_list_spectra = list_spectra(
        srcid, srcid_obsid_mapping, args.data_dir,
        dir_cache=dir_listing_cache,
        subdir=args.subdir)
    nspec = len(srcid_list_spectra)
    if nspec > 0:
        logger.info(
            f'   {nspec} spectra found for SRCID {srcid}')
        logger.info(f'       {srcid_list_spectra}')
        print(f'   {nspec} spectra found for SRCID {srcid}')
    else:
        logger.error(
            f' No spectra found for SRCID {srcid}\n\n')
        print(f'\n\n    ERROR 3: No spectra found for '
              f'SRCID={srcid}\n\n')
        return (3, None)

    # Good spectra
    logger.info(
        '\n\n Finding which of those spectra are suitable '
        'for fitting')
    print('\n\n Finding which of those spectra are suitable '
          'for fitting')
    pn_list, mos_list = check_spectra(
        srcid_list_spectra, args.responses_dir,
        src_dir, log_file)

    pn_good = [t for t in pn_list if t[5] == 0]
    npn = len(pn_list)
    npn_good = len(pn_good)
    logger.info(
        f'    pn: {npn} spectra found, of which '
        f'{npn_good} are suitable for fitting')
    print(
        f'    pn: {npn} spectra found, of which '
        f'{npn_good} are suitable for fitting')

    mos_good = [t for t in mos_list if t[5] == 0]
    nmos = len(mos_list)
    nmos_good = len(mos_good)
    logger.info(
        f'   MOS: {nmos} spectra found, of which '
        f'{nmos_good} are suitable for fitting')
    print(
        f'   MOS: {nmos} spectra found, of which '
        f'{nmos_good} are suitable for fitting')

    if npn_good + nmos_good == 0:
        logger.error(
            f' No spectra suitable for fitting found '
            f'for SRCID {srcid}')
        print(
            f'\n\n    ERROR 4: No spectra suitable for '
            f'fitting found for SRCID {srcid}\n\n')
        # Determine source-level flag from per-spectrum flags.
        # flag 1 = zero/negative counts, flag 2 = zero/negative
        # net counts, anything else = generic pre-fit failure.
        all_spec_flags = [t[5] for t in pn_list + mos_list]
        count_flags = [f for f in all_spec_flags if f in (1, 2)]
        if count_flags:
            src_flag = max(count_flags)  # 2 supersedes 1
        else:
            src_flag = 3  # file access or other issue
        return (4, {"flag": src_flag})

    # Merge
    logger.info(
        '\n\n Selecting which spectra to merge for each '
        'instrument and merging them')
    print(
        '\n\n Selecting which spectra to merge for each '
        'instrument and merging them')
    merged_list = merge_spectra(
        pn_list, mos_list, srcid,
        src_dir, log_file, mincts=1)

    merged_list_good = [
        t for t in merged_list if t[5] == 0]
    if len(merged_list_good) == 0:
        logger.error(
            f' No merged spectra suitable for fitting '
            f'found for SRCID {srcid}')
        print(
            f'\n\n    ERROR 5: No merged spectra suitable '
            f'for fitting found for SRCID {srcid}\n\n')
        return (5, {"flag": 3})

    fit_list = merged_list_good

    # Fits
    logger.info('\n\n Starting the fits')
    print('\n\n Starting the fits')

    if args.use_bxa:
        if len(fit_list) == 0:
            logger.error(
                f"   No good spectra available for "
                f"simultaneous fitting for source {srcid}")
            return (5, {"flag": 3})

        # === Background check (optional) ===
        if not args.skip_bkg_check:
            valid_specs = []
            for spec in fit_list:
                sp_dic = spec[8]
                spectrum_file = sp_dic['SPECFILE']
                bkg_file = sp_dic['BACKFILE']
                rmf_file = sp_dic['RESPFILE']
                arf_file = sp_dic['ANCRFILE']

                inst = spec[7] if len(spec) > 7 else "unknown"

                logger.info(
                    f"   Checking background for instrument "
                    f"{inst} (SRCID {srcid})")
                pval = check_background_fit(
                    spectrum_file, bkg_file,
                    rmf_file, arf_file,
                    args.model_name, args.redshift,
                    logger, srcid=str(srcid))
                if pval is None:
                    logger.warning(
                        f"   Background NOT OK for {inst} "
                        f"→ skipping this instrument")
                else:
                    valid_specs.append(spec)

            if not valid_specs:
                logger.warning(
                    f"   No valid instruments left after "
                    f"background check → flagging SRCID "
                    f"{srcid} as 3")
                print(
                    f"   Skipping source {srcid} (flag=3)")
                return (3, {"flag": 3})
            fit_list = valid_specs
        # else: skip_bkg_check is True, use fit_list as-is

        spectrum_files = []
        background_files = []
        rmf_files = []
        arf_files = []
        instruments = []

        for spec in fit_list:
            sp_dic = spec[8]
            spectrum_files.append(sp_dic['SPECFILE'])
            background_files.append(sp_dic['BACKFILE'])
            rmf_files.append(sp_dic['RESPFILE'])
            arf_files.append(sp_dic['ANCRFILE'])
            instruments.append(
                spec[7] if len(spec) > 7 else "unknown")

        logger.info(
            f"   Performing simultaneous BXA fit for "
            f"source {srcid} using "
            f"{len(spectrum_files)} spectra:")
        print(
            f"   Performing simultaneous BXA fit for "
            f"source {srcid} using "
            f"{len(spectrum_files)} spectra:")

        for sf, inst in zip(spectrum_files, instruments):
            logger.info(f"      Spectrum file: {sf} ({inst})")
            print(f"      Spectrum file: {sf} ({inst})")

        results = fit_spectrum_bxa(
            spectrum_files=spectrum_files,
            background_files=background_files,
            rmf_files=rmf_files,
            arf_files=arf_files,
            redshift=args.redshift,
            model_name=args.model_name,
            srcid=srcid,
            output_base=output_dir,
            log_file=log_file,
            prefit=not getattr(args, 'no_prefit', False),
            instruments=instruments,
            use_iin=not getattr(args, 'no_iin', False),
            iin_bounds=(
                getattr(args, 'iin_lo', 0.5),
                getattr(args, 'iin_hi', 1.5)),
            observed_flux=not getattr(
                args, 'intrinsic_flux', False)
        )

        # Background failure from fit_spectrum_bxa
        if results.get("flag", 0) == 3:
            logger.warning(
                f"   Skipping source {srcid} due to "
                f"background failure (flag=3).")
            print(
                f"   Skipping source {srcid} due to "
                f"background failure (flag=3).")
            return (3, {"flag": 3})

        flag = results.get("flag", 0)

        # Flags 5–11: fit completed with quality warnings
        #   (posteriors are valid, just flagged for review).
        # Flag 4: BXA solver failed (no chain / no posterior).
        if flag == 0 or flag >= 5:
            message = "\n\nFit completed successfully"
            if flag >= 5:
                message += f" (quality flag={flag})"
            for par, med, p16, p84 in zip(
                    results["parameter_names"],
                    results["posterior_median"],
                    results["posterior_p16"],
                    results["posterior_p84"]):
                message += (
                    f"\n   {par} (median [percentiles "
                    f"16,84]): {med:.3e} "
                    f"[{p16:.3e},{p84:.3e}]")
            print(message)
            logger.info(message)

            # In single-source mode, export immediately
            # (in batch mode, bulk export is done after all)
            if args.srcid_file is None:
                if args.export_results_fits:
                    export_bxa_results_to_fits(
                        srcid, output_dir,
                        args.export_filename,
                        log_file=log_file,
                        global_results=True)
                # Clean up after export (single-source mode)
                if args.cleanup_chains:
                    _cleanup_chain_dir(
                        results, srcid, output_dir)

            return (0, results)
        else:
            # Flag 4 or unexpected: BXA failed, no posteriors
            logger.info('\n\n')
            message = (
                f'\n\nFit failed with '
                f'flag={results["flag"]} \n\n ')
            logger.error(message)
            print(
                f'\n\n    ERROR 6: Fit failed with '
                f'flag={results["flag"]} \n\n')
            return (6, {"flag": results["flag"]})

    else:
        perform_spectrum_fitting(
            args, srcid, log_file, fit_list, output_dir)
        return (0, None)


def main():
    parser = argparse.ArgumentParser(
        description="Automated spectral fitting pipeline. "
        "Single-source mode: pass SRCID as first argument. "
        "Batch mode: use --srcid_file instead of SRCID.")

    # SRCID: positional, required for single-source mode.
    # For batch mode (--srcid_file), pass 0 as a placeholder.
    parser.add_argument(
        "srcid", type=int,
        help="SRCID to fit (single-source mode), or 0 when "
             "using --srcid_file for batch mode")
    parser.add_argument(
        "data_dir",
        help="Path to the directory containing the data")
    parser.add_argument(
        "script_dir",
        help="Path to the directory containing the scripts")
    parser.add_argument(
        "responses_dir",
        help="Path to the directory containing the "
             "response matrices")
    parser.add_argument(
        "output_dir",
        help="Path to the output directory")
    parser.add_argument(
        "catalog",
        help="Stacked catalog FITS filename (including path)")
    parser.add_argument(
        "output",
        help="Name of the output file with fit results")

    parser.add_argument(
        "--srcid_file", type=str, default=None,
        help="Text file with one SRCID per line (batch mode)")
    parser.add_argument(
        "--subdir", type=str, default="pps",
        help="Subdirectory name under each OBS_ID folder "
             "containing spectra (default: pps). "
             "Use 'product' for 5XMM production data.")
    parser.add_argument(
        "--skip_bkg_check", action="store_true",
        help="Skip background quality check before BXA "
             "fitting")
    parser.add_argument(
        "--skip_completed", action="store_true",
        help="Skip SRCIDs that already have a chain.fits "
             "for the requested model (useful for restarts)")

    parser.add_argument(
        "--init", action="store_true",
        help="initialize the directory")
    parser.add_argument(
        "--combine", action="store_true",
        help="re-merge the spectra")
    parser.add_argument(
        "--fit_bkg", action="store_true",
        help="fit the background model")
    parser.add_argument(
        "--get_bkg_stat", action="store_true",
        help="get the bkg statistics")

    parser.add_argument(
        '--model_name', type=str, default='powerlaw',
        choices=['powerlaw', 'apec_single',
                 'blackbody', 'bremss'],
        help='Spectral model to use with BXA '
             '(default: powerlaw)')

    parser.add_argument(
        "--suffix", default=None,
        help="directory suffix to add to srcid on dataAGN")
    parser.add_argument(
        "--suffix2", default=None,
        help="directory suffix to add for "
             "det_there != det_use")
    parser.add_argument(
        "--redshift", type=float, default=0.0,
        help="redshift to be used for redshifted models "
             "(default: 0.0 = no redshift correction)")
    parser.add_argument(
        "--modelname",
        help="use a custom modelname "
             "i.e. src0001_{modelname}")
    parser.add_argument(
        "--use_galabs", type=int, default=0,
        help="fix foreground galactic absorption")
    parser.add_argument(
        "--use_tbabs_table", type=int, default=0,
        help="Use tbabs table")
    parser.add_argument(
        "--spinfo", type=int, default=0,
        help="calculate spinfo")
    parser.add_argument(
        "--use_xmmsas_version_new", action="store_true",
        help="use new xmmsas version")
    parser.add_argument(
        "--overwrite", type=int, default=1,
        help="overwrite existing model")
    parser.add_argument(
        "--use_bxa", action="store_true",
        help="use BXA fitting instead of XSPEC")
    parser.add_argument(
        "--cleanup_chains", action="store_true",
        help="Delete chain.fits and corner.png after "
             "extracting summary statistics. Saves disk "
             "space at the cost of losing full posteriors.")
    parser.add_argument(
        "--export_results_fits", action="store_true",
        help="Export BXA fit results to updated FITS file")
    parser.add_argument(
        "--export_filename", default="fit_results.fits",
        help="Optional output filename for exported "
             "FITS results")
    parser.add_argument(
        "--no_prefit", action="store_true",
        help="Disable the Levenberg-Marquardt pre-fit step "
             "that tightens BXA priors. Used for the validation "
             "experiment comparing prefit vs no-prefit "
             "posteriors. WARNING: BXA will be 10-100x slower "
             "without prefit, use only on small samples.")
    parser.add_argument(
        "--intrinsic_flux", action="store_true",
        help="Use phabs*cflux*model (intrinsic/unabsorbed "
             "flux, 5XMM legacy convention). Default is "
             "cflux*phabs*model (observed/absorbed flux, "
             "XMM2ATHENA convention per Viitanen+25).")
    parser.add_argument(
        "--no_iin", action="store_true",
        help="Disable the inter-instrument normalisation "
             "(IIN) constant. By default a multiplicative "
             "constant is included: frozen at 1.0 for pn "
             "(reference) and free [0.5, 1.5] for MOS, "
             "following Viitanen+25.")
    parser.add_argument(
        "--iin_lo", type=float, default=0.5,
        help="Lower bound for the IIN constant (default 0.5)")
    parser.add_argument(
        "--iin_hi", type=float, default=1.5,
        help="Upper bound for the IIN constant (default 1.5)")
    args = parser.parse_args()

    # Validate: single-source vs batch mode
    # In batch mode, srcid=0 is used as a placeholder
    if args.srcid_file is not None and args.srcid != 0:
        print("WARNING: --srcid_file given; ignoring "
              f"positional SRCID={args.srcid}. "
              "Use 0 as placeholder in batch mode.")
    if args.srcid_file is None and args.srcid == 0:
        print("ERROR: SRCID=0 is a placeholder for batch "
              "mode. Provide --srcid_file or a real SRCID.")
        sys.exit(1)

    output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # ============================================
    # BATCH MODE
    # ============================================
    if args.srcid_file is not None:
        # Read the list of SRCIDs
        with open(args.srcid_file, 'r') as f:
            srcids = [int(line.strip())
                      for line in f
                      if line.strip()]
        n_total = len(srcids)
        print(f"\n Batch mode: {n_total} SRCIDs to process")

        # Set up a batch-level log
        batch_log = os.path.join(
            output_dir, "batch_processing.log")
        logging.basicConfig(
            filename=batch_log, level=logging.INFO,
            filemode='w', force=True)
        logger.info(
            f"Batch processing {n_total} SRCIDs")

        # Pre-load catalog for all SRCIDs at once
        print(" Loading catalog (vectorized)...")
        catalog_mapping = read_stacked_catalog_batch(
            args.catalog, srcids)
        print(f" Catalog loaded: {len(catalog_mapping)} "
              f"SRCIDs found")

        # Pre-build directory listing cache
        print(" Building directory listing cache...")
        dir_cache = build_dir_listing_cache(
            catalog_mapping, args.data_dir,
            subdir=args.subdir)
        print(f" Cache built: {len(dir_cache)} directories")

        # Pre-load completed SRCIDs from existing results
        # (works even when chain files have been cleaned up)
        completed = set()
        accumulated_results = []
        n_preloaded = 0
        if getattr(args, 'skip_completed', False):
            print(" Scanning existing fit_results*.fits "
                  "for completed SRCIDs...")
            completed = _load_completed_srcids(output_dir)
            if completed:
                print(f" Found {len(completed)} already-"
                      f"completed SRCIDs")
                # Pre-load previous results so the final export
                # preserves them (otherwise they get overwritten)
                accumulated_results = _load_previous_results(
                    output_dir)
                n_preloaded = len(accumulated_results)
                print(f" Pre-loaded {n_preloaded} previous "
                      f"results into accumulator")

        # Process each source
        n_success = 0
        n_fail = 0
        n_flushed = n_preloaded  # results already written to disk
        FLUSH_EVERY = 50       # save results every N fits

        for i, srcid in enumerate(srcids):
            # Save/restore cwd (BXA + XSPEC chdir)
            original_cwd = os.getcwd()

            # FD fence: snapshot BEFORE processing this source.
            # Anything new after process_one_source() is a leak.
            fds_before = _get_open_fds()

            try:
                _safe_print(f"\n[{i+1}/{n_total}] "
                            f"Processing SRCID {srcid}")
                code, results = process_one_source(
                    srcid, args, output_dir,
                    srcid_obsid_mapping=catalog_mapping,
                    dir_listing_cache=dir_cache,
                    completed_srcids=completed)

                # Accumulate ALL sources that returned a
                # results dict (successful fits AND failure
                # stubs with flag + NaN fit columns).
                # Only errors 1-2 (source not in catalogue)
                # return None and are excluded from the output.
                if results is not None:
                    accumulated_results.append(
                        (srcid, results))

                if code == 0:
                    n_success += 1
                    # Clean up chain files to save disk
                    if args.cleanup_chains:
                        _cleanup_chain_dir(
                            results, srcid, output_dir)
                    # Periodic flush: save results to FITS
                    # so they survive restarts and are
                    # visible to --skip_completed
                    if (args.export_results_fits and
                            n_success % FLUSH_EVERY == 0):
                        _safe_print(f"   Flushing {len(accumulated_results)} "
                                    f"results to {args.export_filename}")
                        export_bxa_results_to_fits_bulk(
                            accumulated_results,
                            output_base=output_dir,
                            fits_filename=args.export_filename)
                        n_flushed = len(accumulated_results)
                else:
                    n_fail += 1
                    logger.info(
                        f"SRCID {srcid}: error code {code}")
            except Exception as e:
                n_fail += 1
                logger.error(
                    f"SRCID {srcid}: unhandled exception: "
                    f"{e}")
                _safe_print(f"   ERROR: SRCID {srcid} failed "
                            f"with exception: {e}")
            finally:
                os.chdir(original_cwd)

                # 1. Close the per-source FileHandler so its
                #    fd is released before the fence runs
                root_logger = logging.getLogger()
                for h in root_logger.handlers[:]:
                    if isinstance(h, logging.FileHandler):
                        root_logger.removeHandler(h)
                        try:
                            h.close()
                        except Exception:
                            pass

                # 2. Garbage-collect to release ref-counted fds
                gc.collect()

                # 3. FD fence: force-close any leaked fds.
                #    Catches leaks from XSPEC/cfitsio, BXA,
                #    ultranest, h5py, matplotlib — everything.
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
                        _safe_print(
                            f"   [fd-fence] Closed {n_closed} "
                            f"leaked fds (of {len(leaked)} new)")

        # Final export (writes all results, including
        # any accumulated since the last flush)
        if args.export_results_fits and accumulated_results:
            if len(accumulated_results) > n_flushed:
                _safe_print(f"\n Exporting {len(accumulated_results)} "
                            f"fit results to FITS...")
            export_bxa_results_to_fits_bulk(
                accumulated_results,
                output_base=output_dir,
                fits_filename=args.export_filename)

        total_fits = n_success + n_preloaded
        _safe_print(f"\n\n Batch processing complete: "
                    f"{n_success} new fits, {n_fail} failed "
                    f"out of {n_total} total"
                    f"{f' (+{n_preloaded} pre-loaded from previous run)' if n_preloaded else ''}.\n"
                    f" Total results in FITS: {total_fits}\n")
        logger.info(
            f"Batch complete: {n_success} new + "
            f"{n_preloaded} pre-loaded = {total_fits} total, "
            f"{n_fail} failed out of {n_total}")

    # ============================================
    # SINGLE-SOURCE MODE (fully backward compatible)
    # ============================================
    else:
        srcid = int(args.srcid)

        # Set up logging (same as before)
        src_dir = os.path.join(output_dir, f'{srcid}')
        if not os.path.exists(src_dir):
            os.makedirs(src_dir, exist_ok=True)

        log_file = os.path.join(
            src_dir,
            f"{srcid}_process_log_{args.model_name}.txt")
        open(log_file, 'w').close()
        logging.basicConfig(
            filename=log_file, level=logging.INFO,
            filemode='w', force=True)

        logger.info('\n\n Values of all the input arguments '
                    'for this run')
        print('\n\n Values of all the input arguments '
              'for this run')
        for arg in vars(args):
            logger.info(
                f'    {arg} = {getattr(args, arg)} ')
            print(f'    {arg} = {getattr(args, arg)} ')

        code, results = process_one_source(
            srcid, args, output_dir)
        if code != 0:
            sys.exit(code)

        print('\n\n automated_fits.py finished '
              'successfully \n\n')
        logger.info('\n\n automated_fits.py finished '
                    'successfully \n\n')


if __name__ == "__main__":
    main()
