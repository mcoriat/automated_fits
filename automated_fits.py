import argparse
import os
import glob
import numpy as np
from astropy.io import fits
import bxa.xspec as bxa
import xspec

DIR_DR10 = "/dataAGN/xmmdata/DR10_Spec_bin"
DIR_SCRIPTS = "/dataAGN/xmmdata/xmm2athena_wp6/src"
LOG_FILE = "process_log.txt"

def write_log(message):
    with open(LOG_FILE, 'a') as log:
        log.write(message + '\n')
    print(message)

def read_catalog(catalog_file):
    with fits.open(catalog_file) as hdul:
        catalog_data = hdul[1].data
    return catalog_data

def find_obsids(srcid):
    obsid_paths = glob.glob(f"{DIR_DR10}/{srcid}/*SRSPEC*.FTZ")
    obsids = set()
    for path in obsid_paths:
        obsid = os.path.basename(path).split('_')[0]
        obsids.add(obsid)
    return list(obsids)

def check_spectra(srcid, obsids):
    good_spectra = []
    for obsid in obsids:
        spectrum_file = f"{DIR_DR10}/{srcid}/{obsid}_SRSPEC0001.FTZ"
        background_file = f"{DIR_DR10}/{srcid}/{obsid}_BGSPEC0001.FTZ"
        
        if not os.path.exists(spectrum_file):
            write_log(f"SRCID {srcid}: OBSID {obsid} - Missing source spectrum")
            continue
        
        if not os.path.exists(background_file):
            write_log(f"SRCID {srcid}: OBSID {obsid} - Missing background spectrum")
            continue
        
        with fits.open(background_file) as bg_hdul:
            bg_counts = bg_hdul[1].data['COUNTS'].sum()
        
        if bg_counts > 0:
            with fits.open(spectrum_file) as sp_hdul:
                sp_counts = sp_hdul[1].data['COUNTS'].sum()
            
            if sp_counts > 0:
                snr = sp_counts / np.sqrt(sp_counts + bg_counts)
                good_spectra.append((obsid, sp_counts, bg_counts, snr))
                write_log(f"SRCID {srcid}: OBSID {obsid} - SNR {snr:.2f} (Good for fitting)")
            else:
                write_log(f"SRCID {srcid}: OBSID {obsid} - Source spectrum has zero counts")
        else:
            write_log(f"SRCID {srcid}: OBSID {obsid} - Background spectrum has zero counts")
            
    return good_spectra

def fit_with_bxa(srcid, spectrum_file, background_file, model_name, redshift, use_galabs, use_tbabs_table, output_dir="bxa"):
    try:
        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Load the spectrum and background
        xspec.AllData.clear()
        xspec.AllData(spectrum_file)
        spectrum = xspec.AllData(1)
        spectrum.background = background_file

        # Define the model
        if model_name == "powerlaw":
            model = xspec.Model("powerlaw")
        elif model_name == "blackbody":
            model = xspec.Model("bbody")
        elif model_name == "apec_single":
            model = xspec.Model("apec")
        elif model_name == "apec_apec":
            model = xspec.Model("apec + apec")
        elif model_name == "apec_apec_const":
            model = xspec.Model("apec + apec + constant")
        elif model_name == "bremss":
            model = xspec.Model("bremss")
        elif model_name == "powerlaw_blackbody":
            model = xspec.Model("powerlaw + bbody")
        elif model_name == "powerlaw_blackbody_const":
            model = xspec.Model("powerlaw + bbody + constant")
        elif model_name == "zpowlaw":
            model = xspec.Model(f"zpowerlw,zpowerlw.redshift={redshift}")
        elif model_name == "double_zpowlaw":
            model = xspec.Model(f"zpowerlw + zpowerlw,zpowerlw.redshift={redshift}")
        else:
            raise ValueError(f"Unknown model name: {model_name}")

        # Set additional options
        if use_galabs:
            model.setPars(galabs=True)
        if use_tbabs_table:
            model.setPars(tbabs=True)

        # Configure BXA
        fit = bxa.Fit(model, output_dir)
        fit.run()

        # Save the results
        fit.results(output_dir + "/fit_results.fits")
    except Exception as e:
        write_log(f"SRCID {srcid}: BXA fitting failed for OBSID {obsid} with model {model_name}. Error: {str(e)}")

def fit_spectrum(srcid, spectra, args, model_name):
    best_spectrum = max(spectra, key=lambda x: x[3])  # Choosing based on highest SNR
    obsid = best_spectrum[0]
    spectrum_file = f"{DIR_DR10}/{srcid}/{obsid}_SRSPEC0001.FTZ"
    background_file = f"{DIR_DR10}/{srcid}/{obsid}_BGSPEC0001.FTZ"
    
    write_log(f"Using BXA to fit {model_name} model for SRCID {srcid} OBSID {obsid}")
    fit_with_bxa(srcid, spectrum_file, background_file, model_name, args.redshift, args.use_galabs, args.use_tbabs_table)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog", help="Path to the catalog FITS file")
    parser.add_argument("output", help="Path to the output file for results")
    parser.add_argument("--init", action="store_true", help="initialize the directory")
    parser.add_argument("--combine", action="store_true", help="re-merge the spectra")
    parser.add_argument("--fit_bkg", action="store_true", help="fit the background model")
    parser.add_argument("--get_bkg_stat", action="store_true", help="get the bkg statistics")

    parser.add_argument("--fit_pl", action="store_true", help="fit powerlaw")
    parser.add_argument("--fit_bb", action="store_true", help="fit blackbody")
    parser.add_argument("--fit_apec_single", action="store_true", help="fit apec_singe")
    parser.add_argument("--fit_apec_apec", action="store_true", help="fit apec_apec")
    parser.add_argument("--fit_apec_apec_const", action="store_true", help="fit apec_apec_const")
    parser.add_argument("--fit_bremss", action="store_true", help="fit bremss")
    parser.add_argument("--fit_bbpl", action="store_true", help="fit blackbody powerlaw")
    parser.add_argument("--fit_bbpl_const", action="store_true", help="fit blackbody powerlaw with constant normalization")
    parser.add_argument("--fit_bbpl_const2", action="store_true", help="fit blackbody powerlaw with constant normalization")
    parser.add_argument("--fit_zpl", action="store_true", help="fit redshifted powerlaw")
    parser.add_argument("--fit_zplpl", action="store_true", help="fit redshifted double powerlaw")

    parser.add_argument("--suffix", default=None, help="directory suffix to add to srcid on dataAGN. If left out, it is found automatically based on the best-fit parameters. Use --suffix='' to override this behavior.")
    parser.add_argument("--suffix2", default=None, help="directory suffix to add for det_there != det_use")
    parser.add_argument("--redshift", type=float, help="redshift to be used for redshifted models")
    parser.add_argument("--modelname", help="use a custom modelname i.e. src0001_{modelname}")
    parser.add_argument("--use_galabs", type=int, default=0, help="fix foreground galactic absorption for extragalactic sources. The galactic absorption is estimated automatically according to the source coordinates.")
    parser.add_argument("--use_tbabs_table", type=int, default=0, help="Use tbabs table instead of xspec tbabs")
    parser.add_argument("--spinfo", type=int, default=0, help="calculate spinfo")
    parser.add_argument("--use_xmmsas_version_new", action="store_true", help="use the new (wrong) xmmsas version")
    parser.add_argument("--overwrite", type=int, default=1, help="overwrite existing model")
    args = parser.parse_args()

    catalog_data = read_catalog(args.catalog)
    results = []

    for src in catalog_data:
        srcid = src['SRCID']
        write_log(f"Processing SRCID {srcid}")
        obsids = find_obsids(srcid)
        if not obsids:
            write_log(f"SRCID {srcid}: No OBSIDs found")
            continue
        good_spectra = check_spectra(srcid, obsids)
        if not good_spectra:
            write_log(f"SRCID {srcid}: No good spectra found for fitting")
            continue

        # Fit each model type requested
        if args.fit_pl:
            fit_spectrum(srcid, good_spectra, args, "powerlaw")
        if args.fit_bb:
            fit_spectrum(srcid, good_spectra, args, "blackbody")
        if args.fit_apec_single:
            fit_spectrum(srcid, good_spectra, args, "apec_single")
        if args.fit_apec_apec:
            fit_spectrum(srcid, good_spectra, args, "apec_apec")
        if args.fit_apec_apec_const:
            fit_spectrum(srcid, good_spectra, args, "apec_apec_const")
        if args.fit_bremss:
            fit_spectrum(srcid, good_spectra, args, "bremss")
        if args.fit_bbpl:
            fit_spectrum(srcid, good_spectra, args, "powerlaw_blackbody")
        if args.fit_bbpl_const:
            fit_spectrum(srcid, good_spectra, args, "powerlaw_blackbody_const")
        if args.fit_bbpl_const2:
            fit_spectrum(srcid, good_spectra, args, "powerlaw_blackbody_const2")
        if args.fit_zpl:
            fit_spectrum(srcid, good_spectra, args, "zpowlaw")
        if args.fit_zplpl:
            fit_spectrum(srcid, good_spectra, args, "double_zpowlaw")

        results.append((srcid, good_spectra))
    
    with open(args.output, 'w') as f:
        for result in results:
            srcid, spectra = result
            f.write(f"SRCID: {srcid}\n")
            for spectrum in spectra:
                f.write(f"  OBSID: {spectrum[0]}, Source Counts: {spectrum[1]}, Background Counts: {spectrum[2]}, SNR: {spectrum[3]:.2f}\n")

if __name__ == "__main__":
    # Initialize log file
    with open(LOG_FILE, 'w') as log:
        log.write("Process Log\n")
        log.write("="*40 + "\n")

    main()

