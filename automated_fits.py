import argparse
import os
import sys
import numpy as np
import logging
from read_stacked_catalog import read_stacked_catalog
from list_spectra import list_spectra
from check_spectra import check_spectra
from merge_spectra import merge_spectra
from spectral_fitting import perform_spectrum_fitting
from spectral_fitting_bxa_adapted import fit_spectrum_bxa, export_bxa_results_to_fits

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("srcid", type=int, help="SRCID of the source whose spectra are to be fitted")
    parser.add_argument("data_dir", help="Path to the directory containing the data")
    parser.add_argument("script_dir", help="Path to the directory containing the scripts")
    parser.add_argument("responses_dir", help="Path to the directory containing the response matrices")
    parser.add_argument("output_dir", help="Path to the output directory")
    parser.add_argument("catalog", help="Stacked catalog FITS filename (including path)")
    parser.add_argument("output", help="Name of the output file with fit results")

    parser.add_argument("--init", action="store_true", help="initialize the directory")
    parser.add_argument("--combine", action="store_true", help="re-merge the spectra")
    parser.add_argument("--fit_bkg", action="store_true", help="fit the background model")
    parser.add_argument("--get_bkg_stat", action="store_true", help="get the bkg statistics")

    parser.add_argument('--model_name', type=str, default='powerlaw',
                        choices=['powerlaw', 'apec_single', 'blackbody', 'bremss'],
                        help='Spectral model to use with BXA (default: powerlaw)')

    parser.add_argument("--suffix", default=None, help="directory suffix to add to srcid on dataAGN")
    parser.add_argument("--suffix2", default=None, help="directory suffix to add for det_there != det_use")
    parser.add_argument("--redshift", type=float, help="redshift to be used for redshifted models")
    parser.add_argument("--modelname", help="use a custom modelname i.e. src0001_{modelname}")
    parser.add_argument("--use_galabs", type=int, default=0, help="fix foreground galactic absorption")
    parser.add_argument("--use_tbabs_table", type=int, default=0, help="Use tbabs table")
    parser.add_argument("--spinfo", type=int, default=0, help="calculate spinfo")
    parser.add_argument("--use_xmmsas_version_new", action="store_true", help="use new xmmsas version")
    parser.add_argument("--overwrite", type=int, default=1, help="overwrite existing model")
    parser.add_argument("--use_bxa", action="store_true", help="use BXA fitting instead of XSPEC")
    parser.add_argument("--bxa_output_dir", type=str, default="bxa_fit_results",
                    help="Directory to store BXA results (default: bxa_fit_results)")
    parser.add_argument("--export_results_fits", action="store_true", help="Export BXA fit results to updated FITS file")
    parser.add_argument("--export_filename", default="fit_results.fits", help="Optional output filename for exported FITS results")
    args = parser.parse_args()

    srcid = int(args.srcid)
    message = ''
    if not os.path.exists(args.output_dir):
        message = f' Creating directory {args.output_dir}'
        os.mkdir(args.output_dir)
    output_dir = os.path.join(args.output_dir, f'{srcid}')
    if not os.path.exists(output_dir):
        message += f'\n Creating subdirectory {output_dir}'
        os.mkdir(output_dir)

    log_file = os.path.join(output_dir, f"{srcid}_process_log.txt")
    open(log_file, 'w').close()
    logging.basicConfig(filename=log_file, level=logging.INFO)
    if message:
        logger.warning(message)

    logger.info(f'\n\n Working on SRCID {srcid} ')
    srcid_obsid_mapping = read_stacked_catalog(args.catalog, srcid)
    if len(srcid_obsid_mapping) == 0:
        logger.error(f"SRCID {srcid} not found in file {args.catalog}")
        sys.exit(1)

    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, args.data_dir)
    nspec = len(srcid_list_spectra)
    logger.info(f'   {nspec} spectra found for SRCID {srcid}')
    if nspec == 0:
        logger.error(f' No spectra found for SRCID {srcid}')
        sys.exit(2)

    pn_list, mos_list = check_spectra(srcid_list_spectra, args.responses_dir, output_dir, log_file)
    merged_list = merge_spectra(pn_list, mos_list, srcid, output_dir, log_file, mincts=1)

    fit_list = [spec_tuple for spec_tuple in merged_list if spec_tuple[5] == 0]

    all_fit_results = []
    if args.use_bxa:
        for spec in fit_list:
            spectrum_file = str(spec[0])
            bkg_file = str(spec[1])
            rmf_file = str(spec[2])
            arf_file = str(spec[3])
            print("\nStarting BXA fit for:", spectrum_file)
            results = fit_spectrum_bxa(
                spectrum_file, bkg_file, rmf_file, arf_file,
                redshift=args.redshift,
                model_name=args.model_name,
                srcid=srcid,
                output_base=args.bxa_output_dir   
            )

            print("Fit complete.")
            all_fit_results.append(results)
            for par, val in zip(results["parameter_names"], results["posterior_median"]):
                print(f"{par}: {val:.3e}")

        if args.export_results_fits:
            export_bxa_results_to_fits(all_fit_results, args.catalog, srcid, args.export_filename)

    else:
        perform_spectrum_fitting(args, srcid, log_file, fit_list, output_dir)

    logger.info('\n\n automated_fits.py finished \n\n')

if __name__ == "__main__":
    main()

