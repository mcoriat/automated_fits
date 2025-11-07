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
     3 : could not create merged spectrum
     4 : fit failed
  
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
#import numpy as np
import logging
from read_stacked_catalog import read_stacked_catalog
from list_spectra import list_spectra
from check_spectra import check_spectra
from merge_spectra import merge_spectra
from spectral_fitting import perform_spectrum_fitting
from spectral_fitting_bxa_adapted import fit_spectrum_bxa, export_bxa_results_to_fits, check_background_fit


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
    parser.add_argument("--export_results_fits", action="store_true", help="Export BXA fit results to updated FITS file")
    parser.add_argument("--export_filename", default="fit_results.fits", help="Optional output filename for exported FITS results")
    args = parser.parse_args()

    srcid = int(args.srcid)
    message = ''

    output_dir=os.path.abspath(args.output_dir)
    if not os.path.exists(output_dir):
        message = f'\n Creating directory {output_dir}'
        os.mkdir(output_dir)

    src_dir = os.path.join(output_dir, f'{srcid}')
    if not os.path.exists(src_dir):
        message += f'\n Creating subdirectory {src_dir}'
        os.mkdir(src_dir)

    log_file = os.path.join(src_dir, f"{srcid}_process_log_{args.model_name}.txt")
    open(log_file, 'w').close()
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode='w', force=True)
    if message:
        logger.warning(message)

    logger.info(f'\n\n The initial name of this logger file is {log_file}')

    # Arguments
    logger.info('\n\n Values of all the input arguments for this run')
    print('\n\n Values of all the input arguments for this run')
    for arg in vars(args):
        logger.info(f'    {arg} = {getattr(args,arg)} ')
        print(f'    {arg} = {getattr(args,arg)} ')

    # Start
    message=f'\n\n\n Working on SRCID {srcid} '
    logger.info(message)
    print(message)

    # Catalogue
    message=f'\n\n Getting the OBS_ID and SRCNUM associated to SRCID={srcid} in file {args.catalog}'
    logger.info(message)
    print(message)
    srcid_obsid_mapping = read_stacked_catalog(args.catalog, srcid, log_file)
    n_mapping=len(srcid_obsid_mapping)
    if  n_mapping== 0:
        logger.error(f"   SRCID {srcid} not found \n\n")
        print(f"\n\n   ERROR 1: SRCID={srcid} not found \n\n")
        sys.exit(1)
    else:
        n_pairs=len(srcid_obsid_mapping[srcid])
        if (n_pairs==0):
            logger.error(f'    There are not OBS_ID,SRC_NUM pairs for SRCID={srcid} \n\n')
            print(f'\n\n    ERROR 2: There are not OBS_ID,SRC_NUM pairs for SRCID={srcid} \n\n')
            sys.exit(2)
        else:
            logger.info(f'    {n_pairs} OBS_ID,SRC_NUM pairs found for SRCID={srcid} ')
            logger.info(f'        {srcid_obsid_mapping}')
            print(f'    {n_pairs} OBS_ID,SRC_NUM pairs found for SRCID={srcid} ')
            print(f'        {srcid_obsid_mapping}')

    # Spectra
    logger.info('\n\n Finding which of those combinations correspond to existing spectra on disk')
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, args.data_dir, log_file)
    nspec = len(srcid_list_spectra)
    if (nspec>0):
        logger.info(f'   {nspec} spectra found for SRCID {srcid}')
        logger.info(f'       {srcid_list_spectra}')
        print(f'   {nspec} spectra found for SRCID {srcid}')
        print(f'       {srcid_list_spectra}')
    else:
        logger.error(f' No spectra found for SRCID {srcid}\n\n')
        print(f'\n\n    ERROR 3: No spectra found for SRCID={srcid}\n\n')
        sys.exit(3)

    # Good spectra
    logger.info('\n\n Finding which of those spectra are suitable for fitting')
    print('\n\n Finding which of those spectra are suitable for fitting')
    pn_list, mos_list = check_spectra(srcid_list_spectra, args.responses_dir, src_dir, log_file)

    pn_good=[spec_tuple for spec_tuple in pn_list if spec_tuple[5] == 0]
    npn=len(pn_list); npn_good=len(pn_good)
    logger.info(f'    pn: {npn} spectra found, of which {npn_good} are suitable for fitting')
    print(f'    pn: {npn} spectra found, of which {npn_good} are suitable for fitting')

    mos_good=[spec_tuple for spec_tuple in mos_list if spec_tuple[5] == 0]
    nmos=len(mos_list); nmos_good=len(mos_good)
    logger.info(f'   MOS: {nmos} spectra found, of which {nmos_good} are suitable for fitting')
    print(f'   MOS: {nmos} spectra found, of which {nmos_good} are suitable for fitting')

    if (npn_good+nmos_good==0):
        logger.error(f' No spectra suitable for fitting found for SRCID {srcid}')
        print(f'\n\n    ERROR 4: No spectra suitable for fitting found for SRCID {srcid}\n\n')
        sys.exit(4)

    # Merge
    logger.info('\n\n Selecting which spectra to merge for each instrument and merging them')
    print('\n\n Selecting which spectra to merge for each instrument and merging them')
    merged_list = merge_spectra(pn_list, mos_list, srcid, src_dir, log_file, mincts=1)

    merged_list_good=[spec_tuple for spec_tuple in merged_list if spec_tuple[5]==0]
    if len(merged_list_good)==0:
        logger.error(f' No merged spectra suitable for fitting found for SRCID {srcid}')
        print(f'\n\n    ERROR 5: No merged spectra suitable for fitting found for SRCID {srcid}\n\n')
        sys.exit(5)

    fit_list = merged_list_good

    # Fits
    logger.info('\n\n Starting the fits')
    print('\n\n Starting the fits')
    all_fit_results = []

    if args.use_bxa:
        if len(fit_list) == 0:
            logger.error(f"   No good spectra available for simultaneous fitting for source {srcid}")
            sys.exit(5)

        spectrum_files = []
        background_files = []
        rmf_files = []
        arf_files = []
        
        
        # === NEW: pre-filter instruments by unified background gate (flag=3 logic) ===
        valid_specs = []
        for spec in fit_list:
            sp_dic = spec[8]
            spectrum_file = sp_dic['SPECFILE']
            bkg_file      = sp_dic['BACKFILE']
            rmf_file      = sp_dic['RESPFILE']
            arf_file      = sp_dic['ANCRFILE']

            inst = spec[7] if len(spec) > 7 else "unknown"

            logger.info(f"   Checking background for instrument {inst} (SRCID {srcid})")
            pval = check_background_fit(spectrum_file, bkg_file, rmf_file, arf_file,
                                    args.model_name, args.redshift, logger, srcid=str(srcid))
            if pval is None:
                logger.warning(f"   Background NOT OK for {inst} → skipping this instrument")
            else:
                valid_specs.append(spec)

        # Decision logic: if nothing survives, flag source as 3 and stop
        if not valid_specs:
            logger.warning(f"   No valid instruments left after background check → flagging SRCID {srcid} as 3")
            print(f"   Skipping source {srcid} (flag=3)")
            # keep your existing behavior: stop this srcid but continue overall
            return
        # Otherwise, use only valid_specs for fitting
        fit_list = valid_specs


        for spec in fit_list:
            sp_dic = spec[8]
            spectrum_files.append(sp_dic['SPECFILE'])
            background_files.append(sp_dic['BACKFILE'])
            rmf_files.append(sp_dic['RESPFILE'])
            arf_files.append(sp_dic['ANCRFILE'])

        logger.info(f"   Performing simultaneous BXA fit for source {srcid} using {len(spectrum_files)} spectra:")
        print(f"   Performing simultaneous BXA fit for source {srcid} using {len(spectrum_files)} spectra:")

        for sf in spectrum_files:
            logger.info(f"      Spectrum file: {sf}")
            print(f"      Spectrum file: {sf}")


        results = fit_spectrum_bxa(
            spectrum_files=spectrum_files,
            background_files=background_files,
            rmf_files=rmf_files,
            arf_files=arf_files,
            redshift=args.redshift,
            model_name=args.model_name,
            srcid=srcid,
            output_base=output_dir,
            log_file=log_file
        )

        # --- Unified background failure handling ---
        if results.get("flag", 0) == 3:
            logger.warning(f"   Skipping source {srcid} due to background failure (unified flag=3).")
            print(f"   Skipping source {srcid} due to background failure (flag=3).")
            return

        if results.get("flag", 0) == 0:
            message = "\n\nFit completed successfully "
            all_fit_results.append(results)
            for par, med, p16, p84 in zip(results["parameter_names"],
                                          results["posterior_median"],
                                          results["posterior_p16"],
                                          results["posterior_p84"]):
                message += f"\n   {par} (median [percentiles 16,84]): {med:.3e} [{p16:.3e},{p84:.3e}]"
            print(message)
            logger.info(message)
        else:
            logger.info('\n\n')
            message = f'\n\nFit failed with flag={results["flag"]} \n\n '
            logger.error(message)
            print(f'\n\n    ERROR 6: Fit failed with flag={results["flag"]} \n\n')
            sys.exit(6)

        if args.export_results_fits:
            export_bxa_results_to_fits(srcid, output_dir, args.export_filename, log_file=log_file, global_results=True)

    else:
        perform_spectrum_fitting(args, srcid, log_file, fit_list, output_dir)

    print('\n\n automated_fits.py finished successfully \n\n')
    logger.info('\n\n automated_fits.py finished successfully \n\n')

if __name__ == "__main__":
    main()

