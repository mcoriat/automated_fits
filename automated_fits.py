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
#    parser.add_argument("--bxa_output_dir", type=str, default="bxa_fit_results",
#                    help="Directory to store BXA results (default: bxa_fit_results)")
    parser.add_argument("--export_results_fits", action="store_true", help="Export BXA fit results to updated FITS file")
    parser.add_argument("--export_filename", default="fit_results.fits", help="Optional output filename for exported FITS results")
    args = parser.parse_args()
    #

    srcid = int(args.srcid)
    message = ''
    #
    output_dir=os.path.abspath(args.output_dir)
    if not os.path.exists(output_dir):
        message = f'\n Creating directory {output_dir}'
        os.mkdir(output_dir)
    #
    src_dir = os.path.join(output_dir, f'{srcid}')
    if not os.path.exists(src_dir):
        message += f'\n Creating subdirectory {src_dir}'
        os.mkdir(src_dir)
    #
    log_file = os.path.join(src_dir, f"{srcid}_process_log_{args.model_name}.txt")
    open(log_file, 'w').close()
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode='w', force=True)
    if message:
        logger.warning(message)
    #
    logger.info(f'\n\n The initial name of this logger file is {log_file}')
    
    
    # Writing out to log file the status of all the input arguments
    logger.info('\n\n Values of all the input arguments for this run')
    print('\n\n Values of all the input arguments for this run')
    for arg in vars(args):
        logger.info(f'    {arg} = {getattr(args,arg)} ')
        print(f'    {arg} = {getattr(args,arg)} ')


    # Starting processing
    logger.info(f'\n\n\n Working on SRCID {srcid} ')
    
    
    # getting from the catalogue the list of obs_id,srcnum that correspond to that srcid
    logger.info(f'\n\n Getting the OBS_ID and SRCNUM associated to SRCID={srcid} in file {args.catalog}')
    srcid_obsid_mapping = read_stacked_catalog(args.catalog, srcid)
    n_mapping=len(srcid_obsid_mapping)
    print(f'\n\n   srcid_obsid_mapping={srcid_obsid_mapping}')
    if  n_mapping== 0:
        logger.error(f"   SRCID {srcid} not found ")
        sys.exit(1)
    else:
        n_pairs=len(srcid_obsid_mapping[srcid])
        if (n_pairs==0):
            logger.error(f'    There are not OBS_ID,SRC_NUM pairs for SRCID={srcid}')
            sys.exit(2)
        else:
            logger.info(f'    {n_pairs} OBS_ID,SRC_NUM pairs found for SRCID={srcid} ')
            logger.info(f'        {srcid_obsid_mapping}')


    # finding which of those combinations actually correspond to existing spectra on disk
    logger.info('\n\n Finding which of those combinations correspond to existing spectra on disk')
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, args.data_dir, log_file)
    nspec = len(srcid_list_spectra)
    if (nspec>0):
        logger.info(f'   {nspec} spectra found for SRCID {srcid}')
        logger.info(f'       {srcid_list_spectra}')
    else:
        logger.error(f' No spectra found for SRCID {srcid}')
        sys.exit(3)


    # finding which of those spectra correspond to each instrument, and checking their counts
    logger.info('\n\n Finding which of those spectra are suitable for fitting')
    pn_list, mos_list = check_spectra(srcid_list_spectra, args.responses_dir, src_dir, log_file)
    #
    pn_good=[spec_tuple for spec_tuple in pn_list if spec_tuple[5] == 0]
    npn=len(pn_list)
    npn_good=len(pn_good)
    logger.info(f'    pn: {npn} spectra found, of which {npn_good} are suitable for fitting')
    #
    mos_good=[spec_tuple for spec_tuple in mos_list if spec_tuple[5] == 0]
    nmos=len(mos_list)
    nmos_good=len(mos_good)
    logger.info(f'   MOS: {nmos} spectra found, of which {nmos_good} are suitable for fitting')
    #
    if (npn==0):
        # no pn spectra
        pn_flag=-1
    elif (npn==1 and npn_good==0):
        # no good spectra, getting the flag from the "best" (but still bad) spectrum
        pn_best=pn_list[0]
        pn_flag=pn_best[5]       
    else:
        # all good so far
        pn_flag=0
    #
    if (nmos==0):
        # no mos spectra
        mos_flag=-1
    elif (nmos==1 and nmos_good==0):
        # no good spectra, getting the flag from the "best" (but still bad) spectrum
        mos_best=mos_list[0]
        mos_flag=mos_best[5]
    else:
        # all good so far
        mos_flag=0
    #
    if (npn_good+nmos_good==0):
        logger.error(f' No spectra suitable for fitting found for SRCID {srcid}')
        logger.error(f'    pn_flag={pn_flag}  mos_flag={mos_flag}')
        # provisional solution
        # ideally a line for this ID with the data for this spectrum should be written
        #    to an output file, to be merged into an output catalogue of all spectral fits
        sys.exit(4)
        
    
    # merging the spectra, if there is more than one for any instrument
    #    merge_spectra checks for empty lists, or lists only with not good spectra
    logger.info('\n\n Selecting which spectra to merge for each instrument and merging them')
    merged_list = merge_spectra(pn_list, mos_list, srcid, src_dir, log_file, mincts=1)
    #
    merged_list_pn=[spec_tuple for spec_tuple in merged_list if spec_tuple[7]=='pn']
    merged_list_pn_good=[spec_tuple for spec_tuple in merged_list_pn if spec_tuple[5]==0]
    npn_merged=len(merged_list_pn)
    npn_merged_good=len(merged_list_pn_good)
    logger.info(f'    pn: {npn_merged} spectra merged, of which {npn_merged_good} are suitable for fitting')
    if (npn_merged_good==0 and pn_flag==0): pn_flag==3
    #
    merged_list_mos=[spec_tuple for spec_tuple in merged_list if spec_tuple[7]=='MOS']
    merged_list_mos_good=[spec_tuple for spec_tuple in merged_list_mos if spec_tuple[5]==0]
    nmos_merged=len(merged_list_mos)
    nmos_merged_good=len(merged_list_mos_good)
    logger.info(f'   MOS: {nmos_merged} spectra merged, of which {nmos_merged_good} are suitable for fitting')
    if (nmos_merged_good==0 and mos_flag==0): mos_flag==3
    #
    #
    if (npn_merged_good+nmos_merged_good==0):
        logger.error(f' No merged spectra suitable for fitting found for SRCID {srcid}')
        logger.error(f'    pn_flag={pn_flag}  mos_flag={mos_flag}')
        # provisional solution
        # ideally a line for this ID with the data for this spectrum should be written
        #    to an output file, to be merged into an output catalogue of all spectral fits
        sys.exit(5)
    


    # selecting only the good spectra, after merging
    fit_list = [spec_tuple for spec_tuple in merged_list if spec_tuple[5] == 0]


    # now going for the fits
    #
    logger.info('\n\n Starting the fits')
    all_fit_results = []
    if args.use_bxa:
        for spec in fit_list:
            # getting from the dictionary in the tuple the filenames
            sp_dic=spec[8]
            spectrum_file = sp_dic['SPECFILE']
            bkg_file = sp_dic['BACKFILE']
            rmf_file = sp_dic['RESPFILE']
            arf_file = sp_dic['ANCRFILE']
            #
            #print("\n   Starting BXA fit for:", spectrum_file)
            results = fit_spectrum_bxa(
                spectrum_file, bkg_file, rmf_file, arf_file,
                redshift=args.redshift,
                model_name=args.model_name,
                srcid=srcid,
                output_base=output_dir,
                log_file=log_file
            )
            # analysing results
            if (results['flag']==0):
                message="\n\nFit completed successfully "
                all_fit_results.append(results)
                for par, med, p16, p84 in zip(results["parameter_names"], results["posterior_median"], results["posterior_p16"], results["posterior_p84"]):
                    message+=f"\n   {par} (median [percentiles 16,84]): {med:.3e} [{p16:.3e},{p84:.3e}]"
                print(message)
                logger.info(message)
            else:
                message=f'\n\nFit failed with flag={results["flag"]} '
                print(message)
                logger.error(message)
                sys.exit(6)


        # already called from within fit_spectru_bxa
        #if args.export_results_fits:
            #print('   Calling export_bxa_results... from automated_fits.py')
            #export_bxa_results_to_fits(srcid, output_dir, args.export_filename)

    else:
        perform_spectrum_fitting(args, srcid, log_file, fit_list, output_dir)

    print('\n\n automated_fits.py finished \n\n')
    logger.info('\n\n automated_fits.py finished \n\n')

if __name__ == "__main__":
    main()

