import os
import numpy as np
from check_spectra import check_spectra
from merge_spectra import merge_spectra
from spectral_fitting import perform_spectrum_fitting
import logging

logger = logging.getLogger(__name__)


def main():
    # Placeholder for input arguments (replace with argparse or similar in real implementation)
    srcid = "3067718060100029"
    test_data_dir = "./test_data"
    responses_dir = f"{test_data_dir}/RESPONSES"
    output_dir = f"{test_data_dir}/tests"
    log_file = f"{output_dir}/tests.log"
    catalogue_file = f"{test_data_dir}/test_catalogue.fits"
    redshift = 1.0
    overwrite = True

    logging.basicConfig(filename=log_file, level=logging.INFO)
    logger.info("Starting automated fits script...")

    # Step 1: Check the spectra
    logger.info("Checking spectra...")
    srcid_list_spectra = [
        f"{test_data_dir}/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ",
        f"{test_data_dir}/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ",
        f"{test_data_dir}/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ",
    ]
    pn_list, mos_list = check_spectra(srcid_list_spectra, responses_dir, output_dir, log_file)

    # Step 2: Merge spectra
    logger.info("Merging spectra...")
    merged_list = merge_spectra(pn_list, mos_list, srcid, output_dir, log_file)

    # Step 3: Process merged spectra
    logger.info("Processing merged spectra...")

    if len(merged_list) == 0:
        logger.error("No spectra available for processing.")
        return

    elif len(merged_list) == 1:
        # Only one instrument available
        nan_dict = {
            "spectrum_file": "",
            "sp_counts": np.nan,
            "bg_counts": np.nan,
            "sp_netcts": np.nan,
            "sp_exp": np.nan,
            "flag": -1,
            "snr": np.nan,
            "instrument": ""
        }

        if merged_list[0]["instrument"] == "pn":
            pn_tuple = merged_list[0]
            mos_tuple = {**nan_dict, "instrument": "MOS"}
            det_there = 0
            det_use = 0 if pn_tuple["flag"] == 0 else -1
        else:
            mos_tuple = merged_list[0]
            pn_tuple = {**nan_dict, "instrument": "pn"}
            det_there = 1
            det_use = 1 if mos_tuple["flag"] == 0 else -1
    else:
        # Both instruments available
        pn_tuple = merged_list[0]
        mos_tuple = merged_list[1]
        det_there = 2
        if pn_tuple["flag"] == 0 and mos_tuple["flag"] == 0:
            det_use = 2
        elif pn_tuple["flag"] == 0:
            det_use = 0
        elif mos_tuple["flag"] == 0:
            det_use = 1
        else:
            det_use = -1

    logger.info(f"Detection status: det_there={det_there}, det_use={det_use}")

    # Perform spectrum fitting if applicable
    if det_use >= 0:
        logger.info("Performing spectrum fitting...")
        perform_spectrum_fitting(pn_tuple, mos_tuple, redshift, output_dir, overwrite)
    else:
        logger.warning("No suitable spectra for fitting.")

    logger.info("Automated fits script completed.")


if __name__ == "__main__":
    main()

