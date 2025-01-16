import os
import logging
import numpy as np
from list_spectra import list_spectra
from check_spectra import check_spectra
from merge_spectra import merge_spectra
from spectral_fitting import perform_spectrum_fitting

def main():
    # Configuration
    srcid = "3067718060100029"
    data_dir = "./test_data"
    output_dir = "./test_data/tests"
    log_file = f"{output_dir}/tests.log"
    redshift = 1.0
    overwrite = True

    # Set up logging
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file)
        ]
    )
    logger = logging.getLogger(__name__)

    logger.info("Starting automated fits script...")

    # Step 1: List spectra
    logger.info("Listing spectra...")
    srcid_obsid_mapping = {
        srcid: [
            {"OBS_ID": "0760940101", "SRC_NUM": 11}
        ]
    }
    spectra_details = list_spectra(srcid, srcid_obsid_mapping, data_dir)
    if not spectra_details:
        logger.error("No spectra found for the given SRCID.")
        return

    # Step 2: Check spectra
    logger.info("Checking spectra...")
    pn_spectra, mos_spectra = check_spectra(spectra_details, log_file)

    # Step 3: Merge spectra
    logger.info("Merging spectra...")
    merged_spectra = merge_spectra(pn_spectra, mos_spectra, srcid, output_dir, log_file)
    logger.info(f"Merged spectra: {merged_spectra}")

    # Step 4: Process merged spectra
    logger.info("Processing merged spectra...")

    if len(merged_spectra) == 0:
        logger.error("No spectra available for processing.")
        return

    # Initialize variables for detection status
    if len(merged_spectra) == 1:
        # Handle single instrument
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

        if merged_spectra[0]["instrument"] == "pn":
            pn_dict = merged_spectra[0]
            mos_dict = {**nan_dict, "instrument": "MOS"}
            det_there = 0
            det_use = 0 if pn_dict["flag"] == 0 else -1
        else:
            mos_dict = merged_spectra[0]
            pn_dict = {**nan_dict, "instrument": "pn"}
            det_there = 1
            det_use = 1 if mos_dict["flag"] == 0 else -1
    else:
        # Handle both instruments
        pn_dict = merged_spectra[0]
        mos_dict = merged_spectra[1]
        det_there = 2
        if pn_dict["flag"] == 0 and mos_dict["flag"] == 0:
            det_use = 2
        elif pn_dict["flag"] == 0:
            det_use = 0
        elif mos_dict["flag"] == 0:
            det_use = 1
        else:
            det_use = -1

    logger.info(f"Detection status: det_there={det_there}, det_use={det_use}")

    # Step 5: Perform spectral fitting if applicable
    if det_use >= 0:
        logger.info("Performing spectrum fitting...")
        perform_spectrum_fitting(pn_dict, mos_dict, redshift, output_dir, overwrite)
    else:
        logger.warning("No suitable spectra for fitting.")

    logger.info("Automated fits script completed.")

if __name__ == "__main__":
    main()

