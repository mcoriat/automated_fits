import os
import numpy as np
from astropy.io import fits
import logging

from get_spectral_counts import get_spectral_counts


logger = logging.getLogger(__name__)



# Function to check which spectra are suitable for fitting
def check_spectra(list_spectra, responses_dir, output_dir, log_file):
    """
    Check the spectra for a given list of spectra and return "good spectra" details in dictionaries.

    Parameters:
    - list_spectra (list): List of spectra file paths.
    - responses_dir (str): Directory with response matrices (rmf) for pn and MOS.
    - output_dir (str): Directory to store symbolic links and good spectra details.
    - log_file (str): Log file to record messages.

    Returns:
    - pn_spectra (list): A list of dictionaries for good pn spectra.
    - mos_spectra (list): A list of dictionaries for good mos spectra.
    """
    pn_spectra = []
    mos_spectra = []

    responses_dir = os.path.abspath(responses_dir)
    output_dir = os.path.abspath(output_dir)

    for spectrum_file in list_spectra:
        logger.info(f"\n\nWorking on file {spectrum_file}")
        spectrum_file = os.path.abspath(spectrum_file)
        name = os.path.basename(spectrum_file)

        # Determine instrument type
        instrument = 'pn' if "PN" in name else 'MOS'

        # Initialize dictionary for storing spectrum details
        spectrum_details = {
            "spectrum_file": spectrum_file,
            "sp_counts": np.nan,
            "bg_counts": np.nan,
            "sp_netcts": np.nan,
            "sp_exp": np.nan,
            "flag": -1,
            "snr": np.nan,
            "instrument": instrument
        }

        # Background file
        background_file = spectrum_file.replace("SRSPEC", "BGSPEC")
        spec_tuple = get_spectral_counts(spectrum_file, log_file, background_file=background_file)

        # Populate the dictionary
        spectrum_details.update({
            "sp_counts": spec_tuple[1],
            "bg_counts": spec_tuple[2],
            "sp_netcts": spec_tuple[3],
            "sp_exp": spec_tuple[4],
            "flag": spec_tuple[5],
            "snr": spec_tuple[6]
        })

        # Check for good spectra
        if spectrum_details["flag"] >= 0 and spectrum_details["sp_netcts"] > 0:
            logger.info(f"Good spectrum found: {spectrum_details}")
            if instrument == 'pn':
                pn_spectra.append(spectrum_details)
            else:
                mos_spectra.append(spectrum_details)
        else:
            logger.info(f"Skipping spectrum: {spectrum_details}")

    return pn_spectra, mos_spectra



def test_check_spectra():
    output_dir = './test_data/test'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    log_file = os.path.join(output_dir, 'test_process_log.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)

    responses_dir = './test_data/RESPONSES'

    try:
        # Test case: only pn spectrum
        print("Testing with only pn spectrum...")
        list_spectra = ['./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ']
        pn_list, mos_list = check_spectra(list_spectra, responses_dir, output_dir, log_file)
        assert len(pn_list) == 1, f"Expected 1 pn spectrum, got {len(pn_list)}"
        assert pn_list[0]["flag"] == 0, "Expected flag 0 for pn spectrum"
        assert len(mos_list) == 0, f"Expected 0 mos spectra, got {len(mos_list)}"
        print("Test 1 passed: pn spectrum only.\n")

        # Test case: only MOS spectra
        print("Testing with only MOS spectra...")
        list_spectra = [
            './test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ',
            './test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ'
        ]
        pn_list, mos_list = check_spectra(list_spectra, responses_dir, output_dir, log_file)
        assert len(pn_list) == 0, f"Expected 0 pn spectra, got {len(pn_list)}"
        assert len(mos_list) == 2, f"Expected 2 mos spectra, got {len(mos_list)}"
        assert mos_list[0]["flag"] == 0, "Expected flag 0 for first MOS spectrum"
        assert mos_list[1]["flag"] == 0, "Expected flag 0 for second MOS spectrum"
        print("Test 2 passed: MOS spectra only.\n")

        # Test case: all spectra
        print("Testing with pn and MOS spectra...")
        list_spectra = [
            './test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ',
            './test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ',
            './test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ'
        ]
        pn_list, mos_list = check_spectra(list_spectra, responses_dir, output_dir, log_file)
        assert len(pn_list) == 1, f"Expected 1 pn spectrum, got {len(pn_list)}"
        assert len(mos_list) == 2, f"Expected 2 mos spectra, got {len(mos_list)}"
        assert pn_list[0]["flag"] == 0, "Expected flag 0 for pn spectrum"
        assert mos_list[0]["flag"] == 0, "Expected flag 0 for first MOS spectrum"
        assert mos_list[1]["flag"] == 0, "Expected flag 0 for second MOS spectrum"
        print("Test 3 passed: pn and MOS spectra.\n")

        print("All tests passed successfully.")

    except AssertionError as e:
        print(f"Test failed: {e}")
        
if __name__ == "__main__":
    test_check_spectra()


