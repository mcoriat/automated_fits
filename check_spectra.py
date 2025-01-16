import os
import numpy as np
from astropy.io import fits
import logging
from get_spectral_counts import get_spectral_counts

logger = logging.getLogger(__name__)

def check_spectra(spectra_details, log_file):
    pn_spectra = []
    mos_spectra = []

    for spectrum in spectra_details:
        spectrum_file = spectrum["spectrum_file"]
        background_file = spectrum["background_file"]

        # Determine instrument type based on file name
        instrument = 'pn' if "PN" in spectrum_file.upper() else 'MOS'

        # Get spectral counts
        spec_dict = get_spectral_counts(spectrum_file, log_file, background_file)
        spec_dict['instrument'] = instrument  # Add the instrument key
        spec_dict['background_file'] = background_file  # Propagate background_file

        # Append to appropriate list if the spectrum is valid
        if spec_dict["flag"] >= 0 and spec_dict["sp_netcts"] > 0:
            if instrument == 'pn':
                pn_spectra.append(spec_dict)
            else:
                mos_spectra.append(spec_dict)

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

