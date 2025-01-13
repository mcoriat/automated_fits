import sys
import os

# Add the directory containing `spectral_fitting.py` to the Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Now import from spectral_fitting
from spectral_fitting import perform_spectrum_fitting, fit_spectrum

# Mock arguments
class MockArgs:
    redshift = 0.1
    use_galabs = True
    use_tbabs_table = False
    fit_pl = True
    fit_bb = True
    fit_apec_single = False
    fit_apec_apec = False
    fit_apec_apec_const = False
    fit_bremss = False
    fit_bbpl = False
    fit_bbpl_const = False
    fit_bbpl_const2 = False
    fit_zpl = False
    fit_zplpl = False


# Mock data for spectra
good_spectra = [
    {
        "spectrum_file": "./test_data/spectrum1.fits",
        "sp_counts": 1000,
        "bg_counts": 50,
        "sp_netcts": 950,
        "sp_exp": 5000,
        "flag": 0,
        "snr": 15.2,
        "instrument": "pn"
    },
    {
        "spectrum_file": "./test_data/spectrum2.fits",
        "sp_counts": 1200,
        "bg_counts": 60,
        "sp_netcts": 1140,
        "sp_exp": 6000,
        "flag": 0,
        "snr": 18.5,
        "instrument": "mos"
    }
]

from unittest.mock import patch, MagicMock

def test_fit_spectrum():
    srcid = 12345
    args = MockArgs()
    output_dir = "./test_results"
    log_file = "./test_results/log.txt"

    # Mock fits.open to bypass the need for real FITS files
    with patch("spectral_fitting.fits.open", return_value=MagicMock()) as mock_fits:
        fit_spectrum(srcid, good_spectra, args, "powerlaw", log_file, output_dir)

    # Assertions to verify the output directory and model folder
    assert os.path.exists(output_dir), "Output directory was not created"
    assert os.path.exists(os.path.join(output_dir, "powerlaw")), "Model output directory missing"

def test_perform_spectrum_fitting():
    srcid = 12345
    args = MockArgs()
    output_dir = "./test_results"
    log_file = "./test_results/log.txt"

    perform_spectrum_fitting(args, srcid, log_file, good_spectra, output_dir)
    assert os.path.exists(output_dir), "Output directory was not created"

