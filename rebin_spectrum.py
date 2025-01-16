import os
import numpy as np
from astropy.io import fits
import logging

logger = logging.getLogger(__name__)

def rebin_spectrum(infile, outfile, background_file, log_file, mincts=1):
    """
    Rebins an input FITS spectrum file to have >=mincts counts per bin.
    """
    logger.info(f"Rebinning file {infile} with a minimum of {mincts} counts per bin")
    outfile = os.path.abspath(outfile)

    # Get spectral counts
    spec_dict = get_spectral_counts(infile, log_file, background_file)
    logger.info(f"Returned spec_dict in rebin_spectrum: {spec_dict}")

    # Validate spec_dict
    if spec_dict["flag"] < 0:
        logger.error(f"Invalid spec_dict received for {infile}: {spec_dict}")
        return spec_dict

    try:
        with fits.open(infile) as hdul:
            data = hdul[1].data
            counts = data['COUNTS']
            quality = data['QUALITY']
            grouping = data['GROUPING']

            rebinned_counts = []
            rebinned_quality = []
            rebinned_grouping = []

            current_bin = []
            for i, count in enumerate(counts):
                current_bin.append(count)
                if sum(current_bin) >= mincts:
                    rebinned_counts.append(sum(current_bin))
                    rebinned_quality.append(0)  # Mark as good
                    rebinned_grouping.append(-1)  # End of bin
                    current_bin = []

            # Handle remaining counts
            if current_bin:
                rebinned_counts.append(sum(current_bin))
                rebinned_quality.append(2)  # Mark as bad
                rebinned_grouping.append(-1)

            # Create a new FITS HDU for the rebinned data
            new_cols = fits.ColDefs([
                fits.Column(name='CHANNEL', format='I', array=np.arange(1, len(rebinned_counts) + 1)),
                fits.Column(name='COUNTS', format='J', array=np.array(rebinned_counts)),
                fits.Column(name='QUALITY', format='I', array=np.array(rebinned_quality)),
                fits.Column(name='GROUPING', format='I', array=np.array(rebinned_grouping))
            ])
            hdu = fits.BinTableHDU.from_columns(new_cols)

            # Write the new file
            hdu.writeto(outfile, overwrite=True)
            logger.info(f"Rebinned spectrum written to {outfile}")

    except Exception as e:
        logger.error(f"Failed to rebin spectrum file {infile}: {e}")
        spec_dict["flag"] = -1

    return spec_dict



# Function to read in a spectrum and its corresponding background file,
# and return the counts and a flag if any incidence occurred
def get_spectral_counts(infile, log_file, background_file=''):
    """
    Reads in input spectral file infile in FITS format, gets the total counts, and the
    header keywords BACKFILE, BACKSCAL and EXPOSURE.
    Then reads in the background file, gets its total counts and the header
    keyword BACKSCAL

    Then calculates the net counts scaling with the backscale values
    Finally, writes out a dictionary with this information and a flag

    Parameters:
    - infile: input spectrum file in FITS format
    - log_file (str): The log file to write the messages.
    - background_file (str): The name of the background file (see above)

    Returns:
    - spec_dict (dict): A dictionary containing the name of the input spectrum, 
      total source counts, total background counts, total net counts, exposure time, 
      flag, and signal-to-noise ratio.
    """

    # Initializing output values
    spec_dict = {
        "spectrum_file": infile,
        "sp_counts": np.nan,
        "bg_counts": np.nan,
        "sp_netcts": np.nan,
        "sp_exp": np.nan,
        "flag": -2,
        "snr": np.nan
    }

    # Logging input details
    logger.info(f"Processing spectrum file: {infile}")

    # Trying to open the input file
    try:
        with fits.open(infile) as hdul:
            logger.info(f"Opened spectrum file: {infile}")
            logger.info(f"Header keys: {hdul[1].header.keys()}")
            logger.info(f"Data columns: {hdul[1].data.columns.names}")

            spec_dict["sp_counts"] = hdul[1].data['COUNTS'].sum()
            spec_dict["sp_exp"] = hdul[1].header['EXPOSURE']
            sp_backscal = hdul[1].header['BACKSCAL']
            bgd_file = hdul[1].header['BACKFILE']

            # Determine the background file
            if background_file == '':
                # Use BACKFILE header key to construct full path relative to infile directory
                background_file = os.path.join(os.path.dirname(infile), bgd_file)

            logger.info(f"Resolved background file: {background_file}")

            # Check if background file exists
            if not os.path.exists(background_file):
                logger.error(f"Background file not found: {background_file}")
                spec_dict["flag"] = -1
                return spec_dict

            try:
                with fits.open(background_file) as bg_hdul:
                    logger.info(f"Opened background file: {background_file}")

                    spec_dict["bg_counts"] = bg_hdul[1].data['COUNTS'].sum()
                    bg_backscal = bg_hdul[1].header['BACKSCAL']

                    # Calculate net counts
                    spec_dict["sp_netcts"] = spec_dict["sp_counts"] - (spec_dict["bg_counts"] * sp_backscal / bg_backscal)

                    # Set flag values
                    if spec_dict["bg_counts"] <= 0:
                        spec_dict["flag"] = 1
                    elif spec_dict["sp_netcts"] <= 0 or spec_dict["sp_counts"] <= 0:
                        spec_dict["flag"] = 2
                    else:
                        spec_dict["flag"] = 0

                    # Calculate SNR if flag is 0
                    if spec_dict["flag"] == 0:
                        spec_dict["snr"] = spec_dict["sp_netcts"] / np.sqrt(2 * spec_dict["sp_counts"] - spec_dict["sp_netcts"])

            except Exception as e:
                logger.error(f"Failed to open or process background file {background_file}: {e}")
                spec_dict["flag"] = -1

    except Exception as e:
        logger.error(f"Failed to open or process spectrum file {infile}: {e}")
        spec_dict["flag"] = -2

    return spec_dict


def perform_spectrum_fitting(pn_dict, mos_dict, redshift, output_dir, overwrite, args):
    """
    Performs spectral fitting using the input data and arguments.

    Parameters:
    - pn_dict: Dictionary containing pn spectrum data.
    - mos_dict: Dictionary containing MOS spectrum data.
    - redshift: Redshift value for the source.
    - output_dir: Directory to store output.
    - overwrite: Boolean indicating if existing files should be overwritten.
    - args: Dictionary of fitting parameters.
    """
    logger.info("Starting spectrum fitting...")

    if args.get('fit_pl', False):
        logger.info("Performing power-law fitting...")
        # Add power-law fitting logic here

    if args.get('fit_bb', False):
        logger.info("Performing black-body fitting...")
        # Add black-body fitting logic here

    logger.info("Spectrum fitting completed.")


# Function to test get_spectral_counts
def test_get_spectral_counts():
    output_dir = './test_data/test'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Set up logging
    log_file = os.path.join(output_dir, 'test_get_spectral_counts.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)

    # Inexistent spectrum file
    infile = 'dummy_file_FJC.txt'
    spec_dict = get_spectral_counts(infile, log_file)
    assert spec_dict['flag'] == -2

    # Existent spectrum file, getting correct stats
    infile = './test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ'
    spec_dict = get_spectral_counts(infile, log_file, background_file='pps')
    assert spec_dict['sp_counts'] == 766
    assert spec_dict['bg_counts'] == 8474
    assert abs(spec_dict['sp_netcts'] - 271.88) <= 0.01
    assert abs(spec_dict['sp_exp'] - 82181.94) <= 0.01

    print("All tests passed!")

