import os
import logging
import xspec
import bxa.xspec as bxa
from astropy.io import fits as pf

logger = logging.getLogger(__name__)

def decompress_if_needed(path):
    """
    Decompress .FTZ (fpack-compressed) FITS to .fits for XSPEC compatibility.
    """
    if path and path.lower().endswith('.ftz') and os.path.exists(path):
        new_path = path[:-4] + '.fits'
        try:
            with pf.open(path) as hdulist:
                hdulist.writeto(new_path, overwrite=True)
            logger.info(f"Decompressed {path} to {new_path} using Astropy")
        except Exception as e:
            logger.warning(f"Astropy decompression failed for {path}: {e}. Trying funpack...")
            ret = os.system(f"funpack {path}")
            if ret != 0:
                logger.error(f"funpack failed for {path}, using original path")
                return path
            new_path = path[:-4] + '.fits'
        return new_path
    logger.debug(f"No decompression needed for {path}")
    return path


def perform_spectrum_fitting(srcid, pn_spectra, mos_spectra, args_dict):
    """
    Perform spectrum fitting for the given source ID.

    Parameters:
        srcid (int): The source ID.
        pn_spectra (list): List of PN spectra details.
        mos_spectra (list): List of MOS spectra details.
        args_dict (dict): Dictionary of arguments.

    Returns:
        dict: Results of the fitting.
    """
    logger.info(f"Performing spectrum fitting for SRCID {srcid}...")
    if not pn_spectra and not mos_spectra:
        logger.error(f"No valid spectra available for SRCID {srcid}. Skipping fitting.")
        return {}

    good_spectra = pn_spectra + mos_spectra
    output_dir = args_dict.get("output_dir", ".")
    log_file = args_dict.get("log_file", "fitting.log")

    if args_dict.get('fit_pl', False):
        fit_spectrum(srcid, good_spectra, args_dict, "powerlaw", log_file, output_dir)
    if args_dict.get('fit_bb', False):
        fit_spectrum(srcid, good_spectra, args_dict, "blackbody", log_file, output_dir)
    # Add additional model flags as needed


def fit_spectrum(srcid, spectra, args_dict, model_name, log_file, output_dir):
    """
    Fits a spectrum using the BXA method.

    Parameters:
        srcid (int): The source ID.
        spectra (list): List of spectra dicts with 'spectrum_file' and 'snr'.
        args_dict (dict): Pipeline arguments.
        model_name (str): Name of the model to fit.
        log_file (str): Path to the log file.
        output_dir (str): Base output directory.
    """
    logger.info("\n\nStarting spectrum fitting...")

    # Select the spectrum with highest SNR
    best_spectrum = max(spectra, key=lambda x: x['snr'])
    spectrum_file = best_spectrum['spectrum_file']

    # Read header for background, RMF, and ARF filenames
    try:
        hdul = pf.open(spectrum_file)
        raw_bkg = hdul[1].header.get('BACKFILE', '')
        raw_rmf = hdul[1].header.get('RESPFILE', '')
        raw_arf = hdul[1].header.get('ANCRFILE', '')
        hdul.close()
    except Exception as e:
        logger.error(f"Could not open spectrum {spectrum_file} for SRCID {srcid}. Error: {e}")
        return

    # Construct absolute paths and decompress if needed
    bkg_path = os.path.abspath(os.path.join(os.path.dirname(spectrum_file), raw_bkg))
    bkg_path = decompress_if_needed(bkg_path)
    logger.info(f"Passing background to XSPEC: {bkg_path}")

    rmf_path = os.path.abspath(os.path.join(os.path.dirname(spectrum_file), raw_rmf))
    rmf_path = decompress_if_needed(rmf_path)
    logger.info(f"Passing RMF to XSPEC: {rmf_path}")

    arf_path = os.path.abspath(os.path.join(os.path.dirname(spectrum_file), raw_arf))
    arf_path = decompress_if_needed(arf_path)
    logger.info(f"Passing ARF to XSPEC: {arf_path}")

    # Prepare model-specific output directory
    model_out = os.path.join(output_dir, model_name)
    os.makedirs(model_out, exist_ok=True)

    # Run the BXA fitting
    fit_with_bxa(
        srcid=srcid,
        spectrum_file=spectrum_file,
        background_file=bkg_path,
        rmf_file=rmf_path,
        arf_file=arf_path,
        model_name=model_name,
        redshift=args_dict.get("redshift", 0.0),
        use_galabs=args_dict.get("use_galabs", False),
        use_tbabs_table=args_dict.get("use_tbabs_table", False),
        output_dir=model_out,
        log_file=log_file
    )


def xspec_model(model_name, redshift):
    """
    Create an XSPEC model based on the given model name and redshift.

    Parameters:
        model_name (str): Name of the XSPEC model.
        redshift (float): Redshift value.

    Returns:
        xspec.Model: The initialized XSPEC model.
    """
    if model_name == "powerlaw":
        return xspec.Model("powerlaw")
    elif model_name == "blackbody":
        return xspec.Model("bbody")
    # Add other models as needed
    else:
        raise ValueError(f"Unknown model name: {model_name}")


def fit_with_bxa(srcid, spectrum_file, background_file, rmf_file, arf_file,
                 model_name, redshift, use_galabs, use_tbabs_table,
                 output_dir, log_file):
    """
    Execute the BXA fitting workflow in XSPEC.

    Parameters:
        srcid (int): The source ID.
        spectrum_file (str): Path to the PHA spectrum file.
        background_file (str): Path to the background file.
        rmf_file (str): Path to the RMF file.
        arf_file (str): Path to the ARF file.
        model_name (str): Name of the model.
        redshift (float): Redshift for the model.
        use_galabs (bool): Whether to apply Galactic absorption.
        use_tbabs_table (bool): Whether to use TBabs table.
        output_dir (str): Output directory for fit results.
        log_file (str): Path to the log file.
    """
    try:
        logger.info(f"Output directory for SRCID {srcid}: {output_dir}")

        # Clear previous data and load spectrum
        xspec.AllData.clear()
        xspec.AllData(spectrum_file)
        spectrum = xspec.AllData(1)

        # Assign background, RMF, and ARF
        spectrum.background = background_file
        spectrum.response.rmf = rmf_file
        spectrum.response.arf = arf_file

        # Initialize model
        model = xspec_model(model_name, redshift)
        logger.info(f"Initialized model: {model.expression}")

        # Run BXA fit
        fit = bxa.Fit(model, output_dir)
        logger.info(f"Starting BXA fitting for {model_name}...")

        fit.run()
        result_file = os.path.join(output_dir, "fit_results.fits")
        fit.results(result_file)

        logger.info(f"Fitting complete. Results saved in {result_file}")
    except Exception as e:
        logger.error(f"SRCID {srcid}: BXA fitting failed. Error: {e}")

