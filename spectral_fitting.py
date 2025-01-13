import os
import bxa.xspec as bxa
import xspec
import logging
from astropy.io import fits

logger = logging.getLogger(__name__)


# Function to perform the BXA fitting for a given SRCID based on selected model
def perform_spectrum_fitting(args, srcid, log_file, good_spectra, output_dir):
    """
    Perform spectrum fitting for the given source ID and good spectra.

    Parameters:
    - args: The command line arguments.
    - srcid (int): The source ID.
    - log_file (str): The log file to write the fitting results.
    - good_spectra (list): List of dictionaries with good spectra details, including:
        - spectrum_file
        - sp_counts
        - bg_counts
        - sp_netcts
        - sp_exp
        - flag
        - snr
        - instrument
    - output_dir (str): Path to the output directory corresponding to srcid.

    Returns:
    None
    """
    if args.fit_pl:
        fit_spectrum(srcid, good_spectra, args, "powerlaw", log_file, output_dir)
    if args.fit_bb:
        fit_spectrum(srcid, good_spectra, args, "blackbody", log_file, output_dir)
    if args.fit_apec_single:
        fit_spectrum(srcid, good_spectra, args, "apec_single", log_file, output_dir)
    if args.fit_apec_apec:
        fit_spectrum(srcid, good_spectra, args, "apec_apec", log_file, output_dir)
    if args.fit_apec_apec_const:
        fit_spectrum(srcid, good_spectra, args, "apec_apec_const", log_file, output_dir)
    if args.fit_bremss:
        fit_spectrum(srcid, good_spectra, args, "bremss", log_file, output_dir)
    if args.fit_bbpl:
        fit_spectrum(srcid, good_spectra, args, "powerlaw_blackbody", log_file, output_dir)
    if args.fit_bbpl_const:
        fit_spectrum(srcid, good_spectra, args, "powerlaw_blackbody_const", log_file, output_dir)
    if args.fit_bbpl_const2:
        fit_spectrum(srcid, good_spectra, args, "powerlaw_blackbody_const2", log_file, output_dir)
    if args.fit_zpl:
        fit_spectrum(srcid, good_spectra, args, "zpowlaw", log_file, output_dir)
    if args.fit_zplpl:
        fit_spectrum(srcid, good_spectra, args, "double_zpowlaw", log_file, output_dir)


# Function to select the best spectrum based on SNR and perform fitting
def fit_spectrum(srcid, spectra, args, model_name, log_file, output_dir):
    """
    Fits a spectrum using the BXA method.

    Parameters:
        srcid (int): The source ID.
        spectra (list): List of dictionaries containing spectrum details:
                        - spectrum_file
                        - sp_counts
                        - bg_counts
                        - sp_netcts
                        - sp_exp
                        - flag
                        - snr
                        - instrument
        args (Namespace): Command line arguments.
        model_name (str): The name of the model.
        log_file (str): Path to the log file.
        output_dir (str): Path to the output directory.

    Returns:
        None
    """
    message = '\n\n Starting spectrum fitting...'
    logger.info(message)

    # Select the best spectrum based on the highest SNR
    best_spectrum = max(spectra, key=lambda x: x['snr'])
    
    spectrum_file = best_spectrum['spectrum_file']
    try:
        hdul = fits.open(spectrum_file)
        background_file = hdul[1].header['BACKFILE']
        hdul.close()
    except Exception as e:
        message = f"Could not open spectrum {spectrum_file} for SRCID {srcid}. Error: {e}"
        logger.error(message)
        return
    
    # Log the fitting process and perform the BXA fitting
    message = f"Using BXA to fit {model_name} model for SRCID {srcid} to spectrum {spectrum_file}"
    logger.info(message)
    output_dir = os.path.join(output_dir, model_name)
    
    fit_with_bxa(
        srcid=srcid,
        spectrum_file=spectrum_file,
        background_file=background_file,
        model_name=model_name,
        redshift=args.redshift,
        use_galabs=args.use_galabs,
        use_tbabs_table=args.use_tbabs_table,
        output_dir=output_dir,
        log_file=log_file
    )


# Function to define the XSPEC model based on the model name
def xspec_model(model_name, redshift):
    """
    Create an XSPEC model based on the given model name and redshift.

    Parameters:
    - model_name (str): The name of the XSPEC model to create.
    - redshift (float): The redshift value to use for the model.

    Returns:
    - model (xspec.Model): The created XSPEC model.

    Raises:
    - ValueError: If the given model name is unknown.
    """
    if model_name == "powerlaw":
        model = xspec.Model("powerlaw")
    elif model_name == "blackbody":
        model = xspec.Model("bbody")
    elif model_name == "apec_single":
        model = xspec.Model("apec")
    elif model_name == "apec_apec":
        model = xspec.Model("apec + apec")
    elif model_name == "apec_apec_const":
        model = xspec.Model("apec + apec + constant")
    elif model_name == "bremss":
        model = xspec.Model("bremss")
    else:
        raise ValueError(f"Unknown model name: {model_name}")
    return model


# Function to perform the BXA fitting
def fit_with_bxa(srcid, spectrum_file, background_file, model_name, redshift, use_galabs, use_tbabs_table, output_dir, log_file):
    """
    Fits a spectrum using the BXA (Bayesian X-ray Analysis) method.

    Args:
        srcid (int): The source ID.
        spectrum_file (str): Path to the spectrum file.
        background_file (str): Path to the background file.
        model_name (str): Name of the model.
        redshift (float): The redshift value.
        use_galabs (bool): Whether to use galactic absorption.
        use_tbabs_table (bool): Whether to use the tbabs table.
        output_dir (str): Output directory.
        log_file (str): Path to the log file.

    Returns:
        None
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        xspec.AllData.clear()
        xspec.AllData(spectrum_file)
        spectrum = xspec.AllData(1)
        spectrum.background = background_file

        model = xspec_model(model_name, redshift)

        fit = bxa.Fit(model, output_dir)
        fit.run()
        fit.results(os.path.join(output_dir, "fit_results.fits"))
    except Exception as e:
        message = f"SRCID {srcid}: BXA fitting failed. Error: {str(e)}"
        logger.error(message)

