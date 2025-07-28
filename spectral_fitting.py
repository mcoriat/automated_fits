#import argparse
import os
#import sys
#import glob
#import numpy as np
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
    - good_spectra (list): The list of good spectra, made of tuples containing the full path and name for the good pn spectra, total source counts, total background counts, total net counts, exposure time, a flag, the signal-to-noise ratio, and a string with the name of the instrument
    - output_dir (str) full path to the output directory corresponding to srcid

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
        spectra (list): The list of good spectra, made of tuples containing the full path and name for the good pn spectra, total source counts, total background counts, total net counts, exposure time, a flag, the signal-to-noise ratio, and a string with the name of the instrument

        args (Namespace): The command line arguments.
        model_name (str): The name of the model.
        log_file (str): The path to the log file.
        output_dir (str) full path to the output directory corresponding to srcid

    Returns:
        None
    """
    
    message='\n\n Starting'
    logger.info(message)
    
    best_spectrum = max(spectra, key=lambda x: x[6])  # Choosing based on highest SNR
    
    #
    spectrum_file=best_spectrum[0]
    try:
        hdul=fits.open(spectrum_file)
        background_file=hdul[1].header['BACKFILE']
        hdul.close()
        del hdul
    except:
        message = f" Could not open spectrum {spectrum_file} for SRCID {srcid}"
        logger.error(message)
        return
    
    
    # Log the fitting process and perform the BXA fitting
    message = f" Using BXA to fit {model_name} model for SRCID {srcid} to spectrum {spectrum_file}"
    logger.info(message)
    output_dir = os.path.join(output_dir, model_name)
    fit_with_bxa(srcid, spectrum_file, background_file, model_name, args.redshift, args.use_galabs, args.use_tbabs_table, output_dir, log_file)


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
    elif model_name == "powerlaw_blackbody":
        model = xspec.Model("powerlaw + bbody")
    elif model_name == "powerlaw_blackbody_const":
        model = xspec.Model("powerlaw + bbody + constant")
    elif model_name == "zpowlaw":
        model = xspec.Model(f"zpowerlw,zpowerlw.redshift={redshift}")
    elif model_name == "double_zpowlaw":
        model = xspec.Model(f"zpowerlw + zpowerlw,zpowerlw.redshift={redshift}")
    else:
        raise ValueError(f"Unknown model name: {model_name}")
    return model


# Function to perform the BXA fitting for a given SRCID and spectrum
def fit_with_bxa(srcid, spectrum_file, background_file, model_name, redshift, use_galabs, use_tbabs_table, output_dir, log_file):
    """
    Fits a spectrum using the BXA (Bayesian X-ray Analysis) method.

    Args:
        srcid (int): The source ID.
        spectrum_file (str): The path to the spectrum file.
        background_file (str): The path to the background file.
        model_name (str): The name of the model.
        redshift (float): The redshift value.
        use_galabs (bool): Whether to use galactic absorption.
        use_tbabs_table (bool): Whether to use the tbabs table.
        output_dir (str): The output directory.
        log_file (str): The path to the log file.

    Raises:
        Exception: If any error occurs during the fitting process.

    Returns:
        None
    """
    
    message='\n\n Starting'
    logger.info(message)
    
    try:
        # changing focus to directory where the spectra are
        dirname=os.path.dirname(spectrum_file)
        os.chdir(dirname)
        
        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Load the spectrum and background data into XSPEC
        xspec.AllData.clear()
        xspec.AllData(spectrum_file)
        spectrum = xspec.AllData(1)
        spectrum.background = background_file

        # Define the model based on the specified model name
        model = xspec_model(model_name, redshift)

        # Apply additional model options
        if use_galabs:
            model.setPars(galabs=True)
        if use_tbabs_table:
            model.setPars(tbabs=True)

        # Configure and run the BXA fitting process
        ###### bxa.Fit does not exist ########
        fit = bxa.Fit(model, output_dir)
        fit.run()

        # Save the fitting results
        fit.results(output_dir + "/fit_results.fits")
    except Exception as e:
        # Log any errors encountered during the fitting process
        message = f" SRCID {srcid}: BXA fitting failed for spectrum {spectrum_file} with model {model_name}. Error: {str(e)}"
        logger.error(message)

