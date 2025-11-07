import os
import datetime
from xspec import *
Fit.query = "no"  # Disable interactive prompts during fits
import bxa.xspec as bxa
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.table import Table, vstack
import glob
import logging

# python3 automated_fits.py 3067718060100029 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits dummy_output.txt --use_bxa --model_name=powerlaw --redshift=1.0 --overwrite=1 --export_results_fits --export_filename=fit_results.fits --bxa_output_dir=bxa_fit_results

#logger = logging.getLogger(__name__)
logger = logging.getLogger()   # root logger


    
def get_model_and_priors(model_name, redshift=0.0, flux_band=(0.5, 10.0)):
    """
    Construct an XSPEC model wrapped with cflux so that we fit for flux instead of norm.
    Parameters
    ----------
    model_name : str
        Name of the physical model ("powerlaw", "apec_single", "blackbody", "bremss").
    redshift : float
        Redshift for models that need it (e.g., zpowerlw).
    flux_band : tuple
        Energy band (Emin, Emax) in keV for cflux.
    """

    if model_name == "powerlaw":
        model = Model("phabs*cflux*zpowerlw")
        model.zpowerlw.Redshift = redshift
        model.zpowerlw.Redshift.frozen = True

        # set typical values
        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.zpowerlw.PhoIndex.values = "2.0,,1.0,1.0,3.0,3.0"

        # freeze the original norm (cflux will control the flux)
        model.zpowerlw.norm.frozen = True

    elif model_name == "apec_single":
        model = Model("phabs*cflux*apec")

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.apec.kT.values = "1.0,,0.1,0.1,10.0,10.0"
        model.apec.norm.frozen = True

    elif model_name == "blackbody":
        model = Model("phabs*cflux*bbody")

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.bbody.kT.values = "0.1,,0.01,0.01,2.0,2.0"
        model.bbody.norm.frozen = True

    elif model_name == "bremss":
        model = Model("phabs*cflux*bremss")

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.bremss.kT.values = "5.0,,0.1,0.1,20.0,20.0"
        model.bremss.norm.frozen = True

    else:
        raise ValueError(f"Unknown model: {model_name}")

    # Configure cflux energy range
    model.cflux.Emin = flux_band[0]
    model.cflux.Emax = flux_band[1]

    # Typical starting value for flux (erg/cm^2/s)
    model.cflux.Flux.values = "1e-12,,1e-15,1e-15,1e-9,1e-9"

    # Priors: always include flux instead of norm
    priors = [
        bxa.create_uniform_prior_for(model, model.phabs.nH),
    ]

    # Add temperature / index depending on model
    if model_name == "powerlaw":
        priors.append(bxa.create_uniform_prior_for(model, model.zpowerlw.PhoIndex))
    elif model_name == "apec_single":
        priors.append(bxa.create_uniform_prior_for(model, model.apec.kT))
    elif model_name == "blackbody":
        priors.append(bxa.create_uniform_prior_for(model, model.bbody.kT))
    elif model_name == "bremss":
        priors.append(bxa.create_uniform_prior_for(model, model.bremss.kT))

    # Finally, flux prior
    priors.append(bxa.create_loguniform_prior_for(model, model.cflux.Flux))

    return model, priors
    
    
# === NEW: unified background check that always returns flag=3 on any problem ===
def check_background_fit(spectrum_file, background_file, rmf_file, arf_file,
                         model_name, redshift, logger, srcid="unknown"):
    """
    Returns:
      pval (float) if background is OK (>= 0.01),
      None if background cannot be fit or p<0.01 or any exception occurred.
    """
    from xspec import AllData, AllModels, Fit, Spectrum
    import numpy as np
    try:
        # Clean local XSPEC state for the quick test
        AllData.clear()
        AllModels.clear()

        spectra = []
        for i in range(len(spectrum_files)):
            s = Spectrum(spectrum_files[i])
            s.background = background_files[i]
            s.response = rmf_files[i]
            s.response.arf = arf_files[i]
            spectra.append(s)

        AllData.ignore("**-0.3 10.0-**")


        # Build the same model you’ll use (so the set-up is consistent),
        # but *do not* run BXA here—just a quick XSPEC fit to get p-value.
        # Re-use your existing helper:
        try:
            _model, _priors = get_model_and_priors(model_name, redshift)
        except Exception as e:
            logger.warning(f"   Background check: could not build model for {srcid}: {e}")
            return None

        # Perform a quick local fit
        Fit.perform()

        # Convert the fit test statistic to a p-value (approximate; same logic as before)
        try:
            from scipy.stats import chi2
            bg_stat = Fit.testStatistic()
            dof = Fit.dof
            pval = 1.0 - chi2.cdf(bg_stat, dof)
        except Exception as e:
            logger.warning(f"   Background check: cannot compute p-value for {srcid}: {e}")
            return None

        if (pval is None) or (not np.isfinite(pval)) or (pval < 0.01):
            logger.warning(f"   Background check: p={pval:.4f} (<0.01) → NOT OK for {srcid}")
            return None

        logger.info(f"   Background check: p={pval:.4f} → OK for {srcid}")
        return pval

    except Exception as e:
        logger.warning(f"   Background check failed for {srcid}: {e}")
        return None



def fit_spectrum_bxa(spectrum_file, background_file, rmf_file, arf_file,
                     redshift=0.0, model_name="powerlaw",
                     output_base="bxa_fit_results", srcid="unknown", log_file="fit_spectrum_bxa.log"):

    logger.info('\n')
    logger.info(f'Starting BXA fit on spectrum {spectrum_file}')
    dirname = os.path.dirname(spectrum_file)
    os.chdir(dirname)
    logger.info(f'   Changing focus to {dirname}')

    # === NEW: unified background gate (flag=3 on any problem) ===
    pval = check_background_fit(
        spectrum_file, background_file, rmf_file, arf_file,
        model_name, redshift, logger, srcid=str(srcid)
    )
    if pval is None:
        # Any failure (no bkg, could not fit, or p<0.01) → flag 3 and bail out
        return {"flag": 3}

    # === From here on, your original working fitting code ===
    AllData.clear()
    AllModels.clear()

    Fit.statMethod = "cstat"
    Plot.device = "/null"

    spectra = []
    for i in range(len(spectrum_files)):
        s = Spectrum(spectrum_files[i])
        s.background = background_files[i]
        s.response = rmf_files[i]
        s.response.arf = arf_files[i]
        spectra.append(s)

    AllData.ignore("**-0.3 10.0-**")


    # Model + priors
    model, priors_list = get_model_and_priors(model_name, redshift)

    # Output dir
    timestamp = datetime.datetime.now().strftime("%d%m%Y_%H%M")
    model_dirname = f"{model_name}_{timestamp}"
    output_dir = os.path.abspath(os.path.join(output_base, str(srcid), model_dirname))
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f'   Setting the output directory for the fits {output_dir}')

    # Run BXA
    solver = bxa.BXASolver(transformations=priors_list,
                           outputfiles_basename=os.path.join(output_dir))
    solver.run(resume=False)

    # Read chain, make plots, return summaries
    chain_file = os.path.join(output_dir, "chain.fits")
    if os.path.exists(chain_file):
        with fits.open(chain_file) as hdul:
            samples = hdul[1].data
            samples_array = np.column_stack([samples[name] for name in samples.names])
            labels = solver.paramnames
            stds = np.std(samples_array, axis=0)
            valid_cols = stds > 0
            filtered_samples = samples_array[:, valid_cols]
            filtered_labels = [label for i, label in enumerate(labels) if valid_cols[i]]

            if filtered_samples.shape[1] > 0:
                fig = corner.corner(filtered_samples, labels=filtered_labels, show_titles=True, title_fmt=".3e")
                corner_path = os.path.join(output_dir, "corner.png")
                fig.savefig(corner_path)
                logger.info(f'   Saved corner plot to file {corner_path} ')

        posterior_median = np.median(samples_array, axis=0)
        posterior_p16    = np.percentile(samples_array, 16, axis=0)
        posterior_p84    = np.percentile(samples_array, 84, axis=0)
        return {
            "parameter_names": labels,
            "posterior_median": posterior_median,
            "posterior_p16": posterior_p16,
            "posterior_p84": posterior_p84,
            "output_dir": output_dir,
            "flag": 0
        }

    else:
        logger.error(f'   Chain file {chain_file} not found after BXA run ')
        return {"flag": 4}


def export_bxa_results_to_fits(srcid, output_base="bxa_fit_results", fits_filename="fit_results.fits", log_file="fit_spectrum_bxa.log", global_results=False):
    # directory containing the fit results
    src_dir = os.path.join(output_base, str(srcid))
    #
    # generating filename to save the summary fit results to
    if (global_results):
        # adding them to a file including all previous results
        fits_path = os.path.join(output_base, fits_filename)
    else:
        # adding/writing them to a file in the directory with the fit results
        fits_path = os.path.join(src_dir, fits_filename)
    #
    # print('\n\n Inside export...')
    # print(f'    output_base=({output_base})')
    # print(f'    fits_path=({fits_path})')
    # print(f'    src_dir=({src_dir})')
    
    logger.info('\n')
    logger.info(f'Exporting BXA fit results to FITS file {fits_path}')

    os.makedirs(output_base, exist_ok=True)

    # Get all models in SRCID directory
    #    sorted so that, for each model, the last in the list is the latest fit
    model_dirs = sorted( [d for d in os.listdir(src_dir) if os.path.isdir(os.path.join(src_dir, d))] )

    short_map = {"powerlaw": "PL", "blackbody": "BB", "bremss": "BR", "apec_single": "AP"}

    model_data = {}
    for mdir in model_dirs:
        for long_name, short_name in short_map.items():
            if mdir.startswith(long_name):
                model_chain = os.path.join(src_dir, mdir, "chain.fits")
                if os.path.exists(model_chain):
                    with fits.open(model_chain) as hdul:
                        samples = hdul[1].data
                        param_names = hdul[1].columns.names
                        medians = [np.percentile(samples[p], 50) for p in param_names]
                        p16 = [np.percentile(samples[p], 16) for p in param_names]
                        p84 = [np.percentile(samples[p], 84) for p in param_names]
                        model_data[short_name] = {
                            "names": param_names,
                            "medians": medians,
                            "p16": p16,
                            "p84": p84
                        }

    # Create or update table without Models column
    if os.path.exists(fits_path):
        table = Table.read(fits_path)
        src_mask = table["SRCID"] == srcid
        if np.any(src_mask):
            # table exists and data for srcid are already there, updating them
            idx = np.where(src_mask)[0][0]
            for model_short, pdata in model_data.items():
                for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                    col_med = f"{model_short}_{pname}_median"
                    col_p16 = f"{model_short}_{pname}_p16"
                    col_p84 = f"{model_short}_{pname}_p84"
                    for col, val in zip([col_med, col_p16, col_p84], [median, p16, p84]):
                        if col not in table.colnames:
                            table[col] = np.full(len(table), np.nan)
                        table[col][idx] = val
        else:
            # table exists, but data for srcid are not there yet, adding them
            new_row = {"SRCID": srcid}
            for model_short, pdata in model_data.items():
                for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                    new_row[f"{model_short}_{pname}_median"] = median
                    new_row[f"{model_short}_{pname}_p16"] = p16
                    new_row[f"{model_short}_{pname}_p84"] = p84
            for col in new_row.keys():
                if col not in table.colnames:
                    table[col] = np.full(len(table), np.nan)
            table.add_row(new_row)
        table.write(fits_path, overwrite=True)
    else:
        # table does not exist, creating it with the data for srcid
        row_data = {"SRCID": [srcid]}
        for model_short, pdata in model_data.items():
            for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                row_data[f"{model_short}_{pname}_median"] = [median]
                row_data[f"{model_short}_{pname}_p16"] = [p16]
                row_data[f"{model_short}_{pname}_p84"] = [p84]
        Table(row_data).write(fits_path, overwrite=True)

