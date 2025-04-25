
import os
import numpy as np
import logging

logger = logging.getLogger(__name__)

def get_bkg(path):
    class Dummy:
        def __init__(self, path):
            self.name = os.path.splitext(os.path.basename(path))[0]
    return Dummy(path)

def dummy_fit_background_model(path, instrument, galabs, fit):
    logger.info(f"Pretending to fit {instrument.upper()} background model for file {path} (fit={fit})")

    bkg_model = {
        "scaling_factor": 1.0,
        "unit_response": None,
        "background_response": None,
        "galactic_absorption": galabs,
        "fit": fit,
        "fit_statistic": 42.0 if fit else -1.0,
        "fit_dof": 10 if fit else -1,
        "model": f"{instrument}_bkg_model_placeholder",
    }

    return bkg_model

def get_pn_bkg_model(path, galabs, fit=False, fix_response=False):
    return dummy_fit_background_model(path, "pn", galabs, fit)

def get_mos_bkg_model(path, galabs, fit=False, fix_response=False):
    return dummy_fit_background_model(path, "mos", galabs, fit)

def get_bkg_model_cached(path, galabs, instrument="mos"):
    filename = get_bkg(path).name + ".bkgpars"
    getter = get_mos_bkg_model if instrument == "mos" else get_pn_bkg_model
    if os.path.exists(filename):
        model = getter(path, galabs, fit=False)
        try:
            params = np.loadtxt(filename)
            model["fit_statistic"], model["fit_dof"] = params
        except Exception as e:
            logger.warning(f"Could not read parameters from {filename}: {e}")
    else:
        model = getter(path, galabs, fit=True)
        try:
            np.savetxt(filename, [model["fit_statistic"], model["fit_dof"]])
        except Exception as e:
            logger.warning(f"Could not write parameters to {filename}: {e}")
    return model

def get_pn_bkg_model_cached(path, galabs):
    return get_bkg_model_cached(path, galabs, instrument="pn")

def get_mos_bkg_model_cached(path, galabs):
    return get_bkg_model_cached(path, galabs, instrument="mos")
