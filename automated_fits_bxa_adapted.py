import os
from xspec import *
import bxa.xspec as bxa
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import corner


def fit_spectrum_bxa(spectrum_file, background_file, rmf_file, arf_file, redshift=0.0):
    # Clear XSPEC state
    AllData.clear()
    AllModels.clear()

    # Set fit statistic
    Fit.statMethod = "cstat"

    # Load spectrum
    s = Spectrum(spectrum_file)
    s.background = background_file
    s.response = rmf_file
    s.response.arf = arf_file

    # Ignore bad channels
    s.ignore("**-510")

    # Define model: phabs*zpowerlw
    model = Model("phabs*zpowerlw")

    # Set initial values and limits
    model.phabs.nH = 1.0
    model.phabs.nH.frozen = False
    model.phabs.nH.values = (1.0, 0.01, 0.001, 10.0, 10.0, 10.0)

    model.zpowerlw.PhoIndex = 1.9
    model.zpowerlw.PhoIndex.frozen = False
    model.zpowerlw.PhoIndex.values = (1.9, 0.1, 1.0, 3.0, 3.0, 3.0)

    model.zpowerlw.Redshift = redshift
    model.zpowerlw.Redshift.frozen = True

    model.zpowerlw.norm = 1e-4
    model.zpowerlw.norm.frozen = False
    model.zpowerlw.norm.values = (1e-4, 1e-5, 1e-6, 1e-2, 1e-2, 1e-2)

    # Setup BXA priors
    priors = [
        bxa.create_uniform_prior_for(model.phabs.nH),
        bxa.create_uniform_prior_for(model.zpowerlw.PhoIndex),
        bxa.create_jeffreys_prior_for(model.zpowerlw.norm)
    ]

    # Run nested sampling
    solver = bxa.BXASolver(priors)
    results = solver.run(resume=False)

    # Plot corner
    samples_array = solver.posterior_samples
    labels = solver.paramnames

    fig = corner.corner(samples_array, labels=labels, show_titles=True, title_fmt=".3e")
    fig.savefig("bxa_fit_result/corner_plot.png")

    return {
        "parameter_names": labels,
        "posterior_median": np.median(samples_array, axis=0)
    }

