# Automated X-ray Spectral Fitting with BXA

This Python script performs automated spectral fitting using BXA (Bayesian X-ray Analysis) for data from the stacked version of the 4XMM-DR11 catalog. It processes each SRCID and its associated OBSIDs, checks the validity of spectra, and performs the fitting for specified models, logging the results for each source.

## Features

- Reads and processes the stacked 4XMM-DR11 catalog.
- Identifies valid spectra for each source and performs BXA fitting.
- Supports various models (e.g., power law, blackbody, apec, and more).
- Generates detailed logs for each SRCID, stored in the output directory.
- Supports fitting multiple models in a single run.

## Dependencies

To run this script, you will need the following Python packages:

- `astropy` - for handling FITS files
- `numpy` - for numerical operations
- `bxa` - for Bayesian X-ray Analysis integration
- `xspec` - for X-ray spectral fitting

To install the required dependencies, you can use:

```bash
pip install astropy numpy bxa xspec
