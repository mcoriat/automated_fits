import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["HEADASOUTPUT"] = "/dev/null"

import json
import logging
import warnings
from copy import deepcopy
from itertools import count
from pathlib import Path

import bxa.sherpa as bxa
import numpy as np
import sherpa.astro.ui as shp
from scipy.stats import ks_2samp

sys.path.append(r'/dataAGN/hstiele/python')

import logs
import models
import plots
import stats
import utils
import fit
from fluxes import calc_flux_dist, calc_lumin_dist
from priors import BXAPrior
from source import Source
from xmm_backgrounds import get_pn_bkg_model_cached, get_mos_bkg_model_cached

fit.run(1, powerlaw, src1, false, 1)
