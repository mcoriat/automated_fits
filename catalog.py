import json
import logging
from itertools import product

import numpy as np

import stats


def _get_parameters(model_name):
    if model_name == "torus":
        parameters = ["lognh", "gamma", "theta_inc"]
    elif model_name == "powerlaw":
        parameters = ["lognh", "gamma"]
    else:
        raise ValueError("Unknown model!")
        
    return parameters


def _add_columns(table, parameters, ebands):
    table["logz"] = 0.0
    table["cstat"] = 0.0
    table["dof"] = 0
    table["ppp"] = 0.0
    table["ks"] = 0.0
    table["ks_pvalue"] = 0.0
    table["bkg_cstat"] = 0.0
    table["bkg_dof"] = 0
    table["bkg_ks"] = 0.0
    table["bkg_ks_pvalue"] = 0.0

    for p in parameters:
        table[p] = 0.0
        table[f"{p}_errlo"] = 0.0
        table[f"{p}_errhi"] = 0.0
    
    for c, e in product(["flux", "lumin"], ebands):
        table[f"{c}_{e}"] = 0.0
        table[f"{c}_{e}_errlo"] = 0.0
        table[f"{c}_{e}_errhi"] = 0.0

    return table


def _add_src_stats(src, info_path):
    goodness_path = info_path.joinpath("goodness.json")
    goodness = json.load(goodness_path.open("r"))

    src["cstat"] = goodness["src"]["cstat"]
    src["dof"] = goodness["src"]["dof"]
    src["ppp"] = goodness["src"]["ppp"]
    src["ks"] = goodness["src"]["KS"]["stat"]
    src["ks_pvalue"] = goodness["src"]["KS"]["p-value"]

    src["bkg_cstat"] = goodness["bkg"]["cstat"]
    src["bkg_dof"] = goodness["bkg"]["dof"]
    src["bkg_ks"] = goodness["bkg"]["KS"]["stat"]
    src["bkg_ks_pvalue"] = goodness["bkg"]["KS"]["p-value"]
    

def _add_src_parameters(src, parameters, info_path):
    results_path = info_path.joinpath("results.json")
    results = json.load(results_path.open("r"))

    params_median = results["posterior"]["median"]
    params_lo = results["posterior"]["errlo"]
    params_hi = results["posterior"]["errup"]

    src["logz"] = results["logz"]
    for i, p in enumerate(parameters):
        src[p] = params_median[i]
        src[f"{p}_errlo"] = params_median[i] - params_lo[i]
        src[f"{p}_errhi"] = params_hi[i] - params_median[i]


def _add_src_flux_lumin(src, ebands, extra_path):
    for c in ["flux", "lumin"]:
        dist_path = extra_path.joinpath(f"src_{c}_dist.dat")
        dist = np.loadtxt(dist_path)
        perc = np.percentile(dist, [5, 50, 95], axis=0)

        for i, b in enumerate(ebands):
            src[f"{c}_{b}"] = perc[1, i]
            src[f"{c}_{b}_errlo"] = perc[1, i] - perc[0, i]
            src[f"{c}_{b}_errhi"] = perc[2, i] - perc[1, i]

    
def _add_src_flux_lumin_mode(src, ebands, extra_path):
    for c in ["flux", "lumin"]:
        dist_path = extra_path.joinpath(f"src_{c}_dist.dat")
        dist = np.loadtxt(dist_path)
        mode = stats.mode(dist, axis=0)
        hdpi = stats.hdp_intervals(dist)

        for i, b in enumerate(ebands):
            src[f"{c}_{b}"] = mode[i]
            src[f"{c}_{b}_errlo"] = mode[i] - hdpi[0, i]
            src[f"{c}_{b}_errhi"] = hdpi[1, i] - mode[i]


def collect_results(table, model_name, data_path, ebands=("soft", "hard")):
    parameters = _get_parameters(model_name)
    results_table = _add_columns(table, parameters, ebands)

    for src in results_table:
        info_path = data_path.joinpath(src["srcid"], f"bxafit_{model_name}", "info")

        try:
            _add_src_stats(src, info_path)

        except FileNotFoundError:
            logging.warn(f"No data for {src['srcid']}")
            continue

        _add_src_parameters(src, parameters, info_path)
    
        extra_path = data_path.joinpath(src["srcid"], f"bxafit_{model_name}", "extra")
        _add_src_flux_lumin_mode(src, ebands, extra_path)

    return table

