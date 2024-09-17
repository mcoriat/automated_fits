import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["HEADASOUTPUT"] = "/dev/null"

import json
import logging
import warnings
from copy import deepcopy
from itertools import count
from pathlib import Path
from astropy.coordinates import SkyCoord
from scipy.stats import chi2

import bxa.sherpa as bxa
import numpy as np
import sherpa.astro.ui as shp
from scipy.stats import ks_2samp

import logs
import models
import plots
import stats
import utils
from fluxes import calc_flux_dist, calc_lumin_dist
from priors import BXAPrior
from source import Source
from xmm_backgrounds import get_pn_bkg_model_cached, get_mos_bkg_model_cached
from bxa.sherpa.background.pca import auto_background

data_path = Path(".")


def _set_output_path(srcid, model_name=None):
    folder_name = f"bxafit_{model_name}" if model_name else "bxafit"
    return data_path.joinpath(srcid, folder_name)


def _check_fit_exists(output_path):
    results_json_path = output_path.joinpath("info", "results.json")
    return results_json_path.exists()


def _set_sherpa_env(stat="cstat"):
    shp.clean()
    shp.set_stat(stat)
    shp.set_xsabund('wilm')
    shp.set_xsxsect('vern')

    plots.matplotlib_settings()

def _set_sherpa_env2(stat="chi2gehrels"):
    shp.clean()
    shp.set_stat(stat)
    shp.set_xsabund('wilm')
    shp.set_xsxsect('vern')

    plots.matplotlib_settings()


def _set_obs_prefix(coadd):
    if coadd:
        return "coadd"
    else:
        return "P*"


def _load_data(
    srcid,
    coadd=False,
    globstr=''   # 2023-02-07 AV: allow for globbing for different det_there + det_use combinations
):

    source_path = data_path
    #source_path = data_path.joinpath(obsid)
    #obs_prefix = _set_obs_prefix(coadd)

    id = 1
    #for obsid_path in source_path.glob(obs_prefix):
    for spec_path in source_path.glob(f"{globstr}*SRSPEC{srcid}.FTZ"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        shp.load_pha(id, str(spec_path))
        shp.ungroup(id)
        shp.ungroup(id, bkg_id=1)
        #shp.group_counts(id, 1)
        #shp.group_counts(id, 1, bkg_id=1)

        #if coadd:
        #shp.subtract(id)

        shp.ignore_id(id, lo=None, hi=0.2)
        shp.ignore_id(id, lo=12.0, hi=None)
        id += 1

    ids = shp.list_data_ids()

    if not ids:
        raise Exception("No data!!!")

    return ids

def _load_data2(srcid, coadd=False):

    source_path = data_path
    #source_path = data_path.joinpath(obsid)
    #obs_prefix = _set_obs_prefix(coadd)

    id = 1
    #for obsid_path in source_path.glob(obs_prefix):
    for spec_path in source_path.glob(f"*SRSPEC{srcid}.FTZ"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        shp.load_pha(id, str(spec_path))
        detector = utils.get_detector(id)

        if detector == "EPN" or detector == "EPNA":
            shp.group_bins(id, num=53)
        else:
            shp.group_bins(id, num=31)

        #shp.group_bins(id, num=30, bkg_id=1)
        #shp.group_counts(id, 1)
        #shp.group_counts(id, 1, bkg_id=1)

        #if coadd:
        #shp.subtract(id)

        shp.ignore_id(id, lo=None, hi=0.2)
        shp.ignore_id(id, lo=12.0, hi=None)
        id += 1

    ids = shp.list_data_ids()

    if not ids:
        raise Exception("No data!!!")

    return ids

def _load_data3(srcid, coadd=False):

    source_path = data_path
    #source_path = data_path.joinpath(obsid)
    #obs_prefix = _set_obs_prefix(coadd)

    id = 1
    #for obsid_path in source_path.glob(obs_prefix):
    for spec_path in source_path.glob(f"*SRSPEC{srcid}.FTZ"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        shp.load_pha(id, str(spec_path))
        shp.group_counts(id, 20)
        shp.group_counts(id, 20, bkg_id=1)
        #shp.group_bins(id, num=30, bkg_id=1)
        #shp.group_counts(id, 1)
        #shp.group_counts(id, 1, bkg_id=1)

        #if coadd:
        #shp.subtract(id)

        shp.ignore_id(id, lo=None, hi=0.2)
        shp.ignore_id(id, lo=12.0, hi=None)
        id += 1

    ids = shp.list_data_ids()

    if not ids:
        raise Exception("No data!!!")

    return ids

def _load_data4(srcid, coadd=False):

    source_path = data_path
    #source_path = data_path.joinpath(obsid)
    #obs_prefix = _set_obs_prefix(coadd)

    id = 1
    #for obsid_path in source_path.glob(obs_prefix):
    for spec_path in source_path.glob(f"*SRSPEC{srcid}.FTZ"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        shp.load_pha(id, str(spec_path))
        shp.group_counts(id, 10)
        shp.group_counts(id, 10, bkg_id=1)
        #shp.group_bins(id, num=30, bkg_id=1)
        #shp.group_counts(id, 1)
        #shp.group_counts(id, 1, bkg_id=1)

        #if coadd:
        #shp.subtract(id)

        shp.ignore_id(id, lo=None, hi=0.2)
        shp.ignore_id(id, lo=12.0, hi=None)
        id += 1

    ids = shp.list_data_ids()

    if not ids:
        raise Exception("No data!!!")

    return ids

def _load_data5(srcid, coadd=False):

    source_path = data_path
    #source_path = data_path.joinpath(obsid)
    #obs_prefix = _set_obs_prefix(coadd)

    id = 1
    #for obsid_path in source_path.glob(obs_prefix):
    for spec_path in source_path.glob(f"*SRSPEC{srcid}.FTZ"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        shp.load_pha(id, str(spec_path))
        shp.group_counts(id, 20)
        shp.group_counts(id, 20, bkg_id=1)
        #shp.group_bins(id, num=30, bkg_id=1)
        #shp.group_counts(id, 1)
        #shp.group_counts(id, 1, bkg_id=1)

        #if coadd:
        #shp.subtract(id)

        shp.ignore_id(id, lo=None, hi=0.5)
        shp.ignore_id(id, lo=12.0, hi=None)
        id += 1

    ids = shp.list_data_ids()

    if not ids:
        raise Exception("No data!!!")

    return ids


def _get_spec_info(ids, output_path=None):
    counts = []
    expo = []
    detect = []
    rates = []
    bkg_cts = []
    net_cts = []
    bsc_src = []
    bsc_bkg = []
    bscf = []

    for id in ids:
        detect.append(utils.get_detector(id))
        counts.append(shp.calc_data_sum(id=id))
        bkg_cts.append(shp.calc_data_sum(id=id, bkg_id=1))
        bsc_src.append(shp.get_backscal(id=id))
        bsc_bkg.append(shp.get_backscal(id=id, bkg_id=1))
        bscf.append(shp.get_bkg_scale(id=id))
        shp.subtract(id)
        net_cts.append(shp.calc_data_sum(id=id))
        expo.append(shp.get_exposure(id=id))
        rates.append(shp.calc_data_sum(id=id)/shp.get_exposure(id=id))

    a = []
    for i, id in enumerate(ids):
        a.append([detect[i], counts[i], bkg_cts[i], net_cts[i], expo[i], rates[i], bsc_src[i], bsc_bkg[i], bscf[i]])

    #if output_path:
    #    file_path = output_path / "extra"
    #    if not file_path.exists():
    #        file_path.mkdir()

    #    filename = Path(output_path, "extra", f"spec_info.txt")
    #    with open(filename, 'w', 8192) as testfile:
    #        for row in a:
    #            testfile.write(' '.join([str(j) for j in row]) + '\n')

    return a


def _background_models(ids, galabs):
    bkgmodels = []
    for id in ids:
        detector = utils.get_detector(id)

        if detector == "EPN" or detector == "EPNA":
            get_epic_bkg_model = get_pn_bkg_model_cached
        else:
            get_epic_bkg_model = get_mos_bkg_model_cached

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bkgmodels.append(get_epic_bkg_model(id, galabs))

    goodness_stats_bkg = _goodness_fixed(ids, bkg=True)
    _check_goodness_bkgfit(goodness_stats_bkg)

    return bkgmodels, goodness_stats_bkg

def _background_models2(ids, galabs):
    bkgmodels = []
    for id in ids:
        detector = utils.get_detector(id)

        if detector == "EPN" or detector == "EPNA":
            get_epic_bkg_model = get_pn_bkg_model_cached
        else:
            get_epic_bkg_model = get_mos_bkg_model_cached

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bkgmodels.append(get_epic_bkg_model(id, galabs))

    goodness_stats_bkg = _goodness_fixed(ids, bkg=True)
    _check_goodness_bkgfit(goodness_stats_bkg)

    return bkgmodels, goodness_stats_bkg

def _background_models_pca(src, ids, srcmodel, model_name):
    bkgmodels = []
    for id in ids:
        shp.set_model(id, srcmodel[id-1])
        convmodel = shp.get_model(id)
        bkg_model = auto_background(id)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bkgmodels.append(bkg_model * shp.get_bkg_scale(id))

        b = shp.get_bkg_fit_plot(id)
        np.savetxt(src + '_' + model_name + '_pca_bkg_'+str(id)+'.txt.gz', np.transpose([b.dataplot.x, b.dataplot.y, b.modelplot.x, b.modelplot.y]))
        m = shp.get_fit_plot(id)
        np.savetxt(src + '_' + model_name + '_pca_nosrc_'+str(id)+'.txt.gz', np.transpose([m.dataplot.x, m.dataplot.y, m.modelplot.x, m.modelplot.y]))


    goodness_stats_bkg = _goodness(ids, bkg=True)
    _check_goodness_bkgfit(goodness_stats_bkg)

    return bkgmodels, goodness_stats_bkg


def _set_full_model_src_bkg(ids, galabs, srcmodel, bkgmodels, use_galabs=False):
    for i, id in enumerate(ids):
        rsp = shp.get_response(id)
        if not use_galabs:
            shp.set_full_model(id, rsp(srcmodel[i]) + bkgmodels[i])  # bkg_scale included in bkgmodel
        else:
            shp.set_full_model(id, rsp(galabs * srcmodel[i]) + bkgmodels[i])


def _set_full_model_bkg(ids, bkgmodels):
    for i, id in enumerate(ids):
        shp.set_full_model(id, bkgmodels[i])  # bkg_scale included in bkgmodel


def _set_full_model(ids, *args, **kwargs):
    if len(args) == 1:
        _set_full_model_bkg(ids, *args, **kwargs)
    else:
        _set_full_model_src_bkg(ids, *args, **kwargs)


def _check_parameters_different_from_priors(samples, priors, parameters, limit=0.005):
    rng = np.random.default_rng()
    for i, par, prior in zip(count(), parameters, priors):
        prior_sample = [prior(s) for s in rng.random(samples.shape[0])]
        _, pvalue = ks_2samp(samples[:, i], prior_sample)

        if pvalue >= limit:
            logging.warn(
                f"Parameter {par.name} is like prior. "
                f"KS test p-value: {pvalue:0.3f}"
            )


def _fit_all(solver, resume, verbose, **kwargs):
    speed = "safe" # 2*len(parameters)
    results = solver.run(
        resume=resume, verbose=verbose, speed=speed, Lepsilon=0.5, frac_remain=0.1, **kwargs
    )
    best_fit_values = [p.val for p in solver.parameters]

    return results["samples"], best_fit_values


def _results_dir_structure(output_path):
    dirs = ["plots", "info", "extra"]
    paths = {}

    for d in dirs:
        d_path = output_path / d
        if not d_path.exists():
            d_path.mkdir(parents=True)

        paths[d] = d_path

    return paths


def _samples_to_ndarray(samples):
    min_chain_size = min([len(s) for s in samples])
    samples_new = np.zeros((min_chain_size, len(samples)))
    for i, s in enumerate(samples):
        samples_new[:, i] = s[min_chain_size-1::-1, 0]

    return samples_new


def _get_evidence(output_path):
    json_path = output_path.joinpath("info", "results.json")
    return json.load(json_path.open())['logz']


def _save_total_evidence(logZ, output_path):
    json_path = output_path.joinpath("results.json")
    with json_path.open("w") as fp:
        json.dump({"logz": logZ}, fp)


def _fit_independent(ids, parameters, output_path, resume, verbose, **kwargs):
    paths = _results_dir_structure(output_path)

    logZ = 0
    samples, best_fit_values = [], []

    for i, id in enumerate(ids):
        ioutput = output_path / f"bkgid{id}"
        iprior = BXAPrior("background_only", [parameters[i]])
        isolver = bxa.BXASolver(
            id=id,
            prior=iprior.function,
            parameters=[parameters[i]],
            outputfiles_basename=ioutput.as_posix(),
        )
        results = isolver.run(resume=resume, verbose=verbose, **kwargs)
        samples.append(results["samples"])
        best_fit_values.append(parameters[i].val)
        logZ += _get_evidence(ioutput)

    samples = _samples_to_ndarray(samples)
    _save_total_evidence(logZ, paths["info"])

    return samples, best_fit_values


def _run_bxa_fit(
    ids, model_name, parameters, output_path, resume, verbose, **kwargs
):
    prior = BXAPrior(model_name, parameters)

#    if len(parameters) > 3:
#        for i in range(3, len(parameters)):
#            prior.list += [bxa.create_uniform_prior_for(parameters[i])]

    solver = bxa.BXASolver(
        id=ids[0],
        otherids=ids[1:],
        prior=prior.function,
        parameters=parameters,
        outputfiles_basename=output_path.as_posix(),
    )

    if model_name != "background_only":
        samples, best_fit_values = _fit_all(solver, resume, verbose, **kwargs)
    else:
        samples, best_fit_values = _fit_independent(
            ids, parameters, output_path, resume, verbose, **kwargs
        )

    _check_parameters_different_from_priors(samples, prior.list, parameters)

    gof = _goodness_fixed(ids)

#    if model_name != "background_only":
#        ppp = _goodness_ppp(ids, parameters, samples, output_path=output_path)
#        gof.update(ppp)

    solver.set_best_fit()

    return samples, best_fit_values, solver, gof
#    return samples, best_fit_values, solver


def _get_model_disp(id, samples, solver, percentiles, model_bkg_scaled, filename):
    model_disp_path = Path(solver.outputfiles_basename, "extra", filename)

    if not model_disp_path.exists():
        model_disp = utils.calc_model_dispersion(
            id, samples, solver.parameters, model_bkg_scaled, percentiles
        )
        np.savetxt(model_disp_path, model_disp)
    else:
        model_disp = np.loadtxt(model_disp_path)

    return model_disp


def _get_data_for_plots(ids, samples, solver, model_name, best_fit_values, coadd=False):
    # percentiles = [15.87, 50, 84.13]  # -1sigma, median, 1sigma
    percentiles = [3.3, 50, 97.7]  # -2sigma, median, 2sigma
    data_for_plots = []

    for i, id in enumerate(ids):
        shp.notice_id(id)
        shp.ungroup(id)

        data_for_plots.append({})
        current = data_for_plots[i]

        current["model"] = deepcopy(shp.get_fit_plot(id).modelplot)
        if coadd:
            model_bkg_scaled = np.zeros_like(current["model"].y)
            disp_filename = f"model_disp_coadd{id}.dat"
        else:
            current["model_bkg"] = deepcopy(shp.get_bkg_fit_plot(id).modelplot)
            current["backscale"] = shp.get_bkg_scale(id)
            model_bkg_scaled = current["backscale"] * current["model_bkg"].y
            disp_filename = f"model_disp_{id}.dat"

#        model_disp = _get_model_disp(
#            id, samples, solver, percentiles, model_bkg_scaled, disp_filename
#        )
#        current["model_disp"] = model_disp[:, :3]
#        current["model_src_disp"] = model_disp[:, 3:]

        utils.set_model_to_values(solver.parameters, best_fit_values)
        utils.group_snr_dataset(id, snr=3)

        current["data"] = deepcopy(shp.get_fit_plot(id).dataplot)
        current["resid"] = deepcopy(shp.get_resid_plot(id))
        if not coadd:
            current["data_bkg"] = deepcopy(shp.get_bkg_fit_plot(id).dataplot)
            current["resid_bkg"] = deepcopy(shp.get_bkg_resid_plot(id))

        if model_name == "background_only":
            current["backscale"] = 0

    return data_for_plots


def _get_degrees_of_freedom(ids, data_points, bkg=False):
    degrees_of_freedom = data_points

    for id in ids:
        if bkg:
            # When the background fitting process finish, all parameters are frozen,
            # so we count all parameters in the model to calculate the dof
            m = shp.get_bkg_model(id)
            free_parameters = m.pars
        else:
            m = shp.get_model(id)
            free_parameters = [p for p in m.pars if not p.frozen]

        degrees_of_freedom -= len(free_parameters)

    return degrees_of_freedom


def _calc_cstat(ids, data, model, bkg=False):
    dof = _get_degrees_of_freedom(ids, len(data), bkg=bkg)
    cstat = stats.cstat(data, model)

    return {"cstat": cstat, "dof": dof}


def _goodness(ids, bkg=False, add_cstat=True, add_ppp=True):
    model_total, data_total = [], []
    for id in ids:
        data_counts, model_counts = utils.get_data_model_counts(id, bkg=bkg)
        model_total = np.concatenate((model_total, model_counts))
        data_total = np.concatenate((data_total, data_counts))

    gof = stats.goodness(data_total, model_total)

    if add_cstat:
        gof.update(_calc_cstat(ids, data_total, model_total, bkg=bkg))

    return gof

def _goodness_fixed(ids, bkg=False, add_cstat=True, add_ppp=True):
    model_total, data_total = [], []
    for id in ids:
        data_counts, model_counts = utils.get_data_model_counts(id, bkg=bkg)
        model_total = np.concatenate((model_total, model_counts))
        data_total = np.concatenate((data_total, data_counts))

    gof = stats.goodness_fixed(data_total, model_total)

    if add_cstat:
        gof.update(_calc_cstat(ids, data_total, model_total, bkg=bkg))

    return gof

def _goodness_ppp(ids, parameters, samples, nsims=1000, output_path=None):
#    logging.info(
#        "Calculating posterior predictive p-value. "
#        f"Running {nsims} simulations..."
#    )
    ppp = stats.posterior_predictive_pvalue(
        ids, parameters, samples, nsims=nsims, output_path=output_path
    )

    return {"ppp": ppp}


def _check_goodness_bkgfit(gof, limit=0.1):
    pvalue = gof["KS"]["p-value"]

    if  pvalue < limit:
        logging.warn(
            "Background model rejected with 90% confidence. "
            f"KS test p-value: {pvalue:0.3f}"
        )


def _save_goodness_results(gof, gof_bkg, output_path):
    info_path = output_path / "info"
    if not info_path.exists():
        info_path.mkdir()

    json_path = info_path.joinpath("goodness.json")

    with json_path.open("w") as fp:
        json.dump({"src": gof, "bkg": gof_bkg}, fp, indent=2)

def _save_goodness_results2(gof, gof_bkg, output_path):
    info_path = output_path / "info"
    if not info_path.exists():
        info_path.mkdir()

    json_path = info_path.joinpath("goodness_fixed.json")

    with json_path.open("w") as fp:
        json.dump({"src": gof, "bkg": gof_bkg}, fp, indent=2)



def _set_energy_bands():
    return [
        [0.5, 2.0],
        [2.0, 10.0]
    ]


def _get_src_flux(samples, solver, model, ebands):
    flux_dist_path = Path(solver.outputfiles_basename, "extra", "src_flux_dist.dat")

    if not flux_dist_path.exists():
        log_flux_dist_cgs = calc_flux_dist(samples, solver.parameters, model, ebands)
        np.savetxt(flux_dist_path, log_flux_dist_cgs)
        solver.set_best_fit()
    else:
        log_flux_dist_cgs = np.loadtxt(flux_dist_path)

    return log_flux_dist_cgs


def _get_src_lumin(samples, solver, model, ebands, z):
    lumin_dist_path = Path(solver.outputfiles_basename, "extra", "src_lumin_dist.dat")

    if not lumin_dist_path.exists():
        log_lumin_dist_cgs = calc_lumin_dist(samples, solver.parameters, model, ebands, z)
        np.savetxt(lumin_dist_path, log_lumin_dist_cgs)
        solver.set_best_fit()
    else:
        log_lumin_dist_cgs = np.loadtxt(lumin_dist_path)

    return log_lumin_dist_cgs


def _coadded_spectra(srcid, galabs, srcmodel, samples, solver, model_name):
    shp.clean()
    coadd_ids = _load_data(srcid, coadd=True)

    for id in coadd_ids:
        shp.set_source(id, model=galabs*srcmodel)

    solver.set_best_fit()
    best_fit_values = [p.val for p in solver.parameters]

    data_for_plots_coadd = _get_data_for_plots(
        coadd_ids, samples, solver, model_name, best_fit_values, coadd=True
    )

    plots.coadd_spectra(
        coadd_ids, data_for_plots_coadd, output_path=solver.outputfiles_basename
    )

def _run_bkg(src, model_name, output_path, resume_fit, verbose, fix_const, kwargs_load_data={}, **kwargs):
    _set_sherpa_env()
    ids = _load_data(src, **kwargs_load_data) # 2023-02-07 AV: passing kwargs to load data

    #logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models(ids, galabs)

#    dirname = f"bkg_fit_{src}/"
    info_path = Path(".")
#    if not info_path.exists():
#        info_path.mkdir()

    json_path = info_path.joinpath(f"goodness_bkg_{src}.json")

    with json_path.open("w") as fp:
        json.dump({"bkg": goodness_stats_bkg}, fp, indent=2)


def _run(
    src,
    model_name,
    output_path,
    resume_fit,
    verbose,
    fix_back,
    fix_inst,
    sp_info,
    args_model=(), # 20230131 AV: added in order to allow for fixing the model parameters e.g. redshift
    kwargs_load_data={},
    use_galabs=False, # 20230208 AV: added to allow fixing NH for extragalactic sources
    use_tbabs_table=False, # 20230216 AV: added to allow use of Angel's tbabs table
    **kwargs
):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env()
    ids = _load_data(src, **kwargs_load_data) # 2023-02-07 AV: passing kwargs to load data


#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(
            model_name,
            ids,
            fix_inst,
            use_tbabs_table=use_tbabs_table,
            *args_model # 20230131 AV: see above
        )
        _set_full_model(ids, galabs, srcmodel, bkgmodels, use_galabs=use_galabs) # 20230208 AV: see above

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    logger = logging.getLogger("ultranest")
    logger.setLevel(logging.ERROR)

    samples, best_fit_values, solver, goodness_stats = _run_bxa_fit(
        ids, model_name, parameters, output_path, resume_fit, verbose, **kwargs
    )

#    samples, best_fit_values, solver = _run_bxa_fit(
#        ids, model_name, parameters, output_path, resume_fit, verbose, **kwargs
#    )

    _save_goodness_results(goodness_stats, goodness_stats_bkg, output_path)

#    logging.info("Making plots...")
#    plots.qqplots(ids, output_path=output_path)
    obsid = shp.get_data(1).name
    plots.plot_corner(samples, parameters, obsid[1:11], src, output_path=output_path)
    plots.parameters(samples, parameters, output_path=output_path)

    data_for_plots = _get_data_for_plots(
        ids, samples, solver, model_name, best_fit_values
    )
    plots.spectra(ids, data_for_plots, obsid[1:11], src, output_path=output_path)

    if sp_info == 1:
        sp_info = _get_spec_info(ids, output_path)
        filename = Path(".", f"spec_info5_{src}.txt")
        with open(filename, 'w', 8192) as testfile:
            for row in sp_info:
                testfile.write(' '.join([str(j) for j in row]) + '\n')


#    if model_name != "background_only":
#        logging.info("Calculating fluxes and luminosities...")
#        ebands = _set_energy_bands()
#        logflux_dist = _get_src_flux(samples, solver, srcmodel[0], ebands)
        #loglumin_dist = _get_src_lumin(samples, solver, srcmodel[0], ebands, 0)

#        plots.flux(logflux_dist, ebands, output_path=output_path)
        #plots.lumin(loglumin_dist, ebands, output_path=output_path)

        #logging.info("Plotting coadded spectra...")
        #_coadded_spectra(src.id, galabs, srcmodel, samples, solver, model_name)

def _run2(src, model_name, output_path, resume_fit, verbose, fix_back, fix_inst, sp_info, **kwargs):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env()
    ids = _load_data(src)

#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models2(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_inst, output_path)
        _set_full_model(ids, galabs, srcmodel, bkgmodels)

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    logger = logging.getLogger("ultranest")
    logger.setLevel(logging.ERROR)

#    samples, best_fit_values, solver, goodness_stats = _run_bxa_fit(
#        ids, model_name, parameters, output_path, resume_fit, verbose, **kwargs
#    )

#    samples, best_fit_values, solver = _run_bxa_fit(
#        ids, model_name, parameters, output_path, resume_fit, verbose, **kwargs
#    )

    goodness_stats = _goodness_fixed(ids)

    _save_goodness_results2(goodness_stats, goodness_stats_bkg, output_path)

#    logging.info("Making plots...")
#    plots.qqplots(ids, output_path=output_path)
#    obsid = shp.get_data(1).name
#    plots.plot_corner(samples, parameters, obsid[1:11], src, output_path=output_path)
#    plots.parameters(samples, parameters, output_path=output_path)

#    data_for_plots = _get_data_for_plots(
#        ids, samples, solver, model_name, best_fit_values
#    )
#    plots.spectra(ids, data_for_plots, obsid[1:11], src, output_path=output_path)

    if sp_info == 1:
    	_get_spec_info(ids, output_path)

#    if model_name != "background_only":
#        logging.info("Calculating fluxes and luminosities...")
#        ebands = _set_energy_bands()
#        logflux_dist = _get_src_flux(samples, solver, srcmodel[0], ebands)
        #loglumin_dist = _get_src_lumin(samples, solver, srcmodel[0], ebands, 0)

#        plots.flux(logflux_dist, ebands, output_path=output_path)
        #plots.lumin(loglumin_dist, ebands, output_path=output_path)

        #logging.info("Plotting coadded spectra...")
        #_coadded_spectra(src.id, galabs, srcmodel, samples, solver, model_name)

def _run_bkg_stat(src, model_name, output_path, resume_fit, verbose, fix_back, fix_inst, sp_info, **kwargs):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env2()
    ids = _load_data2(src)

#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models2(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_inst, output_path)
        _set_full_model(ids, galabs, srcmodel, bkgmodels)

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    shp.calc_stat_info
    res = shp.get_stat_info()

    a = []
    a.append([utils.get_detector(1), res[1].statval, res[1].numpoints, res[1].dof])

    if len(ids) > 1:
        a.append([utils.get_detector(2), res[3].statval, res[3].numpoints, res[3].dof])

    filename = f"stat_info1_{src}.txt"
    with open(filename, 'w', 8192) as testfile:
        for row in a:
            testfile.write(' '.join([str(j) for j in row]) + '\n')

    if sp_info == 1:
    	_get_spec_info(ids, output_path)

def _run_bkg_stat2(src, model_name, output_path, resume_fit, verbose, fix_back, fix_inst, sp_info, **kwargs):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env2()
    ids = _load_data3(src)

#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models2(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    print("model_name is", model_name)
    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_inst, output_path)
        _set_full_model(ids, galabs, srcmodel, bkgmodels)

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    shp.calc_stat_info
    res = shp.get_stat_info()
    detector = utils.get_detector(1)
    if detector == "EPN" or detector == "EPNA":
        fpar = 12
    else:
        fpar = 7

    pval = chi2.sf(res[1].statval, res[1].numpoints-fpar)

    if (res[1].numpoints-fpar) > 0:
        chir = res[1].statval/(res[1].numpoints-fpar)
    else:
        chir = 0

    a = []
    a.append([detector, res[1].statval, res[1].numpoints, pval, res[1].dof, chir])

    if len(ids) > 1:
        detector = utils.get_detector(2)
        if detector == "EPN" or detector == "EPNA":
            fpar = 12
        else:
            fpar = 7

        pval = chi2.sf(res[3].statval, res[3].numpoints-fpar)

        if (res[3].numpoints-fpar) > 0:
            chir = res[3].statval/(res[3].numpoints-fpar)
        else:
            chir = 0

        a.append([detector, res[3].statval, res[3].numpoints, pval, res[3].dof, chir])

    filename = f"stat_info2_{src}.txt"
    with open(filename, 'w', 8192) as testfile:
        for row in a:
            testfile.write(' '.join([str(j) for j in row]) + '\n')

    if sp_info == 1:
    	_get_spec_info(ids, output_path)

def _run_bkg_stat3(src, model_name, output_path, resume_fit, verbose, fix_back, fix_inst, sp_info, **kwargs):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env2()
    ids = _load_data2(src)

#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models2(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_inst, output_path)
        _set_full_model(ids, galabs, srcmodel, bkgmodels)

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    percentiles = [3.3, 50, 97.7]  # -2sigma, median, 2sigma
    data_for_plots = []

    for i, id in enumerate(ids):
        shp.notice_id(id)
#        shp.ungroup(id)

        data_for_plots.append({})
        current = data_for_plots[i]

        current["model"] = deepcopy(shp.get_fit_plot(id).modelplot)
        current["model_bkg"] = deepcopy(shp.get_bkg_fit_plot(id).modelplot)
        current["backscale"] = shp.get_bkg_scale(id)
        model_bkg_scaled = current["backscale"] * current["model_bkg"].y
#            disp_filename = f"model_disp_{id}.dat"

#        model_disp = _get_model_disp(
#            id, samples, solver, percentiles, model_bkg_scaled, disp_filename
#        )
#        current["model_disp"] = model_disp[:, :3]
#        current["model_src_disp"] = model_disp[:, 3:]

#        utils.set_model_to_values(solver.parameters, best_fit_values)
#        utils.group_snr_dataset(id, snr=3)

        current["data"] = deepcopy(shp.get_fit_plot(id).dataplot)
        current["resid"] = deepcopy(shp.get_resid_plot(id))
        current["data_bkg"] = deepcopy(shp.get_bkg_fit_plot(id).dataplot)
        current["resid_bkg"] = deepcopy(shp.get_bkg_resid_plot(id))

    obsid = shp.get_data(1).name
    plots.spectra2(ids, data_for_plots, obsid[1:11], src, output_path=output_path)

    if sp_info == 1:
    	_get_spec_info(ids, output_path)

def _run_bkg_stat4(src, model_name, output_path, resume_fit, verbose, fix_back, fix_inst, sp_info, **kwargs):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env2()
    ids = _load_data4(src)

#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models2(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_inst, output_path)
        _set_full_model(ids, galabs, srcmodel, bkgmodels)

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    shp.calc_stat_info
    res = shp.get_stat_info()
    detector = utils.get_detector(1)
    if detector == "EPN" or detector == "EPNA":
        fpar = 12
    else:
        fpar = 7

    if fpar >= res[1].numpoints:
        fpar = res[1].numpoints - 1

    pval = chi2.sf(res[1].statval, res[1].numpoints-fpar)

    if (res[1].numpoints-fpar) > 0:
        chir = res[1].statval/(res[1].numpoints-fpar)
    else:
        chir = 0

    a = []
    a.append([detector, res[1].statval, res[1].numpoints, pval, res[1].dof, chir])

    if len(ids) > 1:
        detector = utils.get_detector(2)
        if detector == "EPN" or detector == "EPNA":
            fpar = 12
        else:
            fpar = 7

        if fpar >= res[1].numpoints:
            fpar = res[1].numpoints - 1

        pval = chi2.sf(res[3].statval, res[3].numpoints-fpar)

        if (res[3].numpoints-fpar) > 0:
            chir = res[3].statval/(res[3].numpoints-fpar)
        else:
            chir = 0

        a.append([detector, res[3].statval, res[3].numpoints, pval, res[3].dof, chir])

    filename = f"stat_info3_{src}.txt"
    with open(filename, 'w', 8192) as testfile:
        for row in a:
            testfile.write(' '.join([str(j) for j in row]) + '\n')

    if sp_info == 1:
    	_get_spec_info(ids, output_path)

def _run_bkg_stat5(src, model_name, output_path, resume_fit, verbose, fix_back, fix_inst, sp_info, **kwargs):

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)

    _set_sherpa_env2()
    ids = _load_data5(src)

#    logging.info("Setting models...")
    ra = shp.get_data(1).header["SRC_CRA"]
    dec = shp.get_data(1).header["SRC_CDEC"]
    coords = SkyCoord(ra, dec, unit='deg')
    galabs = models.galactic_absorption(coords)
    #galabs = shp.xstbabs.galabs
    #galabs.nH = 0.01
    #galabs.nH.freeze()

    bkgmodels, goodness_stats_bkg = _background_models2(ids, galabs)

    if fix_back != 1:
        backcons = []
        for i, id in enumerate(ids):
            backcons.append(bkgmodels[i].pars[0])
            backcons[i].thaw()
            backcons[i].max = 10.0

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_inst, output_path)
        _set_full_model(ids, galabs, srcmodel, bkgmodels)

    if fix_back != 1:
        for i, id in enumerate(ids):
            parameters.append(backcons[i])

    shp.calc_stat_info
    res = shp.get_stat_info()
    detector = utils.get_detector(1)
    if detector == "EPN" or detector == "EPNA":
        fpar = 12
    else:
        fpar = 7

    pval = chi2.sf(res[1].statval, res[1].numpoints-fpar)

    if (res[1].numpoints-fpar) > 0:
        chir = res[1].statval/(res[1].numpoints-fpar)
    else:
        chir = 0

    a = []
    a.append([detector, res[1].statval, res[1].numpoints, pval, res[1].dof, chir])

    if len(ids) > 1:
        detector = utils.get_detector(2)
        if detector == "EPN" or detector == "EPNA":
            fpar = 12
        else:
            fpar = 7

        pval = chi2.sf(res[3].statval, res[3].numpoints-fpar)

        if (res[3].numpoints-fpar) > 0:
            chir = res[3].statval/(res[3].numpoints-fpar)
        else:
            chir = 0

        a.append([detector, res[3].statval, res[3].numpoints, pval, res[3].dof, chir])

    filename = f"stat_info5_{src}.txt"
    with open(filename, 'w', 8192) as testfile:
        for row in a:
            testfile.write(' '.join([str(j) for j in row]) + '\n')

    if sp_info == 1:
        sp_info = _get_spec_info(ids, output_path)
        filename = Path(".", f"spec_info5_{src}.txt")
        with open(filename, 'w', 8192) as testfile:
            for row in sp_info:
                testfile.write(' '.join([str(j) for j in row]) + '\n')

def _run_wpca(src, model_name, output_path, resume_fit, verbose, fix_const, **kwargs):
    _set_sherpa_env()
    ids = _load_data(src)

    logging.info("Setting models...")
    #galabs = models.galactic_absorption(src.coords)
    galabs = shp.xstbabs.galabs
    galabs.nH = 0.01
    galabs.nH.freeze()

    if model_name == "background_only":
        srcmodels, parameters = models.background_only(bkgmodels)
        _set_full_model(ids, srcmodels)
    else:
        srcmodel, parameters = models.get_source_model(model_name, ids, fix_const)

    bkgmodels, goodness_stats_bkg = _background_models_pca(src, ids, srcmodel, model_name)

    backcons = []
    for i, id in enumerate(ids):
        backcons.append(bkgmodels[i].pars[0])
        backcons[i].thaw()
        backcons[i].max = backcons[i].val + 2
        backcons[i].min = backcons[i].val - 2
        parameters.append(backcons[i])

    _set_full_model(ids, galabs, srcmodel, bkgmodels)

    samples, best_fit_values, solver, goodness_stats = _run_bxa_fit(ids, model_name, parameters, output_path, resume_fit, verbose, **kwargs)

    _save_goodness_results(goodness_stats, goodness_stats_bkg, output_path)

    logging.info("Making plots...")
    plots.qqplots(ids, output_path=output_path)
    plots.parameters(samples, parameters, output_path=output_path)

    data_for_plots = _get_data_for_plots(ids, samples, solver, model_name, best_fit_values)

    plots.spectra(ids, data_for_plots, output_path=output_path)

    if model_name != "background_only":
        logging.info("Calculating fluxes and luminosities...")
        ebands = _set_energy_bands()
        logflux_dist = _get_src_flux(samples, solver, srcmodel[0], ebands)
        loglumin_dist = _get_src_lumin(samples, solver, srcmodel[0], ebands, 0)

        plots.flux(logflux_dist, ebands, output_path=output_path)
        plots.lumin(loglumin_dist, ebands, output_path=output_path)


def fit_spectra(
    src_row, model_name="torus", scratch=False, resume_fit=True, verbose=False, **kwargs
):
    src = Source(src_row, syspdf_params=(0.21, 0.03))
    output_path = _set_output_path(src.id, model_name)

    logs.set_logger(src.id, model_name, stdout_to_log=True)

    if _check_fit_exists(output_path) and not scratch:
        logging.info("Fit available, skipping source.")
    else:
        try:
            _run(src, model_name, output_path, resume_fit, verbose, **kwargs)

        except Exception as e:
            logs.log_exception(e)
            logging.error(f"Fitting source {src.id} failed!")
