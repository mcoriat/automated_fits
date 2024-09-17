def _set_sherpa_env2(stat="chi2gehrels"):
    shp.clean()
    shp.set_stat(stat)
    shp.set_xsabund('wilm')
    shp.set_xsxsect('vern')

    plots.matplotlib_settings()

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
