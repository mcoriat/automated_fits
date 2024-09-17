def _set_sherpa_env(stat="cstat"):
    shp.clean()
    shp.set_stat(stat)
    shp.set_xsabund('wilm')
    shp.set_xsxsect('vern')

    plots.matplotlib_settings()

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

def _save_goodness_results(gof, gof_bkg, output_path):
    info_path = output_path / "info"
    if not info_path.exists():
        info_path.mkdir()

    json_path = info_path.joinpath("goodness.json")

    with json_path.open("w") as fp:
        json.dump({"src": gof, "bkg": gof_bkg}, fp, indent=2)

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
