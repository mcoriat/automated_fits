from pathlib import Path

#import matplotlib
#matplotlib.use('qt5agg')

import numpy as np
import sherpa.astro.ui as shp
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

import corner
import stats
import utils


def matplotlib_settings():
    plt.rcParams['mathtext.fontset'] = "stix"
    plt.rcParams['mathtext.rm'] = "STIXGeneral"
    plt.rcParams['font.family'] = "STIXGeneral"
    plt.rcParams['font.size'] = 13
    plt.rcParams["axes.formatter.limits"] = (-3, 3)
    plt.rcParams["axes.formatter.use_mathtext"] = True
    plt.rcParams["legend.labelspacing"] = 0.25
    
    
def _ticks_format(value, index):
    """
    get the value and returns the value as:
       integer: [0,99]
       1 digit float: [0.1, 0.99]
       n*10^m: otherwise
    To have all the number of the same size they are all returned as latex strings
    """
    exp = np.floor(np.log10(value))
    base = value / 10 ** exp
    if exp == 0 or exp == 1:
        return r"${0:d}$".format(int(value))
    if exp == -1:
        return r"${0:.1f}$".format(value)
    else:
        if base == 1:
            return r"$10^{{{0:d}}}$".format(int(exp))
        else:
            return r"${0:d}\\times10^{{{1:d}}}$".format(int(base), int(exp))


def _get_title(id):
    #obsid = shp.get_data(id).name.split("/")[2]
    detector = shp.get_data(id).header["INSTRUME"]
    
    #return f"{obsid}: {detector}"
    return f"{detector}"


def _get_title_coadd(id):
    path_str = shp.get_data(id).name.split("/")
    obsid = path_str[2]
    detector = path_str[3][:2]

    return f"{obsid}: {detector}"


def _plot_data(ax, x, y, yerr, with_upper_limits=False):
    ax.loglog(x, y, visible=False)
    
    if with_upper_limits:
        mask = y - yerr > 0
        ax.plot(x[~mask], y[~mask] + yerr[~mask], marker="v", lw=0, color="#132A13")
    else:
        mask = [True] * len(x)

    ax.errorbar(
        x[mask],
        y[mask],
        yerr=yerr[mask],
        lw=0,
        elinewidth=2,
        marker="o",
        color="#132A13"
    )


def _plot_model(ax, x, y, disp, color, ls, area_alpha, label):
    ax.plot(x, y, color=color, ls=ls, lw=2, label=label)
    if disp is not None:
        ax.fill_between(x, disp[:, 0], disp[:, 2], color=color, alpha=area_alpha, lw=0)
        # ax.plot(x, disp[:, 1])


def _plot_src_spec(
    ax,
    d,
    m,
    mbkg,
    backscale,
    xlim,
    title,
    m_disp=None,
    msrc_disp=None,
    add_legend=True,
    with_upper_limits=False,
):
    ax.set_facecolor("#F4FDD9")
    ax.set_title(f"Source [{title}]")

    _plot_data(ax, d.x, d.y, d.yerr, with_upper_limits)
    ymin, ymax = ax.get_ylim()

    # Plot best fit model
    _plot_model(ax, m.x, m.y, m_disp, "#CF5C36", "-", 0.5, "best fit")

    if backscale:
        y_src_model = m.y - backscale * mbkg.y
        _plot_model(ax, m.x, y_src_model, msrc_disp, "#5C415D", ":", 0.3, "src model")
    
        y_bkg_model = backscale * mbkg.y
        _plot_model(ax, m.x, y_bkg_model, None, "#6A8D73", "--", 0.3, "bkg")

    ax.set_xticklabels([])
    ax.set_ylabel(d.ylabel)
    ax.set_xlim(*xlim)
    ax.set_ylim(ymin/10, 10*ymax)
    
    if add_legend:
        ax.legend()
    
    ax.xaxis.set_major_formatter(FuncFormatter(_ticks_format))
    ax.tick_params(axis="x", which="minor", bottom="on")


def _plot_bkg_spec(ax, d, m, xlim, title, m_disp=None):
    ax.set_facecolor("#F4FDD9")
    ax.set_title(f"Background [{title}]")

    _plot_data(ax, d.x, d.y, d.yerr, with_upper_limits=False)
    ymin, ymax = ax.get_ylim()

    _plot_model(ax, m.x, m.y, m_disp, "#6A8D73", "-", 0.5, "bkg")

    ax.set_xticklabels([])
    ax.set_xlim(*xlim)
    ax.set_ylim(ymin, ymax)

    ax.xaxis.set_major_formatter(FuncFormatter(_ticks_format))
    ax.tick_params(axis="x", which="minor", bottom="on")


def _plot_resid(ax, resid, xlim, color="k", ylabel=None):
    ax.set_facecolor("#F4FDD9")

    ax.semilogx(resid.x, resid.y, visible=False)
    ax.errorbar(
        resid.x, resid.y, xerr=resid.xerr/2, yerr=resid.yerr, lw=0, elinewidth=2, marker="o", color=color
    )

    ax.axhline(0, c="#CF5C36", ls="--")
    ax.set_xlabel(resid.xlabel)    
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.set_xlim(*xlim)

    ax.xaxis.set_major_formatter(FuncFormatter(_ticks_format))
    ax.tick_params(axis="x", which="minor", bottom="on")


def spectra(ids, data_for_plots, obsid, src, xlim=(0.3, 10), output_path=None, fmt="png"):
    number_of_plots = len(ids)
    fig = plt.figure(
        figsize=(12, 5*number_of_plots), facecolor="#F4FDD9", constrained_layout=True
    )
    srcd = int(src, 16)
    _ = fig.suptitle(f"{obsid}: {srcd}")
    gs = GridSpec(4*number_of_plots, 2, figure=fig)

    for i, id in enumerate(ids):
        data = data_for_plots[i]["data"]
        data_bkg = data_for_plots[i]["data_bkg"]
        model = data_for_plots[i]["model"]
        model_bkg = data_for_plots[i]["model_bkg"]
#        model_disp = data_for_plots[i]["model_disp"]
#        model_src_disp = data_for_plots[i]["model_src_disp"]
        backscale = data_for_plots[i]["backscale"]
        residuals = data_for_plots[i]["resid"]
        residuals_bkg = data_for_plots[i]["resid_bkg"]
        model_disp = None
        model_src_disp = None        

        if i == 0:
            add_legend = True
        else:
            add_legend = False

        ax_src_spec = fig.add_subplot(gs[4*i:3 + 4*i, 0])
        ax_bkg_spec = fig.add_subplot(gs[4*i:3 + 4*i, 1])
        ax_src_resid = fig.add_subplot(gs[3 + 4*i, 0])
        ax_bkg_resid = fig.add_subplot(gs[3 + 4*i, 1])
    
        title = _get_title(id)
        _plot_src_spec(
            ax_src_spec,
            data,
            model,
            model_bkg,
            backscale,
            xlim,
            title,
            model_disp,
            model_src_disp,
            add_legend,
        )
        _plot_resid(ax_src_resid, residuals, xlim, color="#132A13", ylabel="resid.")

        _plot_bkg_spec(ax_bkg_spec, data_bkg, model_bkg, xlim, title)
        _plot_resid(ax_bkg_resid, residuals_bkg, xlim, color="#132A13")
        
    if output_path:
        filename = Path(output_path, "plots", f"spectra_fit_resid.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()


def coadd_spectra(ids, data_for_plots, xlim=(0.3, 10), output_path=None, fmt="png"):
    number_of_plots = len(ids)
    fig = plt.figure(
        figsize=(6, 5*number_of_plots), facecolor="#F4FDD9", constrained_layout=True
    )
    gs = GridSpec(4*number_of_plots, 1, figure=fig)

    for i, id in enumerate(ids):
        data = data_for_plots[i]["data"]
        model = data_for_plots[i]["model"]
        model_disp = data_for_plots[i]["model_disp"]
        residuals = data_for_plots[i]["resid"]

        if i == 0:
            add_legend = True
        else:
            add_legend = False

        ax_src_spec = fig.add_subplot(gs[4*i:3 + 4*i, 0])
        ax_src_resid = fig.add_subplot(gs[3 + 4*i, 0])
    
        title = _get_title_coadd(id)
        _plot_src_spec(
            ax_src_spec,
            data,
            model,
            None,
            0,
            xlim,
            title,
            model_disp,
            None,
            add_legend,
            with_upper_limits=True,
        )
        _plot_resid(ax_src_resid, residuals, xlim, color="#132A13", ylabel="resid.")
        
    if output_path:
        filename = Path(output_path, "plots", f"coadd_spectra_fit_resid.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()

    
def _plot_histogram(sample, ci, hdpi, xlabel):    
    plt.hist(sample, bins="auto", color="#5C415D", alpha=0.3, density=True)
    
    plt.axvline(np.median(sample), ls=":", color="#132A13", label="median/CI")
    for value in ci:
        plt.axvline(value, ls=":", color="#132A13", alpha=0.6)

    plt.axvline(stats.mode(sample), ls="--", color="#CF5C36", label="mode/HDP")
    for value in hdpi:
        plt.axvline(value, ls="--", color="#CF5C36", alpha=0.6)

    plt.xlabel(xlabel)


def parameters(samples, parameters, output_path=None, fmt="png"):
    ci = stats.credible_intervals(samples)
    hdpi = stats.hdp_intervals(samples)

#    nrows = int(np.ceil(len(parameters) / 3))
#    plt.figure(figsize=(14, 3.5 * nrows), facecolor="#F4FDD9")

    a = []
    for i, p in enumerate(parameters):
 #       plt.subplot(nrows, 3, i + 1, facecolor="#F4FDD9")
 #       _plot_histogram(samples[:, i], ci[:, i], hdpi[:, i], p.name)
        a.append([p.name, np.round(np.median(samples[:, i]), 5), np.round(ci[0, i], 5), np.round(ci[1, i], 5), np.round(stats.mode(samples[:, i]), 5), np.round(hdpi[0, i], 5), np.round(hdpi[1, i], 5)])
        
 #       if i == 0:
 #           plt.legend()

 #   plt.tight_layout()

    if output_path:
  #      filename = Path(output_path, "plots", f"parameters_pdf.{fmt}")
  #      plt.savefig(filename)
        file_path = output_path / "extra"
        if not file_path.exists():
            file_path.mkdir()
            
        filename = Path(output_path, "extra", f"parameters.txt")
        with open(filename, 'w', 8192) as testfile:
            for row in a:
                testfile.write(' '.join([str(j) for j in row]) + '\n')
        
    else:
  #      plt.show()
        print(a)

def plot_corner(samples, parameters, obsid, src, output_path=None, fmt="pdf"):
    hdpi = stats.hdp_intervals(samples)
    m = []
    l = []
    for i, p in enumerate(parameters):
        m.append(stats.mode(samples[:, i]))
        l.append(p.name)

    fig = corner.corner(samples, labels=l, quantiles=[0.05, 0.5, 0.95], show_titles=True, title_kwargs={"fontsize": 12})

    srcd = int(src, 16)
    _ = fig.suptitle(f"{obsid}: {srcd}", weight="bold")
    ndim = len(l)
    axes = np.array(fig.axes).reshape((ndim, ndim))
    for j in range(ndim):
        ax = axes[j, j]
        ax.axvline(m[j], ls=":", color="#132A13", alpha=0.6)
        ax.axvline(hdpi[0, j], ls=":", color="#132A13", alpha=0.6)
        ax.axvline(hdpi[1, j], ls=":", color="#132A13", alpha=0.6)
    
#    corner.overplot_lines(fig, m, ls=":", color="#132A13", alpha=0.6)
#    corner.overplot_lines(fig, hdpi[0, :], ls=":", color="#132A13", alpha=0.6)
#    corner.overplot_lines(fig, hdpi[1, :], ls=":", color="#132A13", alpha=0.6)

    if output_path:
        filename = Path(output_path, "plots", f"corners_pdf.{fmt}")
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()
                
        
def _flux_lumin_histograms(samples, ebands, label="F"):
    ci = stats.credible_intervals(samples)
    hdpi = stats.hdp_intervals(samples)

    nrows = int(np.ceil(len(ebands) / 2))
    plt.figure(figsize=(14, 3.5 * nrows), facecolor="#F4FDD9")

    for i, eband in enumerate(ebands):
        plt.subplot(nrows, 2, i + 1, facecolor="#F4FDD9")
        xlabel = f"log {label}({eband[0]} - {eband[1]} keV)"
        samples = np.array(samples)
        jj = samples[:, i] > -20
        _plot_histogram(samples[jj, i], ci[:, i], hdpi[:, i], xlabel)
                        
        if i == 0:
            plt.legend()

    plt.tight_layout()

    
def flux(samples, ebands, output_path=None, fmt="png"):
    _flux_lumin_histograms(samples, ebands, label="F")

    if output_path:
        filename = Path(output_path, "plots", f"fluxes_pdf.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()


def lumin(samples, ebands, output_path=None, fmt="png"):
    _flux_lumin_histograms(samples, ebands, label="L")

    if output_path:
        filename = Path(output_path, "plots", f"lumin_pdf.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()


def _get_qqplot_data(id, bkg=False):
    data_counts, model_counts = utils.get_data_model_counts(id, bkg=bkg)

    return data_counts.cumsum(), model_counts.cumsum()


def _qqplot(data, model):
    plt.plot(model, data, color="#CF5C36")
    
    cmax = max(model[-1], data[-1])
    plt.plot([0, cmax], [0, cmax], ls="--", color="#5C415D", alpha=0.3)
    plt.xlim(0, cmax)
    plt.ylim(0, cmax)
    
    plt.text(
        cmax/2,
        cmax/2,
        "Data excess",
        rotation=45,
        horizontalalignment="center",
        verticalalignment="bottom",
        alpha=0.5,
    )
    plt.text(
        cmax/2,
        cmax/2,
        "Model excess",
        rotation=45,
        horizontalalignment="center",
        verticalalignment="top",
        alpha=0.5,
    )

    plt.xlabel("Counts [model]")
    plt.ylabel("Counts [data]")
    plt.grid(ls=":")


def qqplots(ids, output_path=None, fmt="png"):
    nrows = len(ids)
    plt.figure(figsize=(9, 4.5 * nrows), facecolor="#F4FDD9")

    for i, id in enumerate(ids):
        title = _get_title(id)
        
        plt.subplot(nrows, 2, 2*i + 1, facecolor="#F4FDD9")
        plt.title(f"Source [{title}]")
        data, model = _get_qqplot_data(id)
        _qqplot(data, model)

        plt.subplot(nrows, 2, 2*i + 2, facecolor="#F4FDD9")
        plt.title(f"Background [{title}]")
        data, model = _get_qqplot_data(id, bkg=True)
        _qqplot(data, model)

    plt.tight_layout()

    if output_path:
        filename = Path(output_path, "plots", f"qqplots.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()


# Functions for plotting the posterior predictive
def _get_data(id):
    channels = shp.get_data(id).x
    counts = shp.get_data(id).counts

    return channels, counts


def _group_data(id, channels, counts):
    grouping = utils.group_snr_dataset(id)
    utils.ungroup_dataset(id)

    channels_grp, groups_width, counts_grp, counts_grp_error = utils.group_counts(
        channels, counts, grouping
    )

    return channels_grp, groups_width, counts_grp, counts_grp_error, grouping


def _group_fake_data(fake_counts, grouping, percentiles):
    channels = range(len(grouping))
    ngroups = sum(grouping > 0)
    fake_counts_grouped = np.zeros((ngroups, fake_counts.shape[1]))

    for i, row in enumerate(fake_counts.T):
        _, _, fake_counts_grouped[:, i], _ = utils.group_counts(channels, row, grouping)

    return np.percentile(fake_counts_grouped, percentiles, axis=1)


def _make_errorboxes(channels_grouped, groups_width, perc_grouped):
    errorboxes = [
        Rectangle((ch - chw/2, cntmin), chw, cntmax - cntmin)
        for ch, chw, cntmin, cntmax in zip(
            channels_grouped, groups_width, perc_grouped[0, :], perc_grouped[2, :]
        )
    ]
    return PatchCollection(errorboxes, facecolor="#CF5C36", alpha=0.5, edgecolor="None")

def _calc_ylims(channels_grp, groups_width, perc_grp):
    mask_good = np.logical_and(
        channels_grp - groups_width >= 0.3, perc_grp[0, :] > 0
    )
    ymin = 0.5 * perc_grp[0, mask_good].min()
    ymax = 2 * perc_grp[2, mask_good].max()
    
    return ymin, ymax


def posterior_predictive(
    id, fake_counts, percentiles=[3.3, 50, 97.7], output_path=None, fmt="png"
):
    channels, counts = _get_data(id)
    channels_grp, groups_width, counts_grp, counts_grp_err, grouping = _group_data(
        id, channels, counts
    )

    perc_grp = _group_fake_data(fake_counts, grouping, percentiles)
    errorbox_patches = _make_errorboxes(channels_grp, groups_width, perc_grp)

    plt.figure(figsize=(7, 4.5), facecolor="#F4FDD9")
    ax = plt.subplot(111, facecolor="#F4FDD9")   
    ax.set_title(_get_title(id))
    
    _plot_data(ax, channels_grp, counts_grp, counts_grp_err)
    ax.add_collection(errorbox_patches)
    
    try:
        ymin, ymax = _calc_ylims(channels_grp, groups_width, perc_grp)

    except ValueError:
        ymin, ymax = ax.get_ylim()

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.3, 10)

    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts/keV")

    ax.xaxis.set_major_formatter(FuncFormatter(_ticks_format))
    ax.tick_params(axis="x", which="minor", bottom="on")
    
    plt.tight_layout()

    if output_path:
        filename = Path(output_path, "plots", f"pp_{id}.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()

    plt.close()
    

def posterior_predictive_stats(cstat_sim, cstat_data, output_path=None, fmt="png"):
    plt.figure(figsize=(7, 4.5), facecolor="#F4FDD9")
    plt.subplot(111, facecolor="#F4FDD9")

    plt.hist(cstat_sim, bins="auto", color="#5C415D", alpha=0.3)
    plt.hist(cstat_data, bins="auto", color="#CF5C36")
    
    plt.xlabel("C-stat")

    plt.tight_layout()

    if output_path:
        filename = Path(output_path, "plots", f"pp_stats.{fmt}")
        plt.savefig(filename)
    else:
        plt.show()
        
    plt.close()
