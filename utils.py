import json
from pathlib import Path

import numpy as np
import sherpa.astro.ui as shp
from astropy.convolution import Gaussian1DKernel, convolve


def get_data_model_counts(id, elo=0.3, ehi=10.0, bkg=False):
    if bkg:
        data = shp.get_bkg(id)
        model = shp.get_bkg_model(id)        
    else:
        data = shp.get_data(id)
        model = shp.get_model(id)

    mask = np.logical_and(data.x >= elo, data.x < ehi)
    
    return data.counts[mask], model(0)[mask]


def calc_model_dispersion(id, samples, parameters, bkg_model_y_scaled, percentiles):
    len_ypreds = len(bkg_model_y_scaled)
    ypreds = np.zeros((len_ypreds, len(samples)))

    for j, row in enumerate(samples):
        set_model_to_values(parameters, row)

        m = shp.get_fit_plot(id)
        ypreds[:, j] = m.modelplot.y

    ypreds_src = ypreds - bkg_model_y_scaled[:, None]

    model_disp = np.zeros((len_ypreds, 6))
    model_disp[:, :3] = np.percentile(ypreds, percentiles, axis=1).T
    model_disp[:, 3:] = np.percentile(ypreds_src, percentiles, axis=1).T

    return model_disp


def set_model_to_values(parameters, values):
    for p, v in zip(parameters, values):
        p.val = v


def group_snr_dataset(id, snr=3, elo=0.3, ehi=10.0):
    shp.ignore_id(id, lo=None, hi=elo)
    shp.ignore_id(id, lo=ehi, hi=None)
    good_channels = shp.get_data(id).mask
    print(id)
    shp.group_snr(id, snr=3, tabStops=~good_channels)
    
    return shp.get_data(id).grouping


def ungroup_dataset(id, elo=0.3, ehi=10.0):
    shp.ungroup(id)
    shp.notice_id(id)
    shp.ignore_id(id, lo=None, hi=elo)
    shp.ignore_id(id, lo=ehi, hi=None)


def group_counts(channels, counts, grouping):
    counts_grp, channel_min, channel_max, nchannels = [], [], [], []
    for cnts, channel, g in zip(counts, channels, grouping):
        if g > 0:
            counts_grp.append(cnts)
            channel_min.append(channel)
            channel_max.append(channel)
            nchannels.append(1)
        else:
            counts_grp[-1] += cnts
            channel_max[-1] = channel
            nchannels[-1] += 1
            
    channel_min, channel_max, nchannels = (
        np.array(channel_min), np.array(channel_max), np.array(nchannels)
    )

    channels_grp = (channel_min + channel_max) / 2
    channels_grp_width = channel_max - channel_min
    
    counts_grp_err = np.sqrt(counts_grp)/nchannels
    counts_grp = counts_grp/nchannels
    
    return channels_grp, channels_grp_width, counts_grp, counts_grp_err


def model_compare(srcid, model_names, data_path, limit=30):
    # limit = 30 # for example, Jeffreys scale for the Bayes factor
    info_path = data_path.joinpath(srcid)
    models = {
        p: json.load(
            info_path.joinpath(f"bxafit_{p}", "info", "results.json").open()
        )['logz'] for p in model_names
    }
    
    best = max(models, key=models.__getitem__)
    Zbest = models[best]
    for m in models:
        models[m] -= Zbest
    Ztotal = np.log(sum(np.exp([Z for Z in models.values()])))
    
    good_models = []
    print()
    print(f'Model comparison (source {srcid})')
    print('**********************************')
    print()
    for m in sorted(models, key=models.__getitem__):
        Zrel = models[m]
        if Zrel >= Ztotal - np.log(limit):
            good_models.append(m)
            
        print('model %-10s: log10(Z) = %7.1f %s' % (m, Zrel / np.log(10),
    		' XXX ruled out' if Zrel < Ztotal - np.log(limit) else '   <-- GOOD' ))
    print()
    print('The last, most likely model was used as normalization.')
    print('Uniform model priors are assumed, with a cut of log10(%s) to rule out models.' % limit)
    print()
    
    return good_models


def syspdfz(z, eta, sigma, pdf):
    # pdf: two columns: z, probability density
    flat = np.ones(len(pdf)) / pdf[:, 0].max()
    # The stddev value for the gaussian kernel must be divided by the
    # size of the z-step in the pdf, since the stddev values is relative
    # to the number of elements in the pdf array, not their actual x-values
    dz = pdf[1, 0] - pdf[0, 0]
    sigma_corr = (1 + z)*sigma/dz
    convolved_pdf = convolve(pdf[:, 1], Gaussian1DKernel(stddev=sigma_corr))
    syspdf = eta * flat + (1 - eta) * convolved_pdf

    # fast? nearest neighbour interpolation
    #zgrid = np.linspace(0, zmax, num=len(pdf))
    #zidx = np.round((len(pdf) - 1) * z/zmax).astype(int)

    return np.interp(z, pdf[:, 0], syspdf)


def get_detector(id):
    detectors = {"PN": "EPN", "M1": "EMOS1", "M2": "EMOS2", "MO": "EMOS", "EPNA": "EPNA"}
    detector = shp.get_data(id).header["INSTRUME"]
    
    if detector in detectors.values():
        return detector
    else:
        filename = Path(shp.get_data(id).name).name
        return detectors[filename[:2]]
