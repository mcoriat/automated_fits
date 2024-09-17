from itertools import repeat, count

import numpy as np

import simulations as sim
from utils import get_data_model_counts
from utils import get_detector

def mode(inputData, axis=None, dtype=None):
    """
    Robust estimator of the mode of a data set using the half-sample mode [1].
    Based on functions from Henry Freudenriech (Hughes STX) statistics library
    (called ROBLIB) from AstroIDL User's Library [2].
    
    [1] Robertson & Cryer, 1974
    [2] https://fornax.phys.unm.edu/lwa/subversion/trunk/lsl/lsl/statistics/robust.py
    .. versionadded: 1.0.3
    """
	
    if axis is not None:
        fnc = lambda x: mode(x, dtype=dtype)
        dataMode = np.apply_along_axis(fnc, axis, inputData)
    else:
        # Create the function that we can use for the half-sample mode
        def _hsm2(data):
            if len(data) < 3:
                return data.mean()

            # Round up to include the middle value, in the case of an odd-length array
            half_idx = int((len(data) + 1)/2) 

            # Calculate all interesting ranges that span half of all data points
            ranges = data[-half_idx:] - data[:half_idx]
            smallest_range_idx = np.argmin(ranges)

            # Now repeat the procedure on the half that spans the smallest range
            data_subset = data[smallest_range_idx : (smallest_range_idx + half_idx)]

            return _hsm2(data_subset)
				
        data = inputData.ravel()
        if type(data).__name__ == "MaskedArray":
            data = data.compressed()
        if dtype is not None:
            data = data.astype(dtype)
			
        # Sort and remove repeated values 
        # (the algorithm doesn't work with unsorted or repeated values)
        data = np.unique(data)
		
        # Find the mode
        dataMode = _hsm2(data)

    return dataMode


def hdp_intervals(samples, alpha=0.1):
    """
    Chen-Shao algorithm for calculating the 100(1 - alpha)% Highest Probability Density
    interval of sample, a chain sampling the posterior probability distribution.
    Implementation of the algorithm as described in [1]. It asumes that sample is
    ergodic and the posterior is unimodal.
    
    TODO: in [1] they suggest a possible extension for multimodal distributions,
    but I don't understand what they propose
    
    [1] M. Chen, Q. Shao & J.G. Ibrahim.
        Monte Carlo Methods in Bayesian Computation, Springer 2000. 
        https://doi.org/10.1007/978-1-4612-1276-8

    Parameters
    ----------
    sample : ndarray or list
        chain sampling the posterior distribution.
    alpha : float
        significance level.
    """
    n, m = samples.shape
    intervals_size = int((1 - alpha) * n)
    
    hdpi = np.zeros((2, m))
    for i in range(m):
        ss = np.sort(samples[:, i])
        R = {
            ss[j + intervals_size] - ss[j]: (ss[j], ss[j + intervals_size])
            for j in range(n - intervals_size)
        }        
        hdp_key = np.min(list(R))
        hdpi[:, i] = R[hdp_key]
        
    return hdpi


def credible_intervals(samples, alpha=0.1):
    return np.percentile(samples, [100*alpha/2, 100*(1 - alpha/2)], axis=0)


def cstat(data, model, trunc_value=1.0e-25):
    D = trunc_value * np.ones_like(data)
    D[data > 0] = data[data > 0]

    M = trunc_value * np.ones_like(model)
    M[model > 0] = model[model > 0]

    return 2 * np.sum(M - D + D*(np.log(D) - np.log(M)))


def ks_2samp(sample1, sample2):
    """
    Kolmogorov-Smirnov test for two samples.
    """
    cdf1 = sample1.cumsum()/sample1.sum()
    cdf2 = sample2.cumsum()/sample2.sum()

    return np.abs(cdf1 - cdf2).max()


def ad_2samp(sample1, sample2):
    """
    Implementation of the Anderson-Darling test for two samples.

    See Scholz, F. W and Stephens, M. A. (1987), K-Sample Anderson-Darling Tests, 
    Journal of the American Statistical Association, Vol. 82, pp. 918-924    
    """
    cdf1 = sample1.cumsum()/sample1.sum()
    cdf2 = sample2.cumsum()/sample2.sum()

    m = len(sample1)
    n = len(sample2)
    N = m + n
    HN = (m*cdf1 + n*cdf2)/N

    mask = np.logical_and(HN > 0, HN < 1)    
    I = (cdf1[mask] - cdf2[mask])**2/(HN[mask]*(1 - HN[mask]))
    dHN = np.diff(HN[mask])
    
    return np.sum(I[:-1]*dHN)*(m*n/N)


def goodness(data, model, niter=1000):
    """
    KS and Anderson-Darling test between data and model. If the null hypothesis 
    is not rejected with high probability, the fit is good. 

    Note that in this kind of tests the parameters of the model should be
    independent of the data. If they are estimated by fitting the data,
    the probabilities of the KS test (and other non-paramatric test based on 
    comparison of the empirical cumulative distributions) are wrong!!!
    Fortunately, we can estimate the p-value using non-parametric bootstrap 
    resampling. 

    To do this we used a permutation test, spliting the data+model
    sample in two equal size sample and estimating the statistic.
    Doing this N times, we can estimate the p-value as how often the
    statistic of the permutation is larger than the observed value in the
    original data/model test. 

    See https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
    and Babu, G. J.  & Rao, C. R. (2004) and
    https://stats.stackexchange.com/questions/59774/test-whether-variables-follow-the-same-distribution/59875#59875
    """
    # Estimate KS and AD stats
    ks = ks_2samp(data, model)
    ad = ad_2samp(data, model)
    
    # Estimate the p-values by bootstraping
    count = np.array([0, 0])
    bs_data = np.concatenate((data, model))

    rng = np.random.default_rng()
    for _ in repeat(None, niter):
        perm1, perm2 = np.split(rng.permutation(bs_data), 2)
        bs_ks = ks_2samp(perm1, perm2)
        bs_ad = ad_2samp(perm1, perm2)

        count[0] += (bs_ks >= ks)
        count[1] += (bs_ad >= ad)

    pvals = 1.*count/niter

    goodness = {
        "KS": {
            "stat": ks,
            "p-value": pvals[0],
        },
        "AD": {
            "stat": ad,
            "p-value": pvals[1],
        },
    }
    return goodness

def goodness_fixed(data, model, niter=1000):
    """
    KS and Anderson-Darling test between data and model. If the null hypothesis 
    is not rejected with high probability, the fit is good. 

    Note that in this kind of tests the parameters of the model should be
    independent of the data. If they are estimated by fitting the data,
    the probabilities of the KS test (and other non-paramatric test based on 
    comparison of the empirical cumulative distributions) are wrong!!!
    Fortunately, we can estimate the p-value using non-parametric bootstrap 
    resampling. 

    To do this we used a permutation test, spliting the data+model
    sample in two equal size sample and estimating the statistic.
    Doing this N times, we can estimate the p-value as how often the
    statistic of the permutation is larger than the observed value in the
    original data/model test. 

    See https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
    and Babu, G. J.  & Rao, C. R. (2004) and
    https://stats.stackexchange.com/questions/59774/test-whether-variables-follow-the-same-distribution/59875#59875
    """
    # Estimate KS and AD stats
    ks = ks_2samp(data, model)
    ad = ad_2samp(data, model)
    
    # Estimate the p-values by bootstraping
    count = np.array([0, 0])

    rng = np.random.default_rng()
    for _ in repeat(None, niter):        
        perm1 = np.zeros_like(data)
        perm2 = np.zeros_like(data)

        mask = rng.choice(a=[False, True], size=len(data))
        perm1[mask] = data[mask]
        perm1[~mask] = model[~mask]
        perm2[mask] = model[mask]
        perm2[~mask] = data[~mask]

        bs_ks = ks_2samp(perm1, perm2)
        bs_ad = ad_2samp(perm1, perm2)

        count[0] += (bs_ks >= ks)
        count[1] += (bs_ad >= ad)

    pvals = 1.*count/niter

    goodness = {
        "KS": {
            "stat": ks,
            "p-value": pvals[0],
        },
        "AD": {
            "stat": ad,
            "p-value": pvals[1],
        },
    }
    return goodness


def _sample_posterior(samples, size=None):
    rng = np.random.default_rng()
    randidx = rng.integers(len(samples), size=size)

    return samples[randidx, :]


def _calc_pp_statistic(id, full_model):
    sim.set_full_model(full_model)
            
    counts_sim, model_sim = get_data_model_counts("sim_src")
    counts_data, model_data = get_data_model_counts(id)

    cstat_sim = cstat(counts_sim, model_sim)
    cstat_data = cstat(counts_data, model_data)
    
    return cstat_sim, cstat_data


def posterior_predictive_pvalue(ids, parameters, samples, nsims=1000, output_path=None):
    import plots

    # See https://stats.stackexchange.com/a/70208
    cstat_sim = np.zeros((len(ids), nsims))
    cstat_data = np.zeros((len(ids), nsims))
    
    parameters_values_for_sims = _sample_posterior(samples, size=nsims)
    
    for i, id in enumerate(ids):    
        m = sim.get_full_model(id)
        detector = get_detector(id)

        if detector == "EMOS" or  detector == "EPNA":
            
            srcm, rmf_src = sim.get_src_model_noarf(m)
            bkgm_particles, rmf_particles = sim.get_bkg_particles_model_noarf(m)
            bkgm_astro, rmf_astro = sim.get_bkg_astro_model_noarf(m)

            texp_src, backscal_src = sim.get_src_parameters(id)
            texp_bkg, backscal_bkg = sim.get_bkg_parameters(id)

            nchannels = sim.get_number_of_channels(id)
            fake_src_counts = np.zeros((nchannels, nsims))
        
            for j, parvals in zip(count(), parameters_values_for_sims):
            	sim.delete_simulation()
            
            	sim.fake_bkg_particles_noarf(
                    bkgm_particles, rmf_particles, texp_bkg, backscal_bkg
            )
            	sim.fake_astro_particles_noarf(
                    bkgm_astro, rmf_astro, texp_bkg, backscal_bkg
            )
            	# Values for the model parameters sampled from
            	# the posterior are assigned inside 'fake_src'
            	fake_src_counts[:, j] = sim.fake_src_noarf(
            	    srcm, rmf_src, texp_src, backscal_src, parameters, parvals
            )
            
            	cstat_sim[i, j], cstat_data[i, j] = _calc_pp_statistic(id, m)

#            plots.posterior_predictive(id, fake_src_counts, output_path=output_path)

 
        else:

            srcm, rmf_src, arf_src = sim.get_src_model(m)
            bkgm_particles, rmf_particles, arf_particles = sim.get_bkg_particles_model(m)
            bkgm_astro, rmf_astro, arf_astro = sim.get_bkg_astro_model(m)

            texp_src, backscal_src = sim.get_src_parameters(id)
            texp_bkg, backscal_bkg = sim.get_bkg_parameters(id)

            nchannels = sim.get_number_of_channels(id)
            fake_src_counts = np.zeros((nchannels, nsims))
        
            for j, parvals in zip(count(), parameters_values_for_sims):
            	sim.delete_simulation()
            
            	sim.fake_bkg_particles(
                    bkgm_particles, rmf_particles, arf_particles, texp_bkg, backscal_bkg
            )
            	sim.fake_astro_particles(
                    bkgm_astro, rmf_astro, arf_astro, texp_bkg, backscal_bkg
            )
            	# Values for the model parameters sampled from
            	# the posterior are assigned inside 'fake_src'
            	fake_src_counts[:, j] = sim.fake_src(
            	    srcm, rmf_src, arf_src, texp_src, backscal_src, parameters, parvals
            )
            
            	cstat_sim[i, j], cstat_data[i, j] = _calc_pp_statistic(id, m)

#            plots.posterior_predictive(id, fake_src_counts, output_path=output_path)

    sim.delete_simulation()
    
    cstat_total_sim = np.sum(cstat_sim, axis=0)
    cstat_total_data = np.sum(cstat_data, axis=0)
#    plots.posterior_predictive_stats(
#        cstat_total_sim, cstat_total_data, output_path=output_path
#    )
    
    return np.sum(cstat_total_sim >= cstat_total_data) / len(cstat_total_sim)
