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
