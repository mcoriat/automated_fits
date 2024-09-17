import bxa.sherpa as bxa
import numpy as np
from scipy.interpolate import interp1d
from scipy import signal
import re


class BXAPrior:
    def __init__(self, model, *args, **kwargs):
        self._priors = self._set_priors(model, *args, **kwargs)
        self._function = self._set_prior_function()

    @property
    def list(self):
        return self._priors

    @property
    def function(self):
        return self._function

    def _set_priors(self, model, parameters, zcdf=None):
        prior_method = self._get_prior_method(model)
        priors = prior_method(parameters)

     #   if not model == "background_only":
     #        priors = self._add_zprior(priors, zcdf)

        # 20230209 AV: what the... commenting this out
        #if len(parameters) > 3:
        #    for i in range(3, len(parameters)):
        #         priors += [bxa.create_uniform_prior_for(parameters[i])]
#                priors += [bxa.create_gaussian_prior_for(parameters[i], 1, 0.25)]

        # 20230209 AV: the priors for constants is not automatically defined in
        # the priors so define them here... I think this is the logic of the
        # above piece of code but for models with more than 3 parameters it
        # will erase loguniform priors..
        for p in parameters:
            if not re.match("constant_[0-9]", p.modelname):
                continue
            priors += [bxa.create_uniform_prior_for(p)]

        return priors

    def _get_prior_method(self, model):
        try:
            prior_method = getattr(self, f"_set_priors_{model}")

        except AttributeError:
            raise ValueError(f"Priors for '{model}' are not defined!!!")

        return prior_method

    def _set_priors_powerlaw(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_uniform_prior_for(parameters[2]),               # PhoIndex
        ]

    def _set_priors_blackbody(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_loguniform_prior_for(parameters[2]),            # Temperature
        ]

    def _set_priors_apec(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_loguniform_prior_for(parameters[2]),            # Temperature
        ]

    def _set_priors_powerlaw2(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_uniform_prior_for(parameters[2]),               # PhoIndex
        ]

    def _set_priors_powerlaw3(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_uniform_prior_for(parameters[2]),               # PhoIndex
        ]


    def _set_priors_zpowlaw(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_uniform_prior_for(parameters[2]),               # PhoIndex
        ]
        #return [
        #    bxa.create_uniform_prior_for(parameters[0]),
        #    bxa.create_gaussian_prior_for(parameters[1], 1.95, 0.15),
        #    bxa.create_uniform_prior_for(parameters[2]),
        #]


    def _set_priors_double_zpowlaw(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux1
            bxa.create_uniform_prior_for(parameters[1]),               # logNH1
            bxa.create_uniform_prior_for(parameters[2]),               # PhoIndex1
            bxa.create_uniform_prior_for(parameters[3]),               # logflux2
            # AV:   commenting out the following since the second obscuration
            #       should not be in use anymore
            #bxa.create_uniform_prior_for(parameters[4]),               # logNH2
            #bxa.create_uniform_prior_for(parameters[3]),               # PhoIndex2
            bxa.create_uniform_prior_for(parameters[4]),               # PhoIndex2
        ]
        
    def _set_priors_torus(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),
            bxa.create_gaussian_prior_for(parameters[1], 1.95, 0.15),
            bxa.create_uniform_prior_for(parameters[2]),
            bxa.create_uniform_prior_for(parameters[3]),
            bxa.create_uniform_prior_for(parameters[4]),
        ]

    def _set_priors_bremss(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_loguniform_prior_for(parameters[2]),            # kT
        ]

    def _set_priors_apec_single(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_loguniform_prior_for(parameters[2]),            # kT
        ]

    def _set_priors_apec_apec(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_loguniform_prior_for(parameters[2]),            # kT1
            bxa.create_loguniform_prior_for(parameters[3]),            # kT2
        ]

    def _set_priors_apec_apec_const(self, parameters):
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_loguniform_prior_for(parameters[2]),            # kT1
            bxa.create_loguniform_prior_for(parameters[3]),            # kT2
            bxa.create_uniform_prior_for(parameters[4]),               # constant
        ]

    def _set_priors_powerlaw_blackbody(self, parameters):
        print(parameters)
        return [
            bxa.create_uniform_prior_for(parameters[0]),               # logflux
            bxa.create_uniform_prior_for(parameters[1]),               # logNH
            bxa.create_uniform_prior_for(parameters[2]),               # PhoIndex
            bxa.create_loguniform_prior_for(parameters[3]),            # kT
        ]

    def _set_priors_powerlaw_blackbody_const(self, parameters):
        print(parameters)
        return [
            bxa.create_uniform_prior_for(parameters[0]),    # logflux1
            bxa.create_uniform_prior_for(parameters[1]),    # PhoIndex
            bxa.create_uniform_prior_for(parameters[2]),    # logflux2
            bxa.create_loguniform_prior_for(parameters[3]), # kT
            bxa.create_uniform_prior_for(parameters[4]),    # logNH
        ]

    def _set_priors_background_only(self, parameters):
        return [bxa.create_uniform_prior_for(p) for p in parameters]

    def _add_zprior(self, priors, zcdf):
        if zcdf is not None:
            zprior = self._inv_zpdf_func(zcdf)
            priors += [zprior]

        return priors

    def _set_prior_function(self):
        return bxa.create_prior_function(priors=self._priors)

    @staticmethod
    def _inv_zpdf_func(zcdf):
        return lambda x: np.interp(x, zcdf[:, 1], zcdf[:, 0])


class BXAPriorFromPosterior:
    def __init__(self, samples, z):
        self._shape = self._get_nbins(z)
        self._cdf, self._parvals = self._calc_cdf(samples)
        self._function = self._set_prior_function()

    @property
    def function(self):
        return self._function

    def _calc_cdf(self, samples):
        posterior_pdf, edges = np.histogramdd(samples, density=True, bins=self._shape)
        middle_points =  self._calc_middle_points(edges)
        dV = self._calc_volume_element(edges)
        posterior_cdf = posterior_pdf.cumsum() * dV

        win = signal.windows.hann(10)
        posterior_cdf = signal.convolve(posterior_cdf, win, mode='same') / sum(win)

        idxs = np.arange(0, len(posterior_cdf), dtype=int)
        itp = interp1d(
            posterior_cdf,
            idxs,
            kind="nearest",
            fill_value=(0, len(posterior_cdf) - 1),
            bounds_error=False
        )

        return itp, middle_points

    def _prior(self, x):
        idx = [np.unravel_index(int(i), self._shape) for i in self._cdf(x)]
        params = np.array(
            [
                [self._parvals[j][i[j]] for j in range(len(self._shape))] for i in idx
            ]
        )
        return params

    def _set_prior_function(self):
        return lambda x: self._prior(x)

    @staticmethod
    def _calc_middle_points(edges):
        return [0.5*(e[:-1] + e[1:]) for e in edges]

    @staticmethod
    def _calc_volume_element(edges):
        # This only work for regularly spaced bins
        return np.prod([e[1] - e[0] for e in edges])

    @staticmethod
    def _get_nbins(z):
        # This is for torus
        # TODO: use model argument to select the correct binning
        nbins = (30, 20, 10, 30, 30, 35)

        if z > 0:
            nbins = nbins[:-1]

        return nbins
