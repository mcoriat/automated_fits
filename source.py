from functools import cached_property
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d

from utils import syspdfz


class Source:
    pdfs_path = Path(".", "zphot_pdfs")

    def __init__(self, src, syspdf_params=None):
        self._ra, self._dec = self._set_radec(src)
        self._zcdf = None
        self._syspdf_params = syspdf_params

        self.id = self._set_srcid(src)
        self.z = self._set_redshift(src)
        self.coords = self._set_coords()

    def _set_srcid(self, src):
        return src["srcid"]

    def _set_radec(self, src):
        return src["Ctpra"], src["Ctpdec"]

    def _set_redshift(self, src):
        return src["Zspec"]

    def _set_coords(self):
        return SkyCoord(self._ra, self._dec, unit="deg")

    @cached_property
    def zcdf(self):
        if self.z < 0:
            self._zcdf = self._calc_zcdf()

        return self._zcdf

    def _calc_zcdf(self):
        zpdf = np.loadtxt(self.pdfs_path.joinpath(f"pdfZ_{self.id}.dat"), delimiter=",")

        if self._syspdf_params is not None:
            zpdf[:, 1] = self._syspdf(zpdf, *self._syspdf_params)

        ic = interp1d(
            zpdf[:, 0], zpdf[:, 1], bounds_error=False, fill_value=0, kind="linear"
        )
        
        step = 0.01
        smooth_zpdfx = np.arange(0, zpdf[-1, 0], step=step)
        smooth_zpdfy = ic(smooth_zpdfx)

        zcdf = np.zeros((len(smooth_zpdfx), 2))
        zcdf[:, 0] = smooth_zpdfx
        zcdf[:, 1] = smooth_zpdfy.cumsum() * step

        return zcdf

    @staticmethod
    def _syspdf(pdf, eta, sigma):
        return [syspdfz(z, eta, sigma, pdf) for z in pdf[:, 0]]
