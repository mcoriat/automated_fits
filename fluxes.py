import numpy as np
import sherpa.astro.ui as shp
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo

from utils import set_model_to_values


def _get_model_flux(model):
    s = shp.get_source_component_plot(1, model=model)
    FE = s.y  # u.s**-1 * u.cm**-2 * u.keV**-1
    dE = (s.xhi - s.xlo)  # u.keV
    E = s.x  # u.keV

    return E, E * FE * dE


def _calc_flux(model, parameters, values, ebands):
    set_model_to_values(parameters, values)

    E, F = _get_model_flux(model)  # in u.keV / u.s / u.cm**2 

    return [np.sum(F[np.logical_and(E >= e[0], E < e[1])]) for e in ebands]


def calc_flux_dist(samples, parameters, model, ebands):
    flux_dist = np.array(
        [
            _calc_flux(model, parameters, row, ebands) for row in samples
        ]
    )
    flux_dist_units = flux_dist << u.keV / u.s / u.cm**2 

    return np.log10(flux_dist_units.to(u.erg / u.s / u.cm**2).value)  # in CGS units


def _calc_lumin(model, parameters, values, ebands, z):
    for p, v in zip(parameters, values):
        if "logNH" in p.name:
            p.val = p.min
        else:
            p.val = v
            
        if "redshift" in p.name:
            z = v

    E, F = _get_model_flux(model)  # in u.keV / u.s / u.cm**2 
    dl = cosmo.luminosity_distance(z).to(u.cm).value

    flux = np.array(
        [
            np.sum(F[np.logical_and(E >= e[0] / (1 + z), E < e[1] / (1 + z))])
            for e in ebands
        ]
    )
    return 4 * np.pi * dl**2 * flux


def calc_lumin_dist(samples, parameters, model, ebands, z):
    lumin_dist = np.array(
        [
            _calc_lumin(model, parameters, row, ebands, z) for row in samples
        ]
    )    
    # remove negative luminosities in case some redshift value is negative
    # weird case, but it happens if you have a zphot pdf close to zero.
    if z < 0:
        mask = samples[:, -1] > 0
        lumin_dist = lumin_dist[mask, :]

    lumin_dist_units = lumin_dist << u.keV / u.s

    return np.log10(lumin_dist_units.to(u.erg / u.s).value)
