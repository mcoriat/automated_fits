# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:45:53 2021

@author: alnoah
"""
import sherpa.astro.ui as shp

from utils import set_model_to_values


def set_full_model(full_model):
    shp.set_full_model("sim_src", full_model)


def get_number_of_channels(id):
    return len(shp.get_data(id).counts)


def get_full_model(id):
    return shp.get_model(id)


def get_src_model(full_model):
    srcm = full_model.parts[0]
    rmf_src = srcm.rmf
    arf_src = srcm.arf

    return srcm.parts[0].parts[1], rmf_src, arf_src

def get_src_model_noarf(full_model):
    srcm = full_model.parts[0]
    rmf_src = srcm.rmf

    return srcm.parts[0].parts[1], rmf_src


def get_bkg_particles_model(full_model):
    bkgm = full_model.parts[1]
    bkgm_particles = bkgm.parts[1].parts[0]

    # diagonals rmf and arf
    rmf_particles = bkgm_particles.rmf
    arf_particles = bkgm_particles.arf

    # Some times the background model raises to extremely high values in
    # the lower channels, which causes an error when generting poissonian counts
    # due to the high value of the expectation value. To avoid this, we set to
    # zero the arf values for the first ten channels.
    arf_particles.y[:10] = 0

    return bkgm_particles.parts[0].parts[1], rmf_particles, arf_particles

def get_bkg_particles_model_noarf(full_model):
    bkgm = full_model.parts[1]
    bkgm_particles = bkgm.parts[1].parts[0]

    # diagonals rmf 
    rmf_particles = bkgm_particles.rmf
    
    return bkgm_particles.parts[0].parts[1], rmf_particles


def get_bkg_astro_model(full_model):
    bkgm = full_model.parts[1]
    bkgm_astro = bkgm.parts[1].parts[1]
    rmf_astro = bkgm_astro.rmf
    arf_astro = bkgm_astro.arf

    return bkgm_astro.parts[0].parts[1], rmf_astro, arf_astro


def get_bkg_astro_model_noarf(full_model):
    bkgm = full_model.parts[1]
    bkgm_astro = bkgm.parts[1].parts[1]
    rmf_astro = bkgm_astro.rmf

    return bkgm_astro.parts[0].parts[1], rmf_astro


def get_src_parameters(id):
    return shp.get_exposure(id), shp.get_backscal(id)


def get_bkg_parameters(id):
    return shp.get_exposure(id, bkg_id=1), shp.get_backscal(id, bkg_id=1)


def delete_simulation():
    shp.delete_model('sim_bkg_particles')
    shp.delete_model('sim_bkg_astro')
    shp.delete_model('sim_src')

    shp.delete_data('sim_bkg_particles')
    shp.delete_data('sim_bkg_astro')
    shp.delete_data('sim_src')


def _fake_pha(id, model, rmf, arf, exposure, backscal, bkg=None):
    shp.set_source(id, model)
    shp.fake_pha(id, arf, rmf, exposure, backscal=backscal, bkg=bkg)


def _fake_pha_noarf(id, model, rmf, exposure, backscal, bkg=None):
    shp.set_source(id, model)
    shp.fake_pha(id, None, rmf, exposure, backscal=backscal, bkg=bkg)


def fake_bkg_particles(model, rmf, arf, exposure, backscal):
    _fake_pha("sim_bkg_particles", model, rmf, arf, exposure, backscal)


def fake_bkg_particles_noarf(model, rmf, exposure, backscal):
    _fake_pha_noarf("sim_bkg_particles", model, rmf, exposure, backscal)


def fake_astro_particles(model, rmf, arf, exposure, backscal):
    _fake_pha("sim_bkg_astro", model, rmf, arf, exposure, backscal)


def fake_astro_particles_noarf(model, rmf, exposure, backscal):
    _fake_pha_noarf("sim_bkg_astro", model, rmf, exposure, backscal)


def fake_src(model, rmf, arf, exposure, backscal, parameters, parameters_values):
    counts_particles = shp.get_data('sim_bkg_particles').counts
    counts_astro = shp.get_data('sim_bkg_astro').counts
    counts_bkg = counts_particles + counts_astro

    bkg = shp.get_data('sim_bkg_astro')
    bkg.counts = counts_bkg

    shp.set_source('sim_src', model)
    set_model_to_values(parameters, parameters_values)

    shp.fake_pha('sim_src', arf, rmf, exposure, bkg=bkg, backscal=backscal)
    counts_src = shp.get_data("sim_src").counts

    return counts_src
    

def fake_src_noarf(model, rmf, exposure, backscal, parameters, parameters_values):
    counts_particles = shp.get_data('sim_bkg_particles').counts
    counts_astro = shp.get_data('sim_bkg_astro').counts
    counts_bkg = counts_particles + counts_astro

    bkg = shp.get_data('sim_bkg_astro')
    bkg.counts = counts_bkg

    shp.set_source('sim_src', model)
    set_model_to_values(parameters, parameters_values)

    shp.fake_pha('sim_src', None, rmf, exposure, bkg=bkg, backscal=backscal)
    counts_src = shp.get_data("sim_src").counts

    return counts_src    
