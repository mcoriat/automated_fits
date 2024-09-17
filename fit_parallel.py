#!/usr/bin/env python3
import os
os.environ["OMP_NUM_THREADS"] = "1"

from pathlib import Path

from astropy.table import Table
from multiprocessing import Pool

import xmm_backgrounds


def get_cores():
    return 30


def wraper_fit(src):
    from fit import fit_spectra
    fit_spectra(src, model_name="background_only", scratch=True, verbose=False)


def src_iters():
    #sources_table_path = Path(".", "highz_xxl_candidates_81_tiles.fits")
    sources_table_path = Path(".", "highz_xxl_candidates_nodropouts_lowzdropouts_114_tiles.fits")
    sources_table = Table.read(sources_table_path)
    #mask_specz = sources_table["Zspec"] < 0

    return sources_table#[2:3]#[mask_specz]#[40:41]


if __name__ == "__main__":
    try:
        pool = Pool(get_cores())
        data_outputs = pool.map(wraper_fit, src_iters())

    except KeyboardInterrupt:
        # Close all running processes if CTRL+C
        pool.close()
        pool.join()

    finally:
        # To make sure processes are closed in the end, even if errors happen
        pool.close()
        pool.join()
