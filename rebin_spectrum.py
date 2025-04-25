import os
import numpy as np
from astropy.io import fits
import logging

logger = logging.getLogger(__name__)

def get_spectral_counts(infile, log_file, background_file=''):
    spec_dict = {
        "spectrum_file": infile,
        "sp_counts": np.nan,
        "bg_counts": np.nan,
        "sp_netcts": np.nan,
        "sp_exp": np.nan,
        "flag": -2,
        "snr": np.nan
    }

    logger.info(f"Processing spectrum file: {infile}")

    try:
        with fits.open(infile) as hdul:
            logger.info(f"Opened spectrum file: {infile}")
            spec_dict["sp_counts"] = hdul[1].data['COUNTS'].sum()
            spec_dict["sp_exp"] = hdul[1].header['EXPOSURE']
            sp_backscal = hdul[1].header.get('BACKSCAL', 1.0)
            bgd_file = hdul[1].header.get('BACKFILE', '')

            if background_file == '' and bgd_file:
                background_file = os.path.join(os.path.dirname(infile), bgd_file)

            logger.info(f"Resolved background file: {background_file}")

            if not background_file or not os.path.exists(background_file):
                logger.warning(f"Background file not found: {background_file}, continuing without it.")
                spec_dict["sp_netcts"] = spec_dict["sp_counts"]
                spec_dict["flag"] = 0
                spec_dict["snr"] = np.sqrt(spec_dict["sp_counts"])
                return spec_dict

            try:
                with fits.open(background_file) as bg_hdul:
                    logger.info(f"Opened background file: {background_file}")
                    spec_dict["bg_counts"] = bg_hdul[1].data['COUNTS'].sum()
                    bg_backscal = bg_hdul[1].header.get('BACKSCAL', 1.0)
                    spec_dict["sp_netcts"] = spec_dict["sp_counts"] - (spec_dict["bg_counts"] * sp_backscal / bg_backscal)

                    if spec_dict["bg_counts"] <= 0:
                        spec_dict["flag"] = 1
                    elif spec_dict["sp_netcts"] <= 0 or spec_dict["sp_counts"] <= 0:
                        spec_dict["flag"] = 2
                    else:
                        spec_dict["flag"] = 0

                    if spec_dict["flag"] == 0:
                        spec_dict["snr"] = spec_dict["sp_netcts"] / np.sqrt(2 * spec_dict["sp_counts"] - spec_dict["sp_netcts"])

            except Exception as e:
                logger.warning(f"Failed to open or process background file {background_file}: {e}, continuing with source only.")
                spec_dict["sp_netcts"] = spec_dict["sp_counts"]
                spec_dict["flag"] = 0
                spec_dict["snr"] = np.sqrt(spec_dict["sp_counts"])

    except Exception as e:
        logger.error(f"Failed to open or process spectrum file {infile}: {e}")
        spec_dict["flag"] = -2

    return spec_dict

def rebin_spectrum(infile, outfile, background_file, log_file, mincts=1):
    logger.info(f"Rebinning file {infile} with a minimum of {mincts} counts per bin")
    outfile = os.path.abspath(outfile)

    spec_dict = get_spectral_counts(infile, log_file, background_file)
    logger.info(f"Returned spec_dict in rebin_spectrum: {spec_dict}")

    if spec_dict["flag"] < 0:
        logger.error(f"Invalid spec_dict received for {infile}: {spec_dict}")
        return spec_dict

    try:
        with fits.open(infile) as hdul:
            data = hdul[1].data
            counts = data['COUNTS']
            quality = data['QUALITY']
            grouping = data['GROUPING']

            rebinned_counts = []
            rebinned_quality = []
            rebinned_grouping = []

            current_bin = []
            for i, count in enumerate(counts):
                current_bin.append(count)
                if sum(current_bin) >= mincts:
                    rebinned_counts.append(sum(current_bin))
                    rebinned_quality.append(0)
                    rebinned_grouping.append(-1)
                    current_bin = []

            if current_bin:
                rebinned_counts.append(sum(current_bin))
                rebinned_quality.append(2)
                rebinned_grouping.append(-1)

            new_cols = fits.ColDefs([
                fits.Column(name='CHANNEL', format='I', array=np.arange(1, len(rebinned_counts) + 1)),
                fits.Column(name='COUNTS', format='J', array=np.array(rebinned_counts)),
                fits.Column(name='QUALITY', format='I', array=np.array(rebinned_quality)),
                fits.Column(name='GROUPING', format='I', array=np.array(rebinned_grouping))
            ])
            hdu = fits.BinTableHDU.from_columns(new_cols)
            hdu.writeto(outfile, overwrite=True)
            logger.info(f"Rebinned spectrum written to {outfile}")

    except Exception as e:
        logger.error(f"Failed to rebin spectrum file {infile}: {e}")
        spec_dict["flag"] = -1

    return spec_dict

