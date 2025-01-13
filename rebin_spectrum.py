import os
import numpy as np
from astropy.io import fits
import logging

from get_spectral_counts import get_spectral_counts

logger = logging.getLogger(__name__)

# Function to rebin a spectrum, writing to disk a copy of it and returning some basic info
def rebin_spectrum(infile, outfile, log_file, mincts=1):
    """
    Reads in an input file, rebins it to have >=mincts counts per bin, taking into account
    if there are bad channels at the beginning, and flagging as "bad" empty channels at the end.
    """
    message = f'\n\n Rebinning file {infile} with a minimum of {mincts} counts per bin'
    logger.info(message)

    outfile = os.path.abspath(outfile)

    # Getting the values for the output dictionary
    spec_info = list(get_spectral_counts(infile, log_file, background_file='path'))
    logger.info(f"spec_info from get_spectral_counts: {spec_info}")

    # Check for invalid data in spec_info
    if len(spec_info) < 7 or any(isinstance(val, str) or val is None for val in spec_info):
        logger.error(f"Invalid spec_info data: {spec_info}")
        return {
            "spectrum_file": outfile,
            "sp_counts": np.nan,
            "bg_counts": np.nan,
            "sp_netcts": np.nan,
            "sp_exp": np.nan,
            "flag": -2,
            "snr": np.nan,
        }

    # Populate the dictionary
    spec_dict = {
        "spectrum_file": outfile,
        "sp_counts": spec_info[1],
        "bg_counts": spec_info[2],
        "sp_netcts": spec_info[3],
        "sp_exp": spec_info[4],
        "flag": spec_info[5],
        "snr": spec_info[6],
    }

    # Checking if input file can be opened
    if spec_dict['flag'] == -2:
        logger.error(f"Cannot open FITS spectrum file {infile}")
        return spec_dict

    # (Rest of the code remains unchanged, handling FITS files and binning.)
    ...


    ########
    # Reading the input file header and data
    with fits.open(infile) as hdul:
        data1 = hdul[1].data
        header0 = hdul[0].header.copy()
        header1 = hdul[1].header.copy()

    counts = data1['COUNTS']
    nchan = len(counts)

    # Checking if quality info provided, otherwise fill with 0 (good quality)
    try:
        quality = data1['QUALITY']
        hasQuality = True
        nbad = 0
        for i in range(nchan):
            if quality[i] > 0:
                nbad += 1
            else:
                break
        message = f'    {nbad} initial channels are marked as bad'
        logger.info(message)

        badmask = quality > 0
        nbadmask = sum(badmask)
        message = f'    {nbadmask} total channels are marked as bad'
        logger.info(message)

        if nbadmask > nbad:
            message = '    Additional bad channels not at the beginning'
            logger.warning(message)
    except:
        hasQuality = False
        quality = np.full(nchan, 0)
        nbad = 0

    try:
        grouping = data1['GROUPING']
        hasGrouping = True
    except:
        hasGrouping = False
        grouping = np.full(nchan, 1)

    groupCounts = 0
    first = True
    first_i = 0

    for i in range(nchan):
        if i >= nbad:
            groupCounts += counts[i]
            if groupCounts >= mincts:
                if first:
                    pass
                else:
                    grouping[i] = -1
                groupCounts = 0
                first = True
                first_i = i + 1
            elif groupCounts < 1:
                if not first:
                    grouping[i] = -1
                first = False

    if groupCounts < mincts:
        grouping[first_i:] = 1
        quality[first_i:] = 2
        message = f'    Setting channels {first_i} to last as bad (quality=2)'
        logger.info(message)

    if hasGrouping:
        data1['GROUPING'] = grouping
        hdu = fits.BinTableHDU.from_columns(data1.columns, header=header1)
    else:
        colGrouping = fits.ColDefs([fits.Column(name='GROUPING', format='I', array=grouping)])
        hdu = fits.BinTableHDU.from_columns(data1.columns + colGrouping, header=header1)

    if hasQuality:
        hdu.data['QUALITY'] = quality
        hdu = fits.BinTableHDU.from_columns(hdu.data.columns, header=hdu.header)
    else:
        colQuality = fits.ColDefs([fits.Column(name='QUALITY', format='I', array=quality)])
        hdu = fits.BinTableHDU.from_columns(hdu.columns + colQuality, header=hdu.header)

    hdu.name = 'SPECTRUM'

    hdu0 = fits.PrimaryHDU(header=header0)
    hdu_new = fits.HDUList([hdu0, hdu])

    hdu_new.writeto(outfile, overwrite=True)
    message = f' {nchan} channels written out to file {outfile}'
    logger.info(message)

    del hdul, hdu, hdu_new, hdu0, grouping, quality, counts, data1, header0, header1

    return spec_dict


def test_rebin_spectrum():
    output_dir = './test_data/test'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    log_file = os.path.join(output_dir, 'test_rebin_spectrum.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)

    infile = 'dummy_file_FJC.txt'
    spec_dict = rebin_spectrum(infile, '', log_file)
    assert spec_dict['flag'] == -2

    mincts = 1
    infile = './test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ'
    outfile = os.path.join(output_dir, 'test_rebin_spectrum.grp1')
    spec_dict = rebin_spectrum(infile, outfile, log_file, mincts=mincts)
    assert spec_dict['flag'] == 0
    assert os.path.exists(outfile)

    with fits.open(outfile) as sp_hdul:
        counts = sp_hdul[1].data['COUNTS']
        quality = sp_hdul[1].data['QUALITY']
        grouping = sp_hdul[1].data['GROUPING']
        good = quality == 0
        tot_good_cts = sum(counts[good])
        grouping_good = grouping[good]
        tot_good_bins = sum(grouping_good == 1)
        assert tot_good_cts == 730
        assert tot_good_bins == 546
        good_cts = counts[good]
        good_grouping = grouping[good]
        ngood = len(good_cts)
        tot_cts = good_cts[0]

        for i in range(1, ngood):
            if good_grouping[i] == 1:
                assert tot_cts >= mincts
                tot_cts = good_cts[i]
            else:
                tot_cts += good_cts[i]

        assert tot_cts >= mincts

    net = 271.88
    tot = 766
    snr = net / np.sqrt(2 * tot - net)
    assert abs(spec_dict['snr'] - snr) <= 0.01

