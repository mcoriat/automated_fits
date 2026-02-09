import os
import logging

logger = logging.getLogger(__name__)


# Function to check which spectra are suitable for fitting
def list_spectra(srcid, srcid_obsid_dict, data_dir,
                 dir_cache=None):
    """
    Find the spectra that actually exist associated to a given
    SRCID.

    Parameters:
    - srcid (long): the SRCID to be fitted
    - srcid_obsid_dict (dictionary): a dictionary of
        (OBS_ID,SRC_NUM) tuples associated to the SRCID
    - data_dir (str): Top directory where the spectra files
        are located.
    - dir_cache (dict or None): Optional cache of directory
        listings keyed by directory path. When provided,
        os.listdir() results are looked up here first and
        stored after the first real listing. Shared across
        sources in batch mode to avoid re-listing the same
        OBS_ID directories.
    Returns:
    - list_spectra (list): A list of strings containing the
        spectra that are present in data_dir
    """

    list_spectra = []
    try:
        obsid_srcnum_list = srcid_obsid_dict[srcid]
    except KeyError:
        return list_spectra
    #

    for obsid, srcnum in obsid_srcnum_list:
        # path corresponding to that obsid
        path = data_dir + '/' + obsid + '/pps/'
        # list of files in the directory (cached if possible)
        if dir_cache is not None:
            if path not in dir_cache:
                try:
                    dir_cache[path] = os.listdir(path=path)
                except FileNotFoundError:
                    logger.warning(
                        f"Directory {path} not found")
                    dir_cache[path] = []
            file_list = dir_cache[path]
        else:
            try:
                file_list = os.listdir(path=path)
            except FileNotFoundError:
                logger.warning(
                    f"Directory {path} not found")
                continue
        # converting srcnum to hex, with no prefix
        srchex = format(srcnum, '04X')
        # going through them to detect spectra for srcnum
        pattern = 'SRSPEC' + srchex
        for file in file_list:
            if pattern in file:
                list_spectra.append(path + file)
        #
    return list_spectra


def build_dir_listing_cache(catalog_mapping, data_dir):
    """
    Pre-scan all OBS_ID directories referenced by the catalog
    mapping and return a cache dict suitable for passing to
    list_spectra(dir_cache=...).

    This avoids repeated os.listdir() calls across sources
    that share the same observations (~14.6K observations
    for ~818K sources).

    Parameters:
    - catalog_mapping (dict): Output of
        read_stacked_catalog_batch(), i.e.
        {srcid: [(obsid, srcnum), ...], ...}
    - data_dir (str): Top directory with observation data.

    Returns:
    - dir_cache (dict): {dir_path: [filenames], ...}
    """
    # Collect unique OBS_IDs across all sources
    obsid_set = set()
    for pairs in catalog_mapping.values():
        for obsid, _srcnum in pairs:
            obsid_set.add(obsid)

    dir_cache = {}
    for obsid in sorted(obsid_set):
        path = data_dir + '/' + obsid + '/pps/'
        try:
            dir_cache[path] = os.listdir(path=path)
        except FileNotFoundError:
            logger.warning(
                f"Directory {path} not found (pre-cache)")
            dir_cache[path] = []

    logger.info(
        f"Pre-cached {len(dir_cache)} OBS_ID directories")
    return dir_cache

# function to test list_spectra
def test_list_spectra():
    srcid = 3067718060100029
    srcid_obsid_mapping = {
        srcid: [('0677180601', 38), ('0760940101', 23)]
    }
    data_dir = './test_data/'

    # Test without cache (original behaviour)
    result_no_cache = list_spectra(
        srcid, srcid_obsid_mapping, data_dir)
    print(f"Without cache: found {len(result_no_cache)} "
          f"spectra:")
    for s in result_no_cache:
        print(f"  - {s}")

    # Test with cache (batch mode)
    dir_cache = {}
    result_with_cache = list_spectra(
        srcid, srcid_obsid_mapping, data_dir,
        dir_cache=dir_cache)
    print(f"\nWith cache: found {len(result_with_cache)} "
          f"spectra:")
    for s in result_with_cache:
        print(f"  - {s}")
    print(f"Cache now has {len(dir_cache)} entries")

    # Results should be identical
    assert result_no_cache == result_with_cache, \
        "Cached and uncached results differ!"

    # Test build_dir_listing_cache
    cache2 = build_dir_listing_cache(
        srcid_obsid_mapping, data_dir)
    result_prebuilt = list_spectra(
        srcid, srcid_obsid_mapping, data_dir,
        dir_cache=cache2)
    assert result_no_cache == result_prebuilt, \
        "Pre-built cache results differ!"
    print("\nAll list_spectra tests passed.")


if __name__ == "__main__":
    test_list_spectra()
    
