import numpy as np
from astropy.io import fits


# Function to read the stacked 4XMM-DR11 catalog and map SRCID to corresponding OBSIDs
def read_stacked_catalog(catalog_file,srcid_ref):
    """
    Read a stacked catalog file and a SRCID and create a dictionary mapping each SRCID to its associated list of OBS_ID and SRC_NUM .

    Parameters:
    catalog_file (str): The path to the stacked catalog file.

    srcid_ref (long): the SRCID to be fitted

    Returns:
    dict: A dictionary associating the SRCID to its list of OBS_ID and SRC_NUM.
    """

    with fits.open(catalog_file) as hdul:
        catalog_data = hdul[1].data

    # Create a dictionary to map each SRCID to its list of OBS_ID/SRC_NUM
    srcid_obsid_mapping = {}
    # Flag to see if we have reached the SRCID yet
    found=False
    # loop over the rows in the input file
    for i in range(len(catalog_data)):
        srcid = catalog_data['SRCID'][i]
        obsid = catalog_data['OBS_ID'][i]
        srcnum=catalog_data['SRC_NUM'][i]

        # checking if this row corresponds to the input SRCID
        if srcid==srcid_ref:    
            if srcid in srcid_obsid_mapping:
                # second and consecutive rows appended
                srcid_obsid_mapping[srcid].append((obsid,srcnum))
            else:
                # ignoring the first row for each SRCID, because no OBS_ID on it
                # initializing the list of tuples (OBS_ID,SRC_NUM)
                srcid_obsid_mapping[srcid] = []
                # setting the flag
                found=True
        elif found:
            # all rows for the same SRCID are consecutive so, once SRCID has been found
            #     all following rows with different SRCID can be safely skipped
            break

    return srcid_obsid_mapping


def read_stacked_catalog_batch(catalog_file, srcid_list):
    """
    Read a stacked catalog file once and return mappings for
    all requested SRCIDs.  Uses numpy vectorized filtering
    instead of a Python for-loop, so the cost is ~30s total
    for the full 2.5M-row catalog rather than ~12s per source.

    The first row for each SRCID group is skipped (no OBS_ID),
    replicating the convention in read_stacked_catalog().

    Parameters:
    - catalog_file (str): Path to the stacked catalog FITS file.
    - srcid_list (list of int): List of SRCIDs to look up.

    Returns:
    - dict: {srcid: [(obsid, srcnum), ...], ...}
      SRCIDs not found in the catalog are omitted from the dict.
      SRCIDs whose only row is the header row (no OBS_ID pairs)
      map to an empty list.
    """
    with fits.open(catalog_file) as hdul:
        data = hdul[1].data
        all_srcids = data['SRCID']
        all_obsids = data['OBS_ID']
        all_srcnums = data['SRC_NUM']

    # Vectorized mask: keep only rows matching requested SRCIDs
    unique_requested = np.array(list(set(srcid_list)),
                                dtype=all_srcids.dtype)
    mask = np.isin(all_srcids, unique_requested)
    sel_srcids = all_srcids[mask]
    sel_obsids = all_obsids[mask]
    sel_srcnums = all_srcnums[mask]

    # Build the output dictionary
    result = {}
    prev_srcid = None
    for i in range(len(sel_srcids)):
        sid = int(sel_srcids[i])
        if sid != prev_srcid:
            # First row of a new SRCID group: skip it (no OBS_ID)
            result[sid] = []
            prev_srcid = sid
        else:
            result[sid].append(
                (str(sel_obsids[i]).strip(),
                 int(sel_srcnums[i]))
            )

    return result


def test_read_stacked_catalog():
    infile = './test_data/test_catalogue.fits'

    srcids = [
        1000000000000000,  # Not in file
        3072415020100239,  # Should return 0 hits
        3040339010100035,  # Should return 1 hit
        3030408050100122   # Should return 5 hits
    ]

    expected = [-1, 0, 1, 5]

    for i, srcid in enumerate(srcids):
        dic = read_stacked_catalog(infile, srcid)
        if len(dic) == 0:
            count = -1
        else:
            count = len(dic[srcid])
        print(f"SRCID: {srcid} → Found: {count}  | Expected: {expected[i]}")
        assert count == expected[i], f"Mismatch for SRCID {srcid}: got {count}, expected {expected[i]}"

def test_read_stacked_catalog_batch():
    infile = './test_data/test_catalogue.fits'

    srcid_list = [
        1000000000000000,  # Not in file
        3072415020100239,  # Should return 0 hits
        3040339010100035,  # Should return 1 hit
        3030408050100122   # Should return 5 hits
    ]

    expected = {
        1000000000000000: -1,  # not found at all
        3072415020100239: 0,
        3040339010100035: 1,
        3030408050100122: 5
    }

    batch_result = read_stacked_catalog_batch(infile, srcid_list)

    for srcid in srcid_list:
        if srcid not in batch_result:
            count = -1
        else:
            count = len(batch_result[srcid])
        exp = expected[srcid]
        print(f"SRCID: {srcid} → Found: {count}  "
              f"| Expected: {exp}")
        assert count == exp, (
            f"Batch mismatch for SRCID {srcid}: "
            f"got {count}, expected {exp}")

    # Cross-check: batch results must match single-source results
    for srcid in srcid_list:
        single = read_stacked_catalog(infile, srcid)
        if srcid not in batch_result:
            assert len(single) == 0
        else:
            if len(single) == 0:
                assert srcid not in batch_result or \
                    len(batch_result[srcid]) == 0
            else:
                assert batch_result[srcid] == single[srcid]

    print("\nAll batch catalog tests passed.")


if __name__ == "__main__":
    test_read_stacked_catalog()
    test_read_stacked_catalog_batch()



          
