from astropy.io import fits


def read_stacked_catalog(catalog_file, srcid_ref):
    """
    Read a stacked catalog file and a SRCID and create a dictionary mapping each SRCID to its associated list of OBS_ID and SRC_NUM.

    Parameters:
    catalog_file (str): The path to the stacked catalog file.
    srcid_ref (long): the SRCID to be fitted

    Returns:
    dict: A dictionary associating the SRCID to its list of dictionaries containing OBS_ID and SRC_NUM.
    """
    with fits.open(catalog_file) as hdul:
        catalog_data = hdul[1].data

    # Create a dictionary to map each SRCID to its list of OBS_ID/SRC_NUM
    srcid_obsid_mapping = {}
    found = False

    for i in range(len(catalog_data)):
        srcid = catalog_data['SRCID'][i]
        obsid = catalog_data['OBS_ID'][i]
        srcnum = catalog_data['SRC_NUM'][i]

        if srcid == srcid_ref:
            if srcid in srcid_obsid_mapping:
                # Add to the list of dictionaries
                srcid_obsid_mapping[srcid].append({"OBS_ID": obsid, "SRC_NUM": srcnum})
            else:
                # Initialize with a list of dictionaries
                srcid_obsid_mapping[srcid] = []
                found = True
        elif found:
            break

    return srcid_obsid_mapping


def test_read_stacked_catalog():
    infile = './test_data/test_catalogue.fits'

    test_cases = [
        (1000000000000000, -1),
        (3072415020100239, 0),
        (3040339010100035, 1),
        (3030408050100122, 5)
    ]

    for srcid, expected in test_cases:
        try:
            result_dict = read_stacked_catalog(infile, srcid)

            # Calculate the actual result
            if len(result_dict) == 0:
                actual_result = -1
            else:
                actual_result = len(result_dict[srcid])

            # Assert the result
            assert actual_result == expected, f"Test failed for SRCID {srcid}: {actual_result} != {expected}"

            # Print success message
            print(f"Test passed for SRCID {srcid}. Expected: {expected}, Got: {actual_result}")

        except AssertionError as e:
            print(e)

            
if __name__ == "__main__":
    test_read_stacked_catalog()


