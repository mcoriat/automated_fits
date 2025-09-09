import os
from astropy.io import fits
import logging


logger = logging.getLogger(__name__)


# Function to read the stacked 4XMM-DR11 catalog and map SRCID to corresponding OBSIDs
def read_stacked_catalog(catalog_file,srcid_ref, log_file):
    """
    Read a stacked catalog file and a SRCID and create a dictionary mapping each SRCID to its associated list of OBS_ID and SRC_NUM .

    Parameters:
    - catalog_file (str): The path to the stacked catalog file.
    - srcid_ref (long): the SRCID to be fitted
    - log_file (str): The log file to write the messages

    Returns:
    - dict: A dictionary associating the SRCID to its list of OBS_ID and SRC_NUM.
    
    """

    # Create a dictionary to map each SRCID to its list of OBS_ID/SRC_NUM
    srcid_obsid_mapping = {}

    try:
        with fits.open(catalog_file) as hdul:
            catalog_data = hdul[1].data
    except:
        # could not open the file, warning about it and returning
        logger.info('\n\n')
        logger.error(f'    Could not open file {catalog_file}\n\n')
        print(f'\n\n    ERROR: Could not open file {catalog_file}\n\n')
        return srcid_obsid_mapping
    #

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


def test_read_stacked_catalog():
    output_dir = './test_data/test'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Set up logging
    log_file = os.path.join(output_dir, 'test_process_log.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)

    # Test 1: inexistent file
    srcid=1
    infile='test_data/not_here.fits'
    #
    dic = read_stacked_catalog(infile, srcid, log_file)
    # it should capture the error and return an empty dictionary
    assert len(dic)==0
    

    # Other tests with different numbers of hits
    infile = './test_data/test_catalogue.fits'

    srcids = [
        1000000000000000,  # Not in file
        3072415020100239,  # Should return 0 hits
        3040339010100035,  # Should return 1 hit
        3030408050100122   # Should return 5 hits
    ]

    expected = [-1, 0, 1, 5]

    for i, srcid in enumerate(srcids):
        dic = read_stacked_catalog(infile, srcid, log_file)
        if len(dic) == 0:
            count = -1
        else:
            count = len(dic[srcid])
        print(f"SRCID: {srcid} → Found: {count}  | Expected: {expected[i]}")
        assert count == expected[i], f"Mismatch for SRCID {srcid}: got {count}, expected {expected[i]}"

          
if __name__ == "__main__":
    test_read_stacked_catalog()



          
