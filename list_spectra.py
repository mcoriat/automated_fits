import os
import logging


logger = logging.getLogger(__name__)



# Function to check which spectra are suitable for fitting
def list_spectra(srcid,srcid_obsid_dict,data_dir,log_file):
    """
    Find the spectra that actually exist associated to a given SRCID
    
    Parameters:
    - srcid (long): the SRCID to be fitted
    - srcid_obsid_dict (dictionary): a dictionary of (OBS_ID,SRC_NUM) tuples associated to the SRCID 
    - data_dir (str): Top directory where the spectra files are located.
    - log_file (str): The log file to write the messages.
    Returns:
    - list_spectra (list): A list of strings containing the spectra that are present in data_dir
    """

    list_spectra = []
    try:
        obsid_srcnum_list=srcid_obsid_dict[srcid]
    except KeyError:
        logger.warning('f     Could not find SRCID={srcid} in input dictionary')
        return list_spectra
    #
        
    for obsid,srcnum in obsid_srcnum_list:
        # path corresponding to that obsid
        path=data_dir+'/'+obsid+'/pps/'
        # list of files in the directory
        try:
            file_list=os.listdir(path=path)
            # converting srcnum to hex, with no prefix, because spectra listed using it
            srchex=format(srcnum,'04X')
            # going through them to detect spectra asociated to srcnum
            pattern='SRSPEC'+srchex
            for file in file_list:
                if pattern in file:
                    list_spectra.append(path+file)
        except:
            logger.warning(f'   Directory {path} not found, skipping OBS_ID={obsid} ')
            continue
        #            
    return list_spectra


# function to test list_spectra
def test_list_spectra():
    output_dir='./test_data/test'
    # absolute full path
    #out_dir=os.path.abspath(output_dir)
    data_dir = './test_data/'
    # creating output directory
    if not os.path.exists(output_dir) : os.mkdir(output_dir)
    # Set up logging
    log_file = os.path.join(output_dir, 'test_list_spectra.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode='w', force=True)

        
    # === Test 1: empty dictionary ===
    # it should never get to call list_spectra with an empty dictionary, but checking anyway
    print("\n🔹 Test 1: empty dictionary")
    #
    srcid = 0
    srcid_obsid_mapping = {}
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, data_dir, log_file)
    print(f"     Output list={srcid_list_spectra} ")
    # output list should be empty
    assert len(srcid_list_spectra)==0

        
    # === Test 2: dictionary with srcid, but empty list===
    # it should never get to call list_spectra with a dictionary with an empty list, but checking anyway
    print("\n🔹 Test 2: dictionary with srcid, but empty list")
    #
    srcid = 1
    srcid_obsid_mapping = {srcid:[]}
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, data_dir, log_file)
    print(f"     Output list={srcid_list_spectra} ")
    # output list should be empty
    assert len(srcid_list_spectra)==0

        
    # === Test 3: dictionary with srcid, and a list, but the O<BS_ID>/pps directory does not exist===
    print("\n🔹 Test 3: dictionary with srcid, and a list, but inexistent <OBS_ID>/pps directory")
    #
    srcid = 2
    srcid_obsid_mapping = {srcid:[('0',0)]}
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, data_dir, log_file)
    print(f"     Output list={srcid_list_spectra} ")
    # output list should be empty
    assert len(srcid_list_spectra)==0

        
    # === Test 4: dictionary with srcid, and two elements in the input list===
    print("\n🔹 Test 4: dictionary with srcid and two elements in input list")
    #
    srcid = 3067718060100029
    srcid_obsid_mapping = {srcid: [('0677180601', 38), ('0760940101', 23)]}
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, data_dir, log_file)
    print(f"Found {len(srcid_list_spectra)} spectra:")
    for s in srcid_list_spectra:
        print(f"  - {s}")
    assert len(srcid_list_spectra)==3




    
if __name__ == "__main__":
    test_list_spectra()
    
