import os



# Function to check which spectra are suitable for fitting
def list_spectra(srcid,srcid_obsid_dict,data_dir):
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
        return list_spectra
    #
        
    for obsid,srcnum in obsid_srcnum_list:
        # path corresponding to that obsid
        path=data_dir+'/'+obsid+'/pps/'
        # list of files in the directory
        file_list=os.listdir(path=path)
        # converting srcnum to hex, with no prefix, because spectra listed using it
        srchex=format(srcnum,'04X')
        # going through them to detect spectra asociated to srcnum
        pattern='SRSPEC'+srchex
        for file in file_list:
            if pattern in file:
                list_spectra.append(path+file)
        #            
    return list_spectra

# function to test list_spectra
def test_list_spectra():
    srcid = 3067718060100029
    srcid_obsid_mapping = {srcid: [('0677180601', 38), ('0760940101', 23)]}
    data_dir = './test_data/'
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, data_dir)
    print(f"Found {len(srcid_list_spectra)} spectra:")
    for s in srcid_list_spectra:
        print(f"  - {s}")

    
if __name__ == "__main__":
    test_list_spectra()
    
