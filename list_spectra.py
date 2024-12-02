import os

def list_spectra(srcid, srcid_obsid_dict, data_dir):
    list_spectra = []
    try:
        obsid_srcnum_list = srcid_obsid_dict[srcid]
    except KeyError:
        print(f"SRCID {srcid} not found in the dictionary.")
        return list_spectra
        
    for entry in obsid_srcnum_list:
        obsid = entry['OBS_ID']
        srcnum = entry['SRC_NUM']
        path = os.path.join(data_dir, obsid, 'pps')
        print(f"Checking directory: {path}")
        try:
            file_list = os.listdir(path)
        except FileNotFoundError:
            print(f"Directory {path} does not exist.")
            continue

        srchex = format(srcnum, '04X')
        pattern = f'SRSPEC{srchex}'
        print(f"Looking for pattern: {pattern}")
        for file in file_list:
            if pattern in file:
                list_spectra.append(os.path.join(path, file))
                
    print(f"Spectra found: {list_spectra}")
    return list_spectra


def test_list_spectra():
    srcid = 3067718060100029
    srcid_obsid_mapping = {
        srcid: [
            {"OBS_ID": "0677180601", "SRC_NUM": 38},
            {"OBS_ID": "0760940101", "SRC_NUM": 23}
        ]
    }
    data_dir = './test_data/'

    print("Starting test...")
    srcid_list_spectra = list_spectra(srcid, srcid_obsid_mapping, data_dir)
    nlist = len(srcid_list_spectra)
    expected_nlist = 3
    print(f"Spectra found: {srcid_list_spectra}")
    assert nlist == expected_nlist, f"Test failed: {nlist} != {expected_nlist}"
    print("Test passed.")


# Main block
if __name__ == "__main__":
    test_list_spectra()

