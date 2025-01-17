import os

def list_spectra(srcid, srcid_obsid_dict, data_dir):
    """
    Lists the spectra and corresponding background files for a given SRCID.

    Parameters:
        srcid (int): Source ID to query.
        srcid_obsid_dict (dict): Mapping of SRCID to a list of dictionaries containing OBS_ID and SRC_NUM.
        data_dir (str): Base directory where data is stored.

    Returns:
        list: List of dictionaries containing spectrum and background file paths along with OBS_ID and SRC_NUM.
    """
    spectra_details = []
    try:
        obsid_srcnum_list = srcid_obsid_dict[srcid]
    except KeyError:
        print(f"SRCID {srcid} not found in the dictionary.")
        return spectra_details

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
                spectrum_path = os.path.join(path, file)
                background_path = spectrum_path.replace("SRSPEC", "BGSPEC")
                spectra_details.append({
                    "spectrum_file": spectrum_path,
                    "background_file": background_path,
                    "OBS_ID": obsid,
                    "SRC_NUM": srcnum
                })

    print(f"Spectra details: {spectra_details}")
    return spectra_details


def test_list_spectra():
    """
    Tests the list_spectra function with a predefined SRCID and mapping.

    Ensures that the function correctly identifies spectra and retains OBS_ID and SRC_NUM.
    """
    srcid = 3067718060100029
    srcid_obsid_mapping = {
        srcid: [
            {"OBS_ID": "0677180601", "SRC_NUM": 38},
            {"OBS_ID": "0760940101", "SRC_NUM": 23}
        ]
    }
    data_dir = './test_data/'

    print("Starting test...")
    spectra_details = list_spectra(srcid, srcid_obsid_mapping, data_dir)
    assert spectra_details, "No spectra details were found."
    for detail in spectra_details:
        assert "OBS_ID" in detail and "SRC_NUM" in detail, "OBS_ID or SRC_NUM is missing in the output."
    print(f"Test passed. Spectra details found: {spectra_details}")


if __name__ == "__main__":
    test_list_spectra()

