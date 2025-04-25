import os
from pathlib import Path

def list_spectra(srcid, srcid_obsid_dict, data_dir):
    """
    Return paths to real spectrum files (SRSPEC...) associated with a given SRCID.

    Parameters:
    - srcid (str): Source ID as string (to match input arguments)
    - srcid_obsid_dict (dict): Dictionary mapping SRCID to list of (OBS_ID, SRC_NUM) tuples
    - data_dir (str): Base data directory

    Returns:
    - List[str]: Full paths to spectrum files that exist
    """
    result = []
    try:
        pairs = srcid_obsid_dict[srcid]
    except KeyError:
        return result

    for obsid, srcnum in pairs:
        spec_dir = Path(data_dir) / obsid / "pps"
        if not spec_dir.exists():
            continue
        files = os.listdir(spec_dir)
        hex_src = f"{srcnum:04X}"
        for fname in files:
            if f"SRSPEC{hex_src}" in fname:
                fpath = spec_dir / fname
                if fpath.exists():
                    result.append(str(fpath))
    return result

