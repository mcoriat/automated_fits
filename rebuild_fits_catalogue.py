#!/usr/bin/env python3
import os
import numpy as np
from astropy.table import Table, vstack, join
from astropy.io import fits

def process_chain(chain_file, model_name):
    """Extract median, 16th, and 84th percentiles for all parameters in chain.fits"""
    with fits.open(chain_file) as hdul:
        samples = hdul[1].data
        param_names = samples.names
        samples_array = np.column_stack([samples[name] for name in param_names])
        p50 = np.percentile(samples_array, 50, axis=0)
        p16 = np.percentile(samples_array, 16, axis=0)
        p84 = np.percentile(samples_array, 84, axis=0)

    results = {}
    for i, name in enumerate(param_names):
        short_name = f"{model_name}_{name}"
        results[f"{short_name}_50"] = p50[i]
        results[f"{short_name}_16"] = p16[i]
        results[f"{short_name}_84"] = p84[i]
    return results

def rebuild_catalogue(base_dir="bxa_fit_results", input_catalogue="./test_data/test_catalogue.fits", output_fits="fit_results.fits"):
    # Load input catalogue
    if not os.path.exists(input_catalogue):
        raise FileNotFoundError(f"Input catalogue {input_catalogue} not found")
    input_table = Table.read(input_catalogue)

    # Create a dict to accumulate results by SRCID
    results_dict = {str(srcid): {} for srcid in input_table["SRCID"]}

    # Walk through BXA results
    for srcid in os.listdir(base_dir):
        src_path = os.path.join(base_dir, srcid)
        if not os.path.isdir(src_path):
            continue
        
        for model_folder in os.listdir(src_path):
            model_path = os.path.join(src_path, model_folder)
            if not os.path.isdir(model_path):
                continue
            
            chain_file = os.path.join(model_path, "chain.fits")
            if os.path.exists(chain_file):
                model_name = model_folder.split("_")[0]
                model_results = process_chain(chain_file, model_name)
                results_dict[srcid].update(model_results)

    # Add results columns to input table
    for srcid, model_results in results_dict.items():
        idx = np.where(input_table["SRCID"] == int(srcid))[0]
        if idx.size > 0:
            for key, val in model_results.items():
                if key not in input_table.colnames:
                    input_table[key] = np.full(len(input_table), np.nan)
                input_table[key][idx] = val

    # Save updated FITS
    input_table.write(output_fits, overwrite=True)
    print(f"✅ FITS catalogue updated and saved to {output_fits}")

if __name__ == "__main__":
    rebuild_catalogue()

