import os
import datetime
from xspec import *
import bxa.xspec as bxa
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.table import Table, vstack
import glob

# python3 automated_fits.py 3067718060100029 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits dummy_output.txt --use_bxa --model_name=powerlaw --redshift=1.0 --overwrite=1 --export_results_fits --export_filename=fit_results.fits --bxa_output_dir=bxa_fit_results


    
def get_model_and_priors(model_name, redshift=0.0):
    if model_name == "powerlaw":
        model = Model("phabs*zpowerlw")
        model.zpowerlw.Redshift = redshift
        model.zpowerlw.Redshift.frozen = True

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.zpowerlw.PhoIndex.values = "2.0,,1.0,1.0,3.0,3.0"
        model.zpowerlw.norm.values = "1e-4,,1e-6,1e-6,1e-2,1e-2"

        priors = [
            bxa.create_uniform_prior_for(model, model.phabs.nH),
            bxa.create_uniform_prior_for(model, model.zpowerlw.PhoIndex),
            bxa.create_loguniform_prior_for(model, model.zpowerlw.norm)
        ]

    elif model_name == "apec_single":
        model = Model("phabs*apec")

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.apec.kT.values = "1.0,,0.1,0.1,10.0,10.0"
        model.apec.norm.values = "1e-4,,1e-6,1e-6,1e-2,1e-2"

        priors = [
            bxa.create_uniform_prior_for(model, model.phabs.nH),
            bxa.create_uniform_prior_for(model, model.apec.kT),
            bxa.create_loguniform_prior_for(model, model.apec.norm)
        ]

    elif model_name == "blackbody":
        model = Model("phabs*bbody")

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.bbody.kT.values = "0.1,,0.01,0.01,2.0,2.0"
        model.bbody.norm.values = "1e-4,,1e-6,1e-6,1e-2,1e-2"

        priors = [
            bxa.create_uniform_prior_for(model, model.phabs.nH),
            bxa.create_uniform_prior_for(model, model.bbody.kT),
            bxa.create_loguniform_prior_for(model, model.bbody.norm)
        ]

    elif model_name == "bremss":
        model = Model("phabs*bremss")

        model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
        model.bremss.kT.values = "5.0,,0.1,0.1,20.0,20.0"
        model.bremss.norm.values = "1e-4,,1e-6,1e-6,1e-2,1e-2"

        priors = [
            bxa.create_uniform_prior_for(model, model.phabs.nH),
            bxa.create_uniform_prior_for(model, model.bremss.kT),
            bxa.create_loguniform_prior_for(model, model.bremss.norm)
        ]

    else:
        raise ValueError(f"Unknown model: {model_name}")

    return model, priors

def fit_spectrum_bxa(spectrum_file, background_file, rmf_file, arf_file,
                     redshift=0.0, model_name="powerlaw",
                     output_base="bxa_fit_results", srcid="unknown"):
    AllData.clear()
    AllModels.clear()

    Fit.statMethod = "cstat"
    Plot.device = "/null"

    s = Spectrum(spectrum_file)
    s.background = background_file
    s.response = rmf_file
    s.response.arf = arf_file

    s.ignore("**-0.3 10.0-**")

    model, priors_list = get_model_and_priors(model_name, redshift)

    timestamp = datetime.datetime.now().strftime("%d%m%Y_%H%M")
    model_dirname = f"{model_name}_{timestamp}"
    output_dir = os.path.join(output_base, str(srcid), model_dirname)
    os.makedirs(output_dir, exist_ok=True)

    solver = bxa.BXASolver(transformations=priors_list,
                           outputfiles_basename=os.path.join(output_dir))
    solver.run(resume=False)

    chain_file = os.path.join(output_dir, "chain.fits")
    if os.path.exists(chain_file):
        with fits.open(chain_file) as hdul:
            samples = hdul[1].data
            samples_array = np.column_stack([samples[name] for name in samples.names])
            labels = solver.paramnames
            stds = np.std(samples_array, axis=0)
            valid_cols = stds > 0
            filtered_samples = samples_array[:, valid_cols]
            filtered_labels = [label for i, label in enumerate(labels) if valid_cols[i]]

            if filtered_samples.shape[1] > 0:
                fig = corner.corner(filtered_samples, labels=filtered_labels, show_titles=True, title_fmt=".3e")
                fig.savefig(os.path.join(output_dir, "corner.png"))

        posterior_median = np.median(samples_array, axis=0)

        # ===  Export results to FITS right after successful fit ===
        export_bxa_results_to_fits(
            srcid=srcid,
            output_base=output_base,
            fits_filename="fit_results.fits"
        )



        return {
            "parameter_names": labels,
            "posterior_median": posterior_median,
            "output_dir": output_dir
        }

    else:
        raise RuntimeError("chain.fits file not found after BXA run")
        

MODEL_SHORT_NAMES = {
    "powerlaw": "PL",
    "blackbody": "BB",
    "bremss": "BR",
    "apec_single": "AS"
}



def export_bxa_results_to_fits(srcid, output_base="bxa_fit_results", fits_filename="fit_results.fits"):
    fits_path = os.path.join(output_base, fits_filename)
    src_dir = os.path.join(output_base, str(srcid))

    os.makedirs(output_base, exist_ok=True)

    # Get all models in SRCID directory
    model_dirs = [d for d in os.listdir(src_dir) if os.path.isdir(os.path.join(src_dir, d))]
    short_map = {"powerlaw": "PL", "blackbody": "BB", "bremss": "BR", "apec_single": "AP"}
    model_data = {}

    for mdir in model_dirs:
        for long_name, short_name in short_map.items():
            if mdir.startswith(long_name):
                model_chain = os.path.join(src_dir, mdir, "chain.fits")
                if os.path.exists(model_chain):
                    with fits.open(model_chain) as hdul:
                        samples = hdul[1].data
                        param_names = hdul[1].columns.names
                        medians = [np.percentile(samples[p], 50) for p in param_names]
                        p16 = [np.percentile(samples[p], 16) for p in param_names]
                        p84 = [np.percentile(samples[p], 84) for p in param_names]
                        model_data[short_name] = {
                            "names": param_names,
                            "medians": medians,
                            "p16": p16,
                            "p84": p84
                        }

    # Create or update table without Models column
    if os.path.exists(fits_path):
        table = Table.read(fits_path)
        src_mask = table["SRCID"] == srcid
        if np.any(src_mask):
            idx = np.where(src_mask)[0][0]
            for model_short, pdata in model_data.items():
                for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                    col_med = f"{model_short}_{pname}_median"
                    col_p16 = f"{model_short}_{pname}_p16"
                    col_p84 = f"{model_short}_{pname}_p84"
                    for col, val in zip([col_med, col_p16, col_p84], [median, p16, p84]):
                        if col not in table.colnames:
                            table[col] = np.full(len(table), np.nan)
                        table[col][idx] = val
        else:
            new_row = {"SRCID": srcid}
            for model_short, pdata in model_data.items():
                for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                    new_row[f"{model_short}_{pname}_median"] = median
                    new_row[f"{model_short}_{pname}_p16"] = p16
                    new_row[f"{model_short}_{pname}_p84"] = p84
            for col in new_row.keys():
                if col not in table.colnames:
                    table[col] = np.full(len(table), np.nan)
            table.add_row(new_row)
        table.write(fits_path, overwrite=True)
    else:
        row_data = {"SRCID": [srcid]}
        for model_short, pdata in model_data.items():
            for pname, median, p16, p84 in zip(pdata["names"], pdata["medians"], pdata["p16"], pdata["p84"]):
                row_data[f"{model_short}_{pname}_median"] = [median]
                row_data[f"{model_short}_{pname}_p16"] = [p16]
                row_data[f"{model_short}_{pname}_p84"] = [p84]
        Table(row_data).write(fits_path, overwrite=True)

