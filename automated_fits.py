import sys
import os
from read_stacked_catalog import read_stacked_catalog
from list_spectra import get_spectra
from check_spectra import check_spectra
from spectral_fitting import run_spectral_fit


def main():
    srcid = sys.argv[1]
    datapath = sys.argv[2]
    outpath = sys.argv[3]
    logpath = sys.argv[4]
    response_path = sys.argv[5]
    temp_path = sys.argv[6]
    catalogue_file = sys.argv[7]

    extra_args = sys.argv[8:]
    redshift = None
    overwrite = False
    fit_pl = False
    fit_bkg = False

    for arg in extra_args:
        if arg.startswith("--redshift="):
            redshift = float(arg.split("=")[1])
        if arg == "--overwrite=1":
            overwrite = True
        if arg == "--fit_pl":
            fit_pl = True
        if arg == "--fit_bkg":
            fit_bkg = True

    print(f"\n==== Fitting SRCID {srcid} ====")

    catalogue = read_stacked_catalog(catalogue_file)
    print("Available SRCIDs in catalog:", list(catalogue.keys()))

    spectra = get_spectra(srcid, datapath, catalogue)
    print("Spectra details:", [s['spectrum_file'] for s in spectra if isinstance(s, dict) and 'spectrum_file' in s])

    results = []
    for spec in spectra:
        if not isinstance(spec, dict):
            print("WARNING: Skipping non-dict spectrum entry:", spec)
            continue

        result_dict = check_spectra(
            spec,
            datapath=datapath,
            outpath=outpath,
            logpath=logpath,
            response_path=response_path,
            temp_path=temp_path,
            overwrite=overwrite
        )

        if not result_dict:
            continue

        if isinstance(result_dict, dict):
            if result_dict.get("instrument") == "pn":
                print(f"-> PN spectrum detected for {srcid}")
            result = fit_spectra(
                result_dict,
                redshift=redshift,
                fit_pl=fit_pl,
                fit_bkg=fit_bkg,
                overwrite=overwrite
            )
            results.append(result)
        else:
            print(f"WARNING: Result was not a dictionary for {spec}")

if __name__ == "__main__":
    main()

