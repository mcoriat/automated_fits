import os
import logging
import numpy as np
from astropy.io import fits
from list_spectra import list_spectra
from check_spectra import check_spectra
from merge_spectra import merge_spectra
from spectral_fitting import perform_spectrum_fitting
from xmm_backgrounds import get_pn_bkg_model_cached, get_mos_bkg_model_cached
import argparse

def save_fits_results(output_dir, srcid, model_name, pn_bkg_model, mos_bkg_model, source_results):
    """
    Saves background and source fitting results in a single FITS file.
    """
    filename = f"{srcid}_{model_name}.fits"
    filepath = os.path.join(output_dir, filename)

    def extract_model_params(model_dict, prefix):
        model = model_dict.get("model")
        params = {}
        if model:
            for p in model.parameters:
                params[f"{prefix}_{p.name}"] = p.values[0]
        params[f"{prefix}_fit_statistic"] = model_dict.get("fit_statistic", -1.0)
        params[f"{prefix}_fit_dof"] = model_dict.get("fit_dof", -1)
        return params

    pn_params = extract_model_params(pn_bkg_model, "pn")
    mos_params = extract_model_params(mos_bkg_model, "mos")

    all_data = {**pn_params, **mos_params, **source_results}
    cols = [fits.Column(name=k, format='E', array=np.array([v])) for k, v in all_data.items()]

    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto(filepath, overwrite=True)
    logging.info(f"Results saved to {filepath}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("srcid", help="Source ID")
    parser.add_argument("data_dir", help="Data directory")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("log_file", help="Log file")
    parser.add_argument("responses_dir", help="Responses directory")
    parser.add_argument("tests_dir", help="Tests directory")
    parser.add_argument("catalogue_file", help="Path to test catalogue file")
    parser.add_argument("--fit_pl", action="store_true", help="Fit powerlaw model")
    parser.add_argument("--fit_bb", action="store_true", help="Fit blackbody model")
    parser.add_argument("--fit_bkg", action="store_true", help="Fit background models")
    parser.add_argument("--redshift", type=float, help="Redshift value", default=1.0)
    parser.add_argument("--overwrite", type=int, default=1, help="Overwrite existing fits")
    parser.add_argument("--fit_apec_single", action="store_true", help="Fit APEC single model")
    parser.add_argument("--fit_apec_apec", action="store_true", help="Fit double APEC model")
    parser.add_argument("--fit_apec_apec_const", action="store_true", help="Fit double APEC model with constant")
    parser.add_argument("--fit_bremss", action="store_true", help="Fit Bremsstrahlung model")
    parser.add_argument("--fit_bbpl", action="store_true", help="Fit blackbody + powerlaw model")
    parser.add_argument("--fit_bbpl_const", action="store_true", help="Fit blackbody + powerlaw with constant model")
    parser.add_argument("--fit_bbpl_const2", action="store_true", help="Fit blackbody + powerlaw with second constant model")
    parser.add_argument("--fit_zpl", action="store_true", help="Fit redshifted powerlaw model")
    parser.add_argument("--fit_zplpl", action="store_true", help="Fit double redshifted powerlaw model")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(args.log_file)]
    )
    logger = logging.getLogger(__name__)

    logger.info("Starting automated fits script...")

    logger.info("Listing spectra...")
    srcid_obsid_mapping = {
        args.srcid: [
            ("0760940101", 11)
        ]
    }
    spectra_details = list_spectra(args.srcid, srcid_obsid_mapping, args.data_dir)
    print(f"Spectra details: {spectra_details}") 
    if not spectra_details:
        logger.error("No spectra found for the given SRCID.")
        return

    logger.info("Checking spectra...")
    wrapped_spectra = [{"spectrum_file": path, "background_file": None} for path in spectra_details]
    pn_spectra, mos_spectra = check_spectra(wrapped_spectra, args.log_file)

    logger.info("Merging spectra...")
    merged_spectra = merge_spectra(pn_spectra, mos_spectra, args.srcid, args.output_dir, args.log_file)
    logger.info(f"Merged spectra: {merged_spectra}")

    logger.info("Fitting background models...")
    pn_bkg_model, mos_bkg_model = {}, {}

    if args.fit_bkg:
        try:
            pn_path = next((s["spectrum_file"] for s in wrapped_spectra if "PN" in s["spectrum_file"].upper()), None)
            if pn_path:
                pn_bkg_model = get_pn_bkg_model_cached(pn_path, "phabs")
                stat = pn_bkg_model.get("fit_statistic", -1.0)
                dof = pn_bkg_model.get("fit_dof", -1)
                logger.info(f"PN background fit statistic: {stat}, dof: {dof}")
        except Exception as e:
            logger.error(f"PN background fitting failed: {e}")
            pn_bkg_model = {}

        try:
            mos_path = next((s["spectrum_file"] for s in wrapped_spectra if "M1" in s["spectrum_file"] or "M2" in s["spectrum_file"]), None)
            if mos_path:
                mos_bkg_model = get_mos_bkg_model_cached(mos_path, "phabs")
                stat = mos_bkg_model.get("fit_statistic", -1.0)
                dof = mos_bkg_model.get("fit_dof", -1)
                logger.info(f"MOS background fit statistic: {stat}, dof: {dof}")
        except Exception as e:
            logger.error(f"MOS background fitting failed: {e}")
            mos_bkg_model = {}

    logger.info("Processing merged spectra...")

    if len(merged_spectra) == 0:
        logger.error("No spectra available for processing.")
        return

    if args.fit_pl or args.fit_bb:
        logger.info("Performing spectrum fitting...")
        args_dict = vars(args)
        source_results = perform_spectrum_fitting(args.srcid, pn_spectra, mos_spectra, args_dict)

        print(f"PN Background Model: {pn_bkg_model}")
        print(f"MOS Background Model: {mos_bkg_model}")
        print(f"Source Results: {source_results}")

        save_fits_results(args.output_dir, args.srcid, "source_fitting", pn_bkg_model, mos_bkg_model, source_results)
    else:
        logger.warning("No source model selected for fitting.")

    logger.info("Automated fits script completed.")


if __name__ == "__main__":
    main()

