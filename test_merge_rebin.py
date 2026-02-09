import os
import logging
import rebin_spectrum
import merge_spectra


def run_merge_and_rebin_tests():
    test_output_dir = "./test_data/test_merge"
    os.makedirs(test_output_dir, exist_ok=True)

    # Ensure required log files exist
    merge_log_file = os.path.join(test_output_dir, "merge_test.log")
    rebin_log_file = os.path.join(test_output_dir, "rebin_test.log")
    open(merge_log_file, 'a').close()
    open(rebin_log_file, 'a').close()

    logging.basicConfig(
        filename=merge_log_file,
        level=logging.INFO,
        filemode='w',
        force=True
    )

    srcid = 3067718060100029

    # Create symlinks in test_output_dir for the source, background and arf files
    # (as check_spectra would normally do in the real pipeline)
    pn_src = os.path.abspath(
        "./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ"
    )
    pn_bkg = os.path.abspath(
        "./test_data/0760940101/pps/P0760940101PNS003BGSPEC0017.FTZ"
    )
    pn_arf = os.path.abspath(
        "./test_data/0760940101/pps/P0760940101PNS003SRCARF0017.FTZ"
    )
    mos_src = os.path.abspath(
        "./test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ"
    )
    mos_bkg = os.path.abspath(
        "./test_data/0760940101/pps/P0760940101M2S002BGSPEC0017.FTZ"
    )
    mos_arf = os.path.abspath(
        "./test_data/0760940101/pps/P0760940101M2S002SRCARF0017.FTZ"
    )

    for src_file in [pn_src, pn_bkg, pn_arf, mos_src, mos_bkg, mos_arf]:
        link = os.path.join(
            test_output_dir, os.path.basename(src_file)
        )
        if os.path.islink(link) or os.path.exists(link):
            os.remove(link)
        if os.path.exists(src_file):
            os.symlink(src_file, link)

    # Tuples as produced by check_spectra:
    # (spectrum_path, sp_counts, bg_counts, sp_netcts,
    #  sp_exp, flag, sp_snr, instrument)
    pn_spec_in_outdir = os.path.join(
        test_output_dir, "P0760940101PNS003SRSPEC0017.FTZ"
    )
    mos_spec_in_outdir = os.path.join(
        test_output_dir, "P0760940101M2S002SRSPEC0017.FTZ"
    )

    pn_list = [
        (pn_spec_in_outdir,
         766, 8474, 271.88, 82181.94,
         0, 7.659, 'pn')
    ]
    mos_list = [
        (mos_spec_in_outdir,
         236, 19138, 99.07, 105554.51,
         0, 5.130, 'MOS')
    ]

    print("\n Running merge_spectra.merge_spectra() test...")
    merged_results = merge_spectra.merge_spectra(
        pn_list, mos_list, srcid,
        test_output_dir, merge_log_file, mincts=1
    )

    for item in merged_results:
        print("\n Merged Result:")
        print(f"  File:               {item[0]}")
        print(f"  Source counts:      {item[1]}")
        print(f"  Background counts:  {item[2]}")
        print(f"  Net counts:         {item[3]:.2f}")
        print(f"  Exposure (s):       {item[4]:.2f}")
        print(f"  Flag:               {item[5]}")
        print(f"  SNR:                {item[6]:.2f}")
        print(f"  Instrument:         {item[7]}")
        if len(item) > 8:
            sp_dic = item[8]
            print(f"  sp_dic SPECFILE:    {sp_dic['SPECFILE']}")
            print(f"  sp_dic BACKFILE:    {sp_dic['BACKFILE']}")
            print(f"  sp_dic RESPFILE:    {sp_dic['RESPFILE']}")
            print(f"  sp_dic ANCRFILE:    {sp_dic['ANCRFILE']}")

    print("\n🔧 Running rebin_spectrum.test_rebin_spectrum()...")
    rebin_spectrum.test_rebin_spectrum()


if __name__ == "__main__":
    run_merge_and_rebin_tests()

