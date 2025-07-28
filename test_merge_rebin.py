import os
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

    srcid = "3067718060100029"

    pn_list = [
        "./test_data/test_check_chain/P0760940101PNS003SRSPEC0017.FTZ"
    ]
    mos_list = [
        "./test_data/test_check_chain/P0760940101M2S002SRSPEC0017.FTZ"
    ]

    print("\n🔧 Running merge_spectra.merge_spectra() test...")
    merged_results = merge_spectra.merge_spectra(
        pn_list, mos_list, srcid, test_output_dir, merge_log_file, mincts=1
    )

    for item in merged_results:
        print("\n✅ Merged Result:")
        print(f"  File:               {item[0]}")
        print(f"  Source counts:      {item[1]}")
        print(f"  Background counts:  {item[2]}")
        print(f"  Net counts:         {item[3]:.2f}")
        print(f"  Exposure (s):       {item[4]:.2f}")
        print(f"  Flag:               {item[5]}")
        print(f"  SNR:                {item[6]:.2f}")
        print(f"  Instrument:         {item[7]}")

    print("\n🔧 Running rebin_spectrum.test_rebin_spectrum()...")
    rebin_spectrum.test_rebin_spectrum()


if __name__ == "__main__":
    run_merge_and_rebin_tests()

