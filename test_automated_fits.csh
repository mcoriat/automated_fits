#
# Tests on automated_fits.py
#
##############
#
# the formal logs from within the program are sent to <output_dir>/<srcid>/<srcid>_process_log_<model name>.txt
#
###
#
# Cleaning the output directory by removing it completely
#
rm -Rf ./test_data/tests
#
### Test 1: inexistent srcid : OK
#           exit code 1 : OK
#      test_data/tests/1/1_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 1 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_1.txt --use_bxa --model_name=blackbody --redshift=1.0 --overwrite=1 --export_results_fits >& test_automated_fits_1.log &
#
### Test 2: existent srcid but no OBS_ID associated : OK
#           exit code 2 : OK
#      test_data/tests/3072415020100239/3072415020100239_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 3072415020100239 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_2.txt --use_bxa --model_name=blackbody --redshift=1.0 --overwrite=1 --export_results_fits >& test_automated_fits_2.log &
#
### Test 3: existent srcid with OBS_ID associated, but no spectrum on disk : OK
#           exit code 3 : OK
#      test_data/tests/3040339010100035/3040339010100035_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 3040339010100035 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_3.txt --use_bxa --model_name=blackbody --redshift=1.0 --overwrite=1 --export_results_fits >& test_automated_fits_3.log &
#
### Test 4: existent srcid with 1 OBS_ID associated 0700990101 with no pn and 1 MOS2 spectrum SRC_NUM=17 (hex 11)
#           exit code 0 : OK
#      test_data/tests/3070099010100059/3070099010100059_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 3070099010100059 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_4.txt --use_bxa --model_name=blackbody --redshift=0.0 --overwrite=1 --export_results_fits >& test_automated_fits_4.log &
#
### Test 5: existent srcid with 2 OBS_ID associated, 1 with no spectra and the other 0760940101 with no pn and MOS1 and MOS2 spectra SRC_NUM=81 (hex 51)
#                  only MOS2 used because merging with M1 would reduce the SNR
#           exit code 0 : OK
#      test_data/tests/3067718060100132/3067718060100132_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 3067718060100132 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_5.txt --use_bxa --model_name=blackbody --redshift=0.0 --overwrite=1 --export_results_fits >& test_automated_fits_5.log &
#
### Test 6: existent srcid with 2 OBS_ID associated, one with no spectra, and one 0760940101 with 1 pn and 2 MOS spectra SRC_NUM=23 (hex 17):
#                 the software has a feature, and it is that it fits separately the pn and merged-MOS spectra
#                    the values written to the output FITS file are those for the last spectrum fitted
#           exit code 0 : OK
#      test_data/tests/3067718060100029/3067718060100029_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 3067718060100029 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_6.txt --use_bxa --model_name=blackbody --redshift=1.0 --overwrite=1 --export_results_fits >& test_automated_fits_6.log &
#
### Test 7: existent srcid with 1 OBS_ID associated, 0760940101 with no pn and MOS1 and MOS2 spectra SRC_NUM=57 (hex 39)
#               MOS1 and MOS2 spectra merged successfully
#           exit code 0 : OK
#      test_data/tests/3067718060100101/3067718060100101_process_log_blackbody.txt
rehash ; python3 -u automated_fits.py 3067718060100101 ./test_data . ./test_data/RESPONSES ./test_data/tests ./test_data/test_catalogue.fits test_automated_fits_7.txt --use_bxa --model_name=blackbody --redshift=0.0 --overwrite=1 --export_results_fits >& test_automated_fits_7.log &




#
exit
