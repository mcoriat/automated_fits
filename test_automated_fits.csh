#
# 241111 FJC: testing automated_fits.py
#

#
# inexistent srcid: the program terminates with an error OK
rehash ; python -u automated_fits.py 1 ./test_data . ./test_data/RESPONSES  ./test_data/tests/ ./test_data/test_catalogue.fits ./test_data/tests.log --init --combine --fit_bkg --get_bkg_stat --fit_pl --redshift=1.0 --overwrite=1

#
# this srcid exists: it should return three OBS_ID, but only one with spectra: 1 pn and 2 MOS
#
# 241112 issues:
#    --use_tbabs_table=1 produces an error: Error: setPars() got an unexpected keyword argument 'tbabs'
#    The command just below produces an error: Error: module 'bxa.xspec' has no attribute 'Fit'
#rehash ; python -u automated_fits.py 3067718060100029 ./test_data . ./test_data/RESPONSES  ./test_data/tests/ ./test_data/test_catalogue.fits ./test_data/tests.log --fit_pl --redshift=1.0 --overwrite=1
#
rehash ; python3 -u automated_fits.py 3067718060100029 ./test_data . ./test_data/RESPONSES  ./test_data/tests/ ./test_data/test_catalogue.fits ./test_data/tests.log --fit_pl --redshift=1.0 --overwrite=1




#
exit
