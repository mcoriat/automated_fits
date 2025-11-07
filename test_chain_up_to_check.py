from read_stacked_catalog import read_stacked_catalog
from list_spectra import list_spectra
from check_spectra import check_spectra

srcid = 3067718060100029
catalog_file = './test_data/test_catalogue.fits'
data_dir = './test_data'
responses_dir = './test_data/RESPONSES'
output_dir = './test_data/test_check_chain'
log_file = output_dir + '/log.txt'

import os
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Step 1: Catalog lookup
mapping = read_stacked_catalog(catalog_file, srcid)
print(f"\nSRCID {srcid} → {mapping[srcid]}")

# Step 2: Locate available spectra
spectra = list_spectra(srcid, mapping, data_dir)
print(f"\nSpectra found:")
for s in spectra:
    print("  →", s)

# Step 3: Check spectra and categorize
pn_list, mos_list = check_spectra(spectra, responses_dir, output_dir, log_file)

print("\nPN spectra suitable for fitting:")
for p in pn_list:
    print("  →", p)

print("\nMOS spectra suitable for fitting:")
for m in mos_list:
    print("  →", m)

