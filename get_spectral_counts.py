
import os
import numpy as np
from astropy.io import fits
import logging


logger = logging.getLogger(__name__)


# Function to read in a spectrum and its corresponding background file,
#     and return the counts and a flag if any incidence occured
def get_spectral_counts(infile,log_file,background_file=''):
    """
    Reads in input spectral file infile in FITS format, gets the total counts, and the
        header keywords BACKFILE, BACKSCAL and EXPOSURE.
    Then reads in the background file, gets its total counts and the header
        keyword BACKSCAL
        The name of the background file used depends on the value of background_file
            '' : default, gets the name diretly from BACKFILE in the header
            'pps' : gets the name replacing SRSPEC by BGSPEC in infile, as in the XMM-Newton PPS
            'path' : gets the name prepending to BACKFILE the full path to infile
            <filename>: if none of the others apply, using directly the filename
        
    Then calculates the net counts scaling with the backscale values
    Finally, writes out a tuple with this information and a flag
        The meaning of the flag is
            0: all correct
            -2: could not open source spectrum
            -1: cound not open background spectrum
            1: background counts <=0
            2: net counts <=0
    
        
    Parameters:
    - infile: input spectrum file in FITS format
    - log_file (str): The log file to write the messages.
    - spectrum_file (str): The name of the background file (see above)
    Returns:
    - spec_tuple (tuple): A tuple containing the name of the input spectrum, total source counts, total background counts, total net counts, exposure time, and the flag

    """
    
    
    # initialising output values
    flag=-2
    sp_counts=np.nan
    bg_counts=np.nan
    sp_netcts=np.nan
    sp_exp=np.nan
    sp_snr=np.nan
    
    # trying to open the input file
    try:
        hdul=fits.open(infile)
        sp_counts=hdul[1].data['COUNTS'].sum()
        sp_exp=hdul[1].header['EXPOSURE']
        sp_backscal=hdul[1].header['BACKSCAL']
        sp_exp=hdul[1].header['EXPOSURE']
        # background file name in the input file
        bgd_file=hdul[1].header['backfile']
        #
        hdul.close()
        del hdul
        #
        #
        # background file for input spectrum
        if (background_file==''):
            # getting the background file name from the header
            background_file=bgd_file
        elif (background_file=='pps'):
            # getting the background file name from the spectrum filename
            #    according to the rules for the XMM-Newton PPS
            background_file=infile.replace("SRSPEC","BGSPEC")
        elif (background_file=='path'):
            fullname=os.path.abspath(infile)
            fullpath=os.path.dirname(fullname)
            background_file=os.path.join(fullpath,bgd_file)
        else:
            # using the input value
            pass
        #
        try:
            bg_hdul=fits.open(background_file)
            bg_counts = bg_hdul[1].data['COUNTS'].sum()
            bg_backscal=bg_hdul[1].header['BACKSCAL']
            bg_hdul.close()
            del bg_hdul
            #
            sp_netcts=sp_counts-bg_counts*sp_backscal/bg_backscal
            #
            if (bg_counts<=0):
                flag=1
            elif (sp_netcts<=0 or sp_counts<=0):
                flag=2
            else:
                flag=0
            #
            # only calculating the signal-to-noise ratio if flag=0
            #      in part, to avoid square roots of negative numbers
            if (flag==0):sp_snr=sp_netcts/np.sqrt(2*sp_counts-sp_netcts)
        except:
            message=f'    Cannot open background FITS file {background_file}'
            logger.info(message)
            flag=-1
        #        
    except:
        message=f'    Cannot open spectrum FITS file {infile}'
        logger.info(message)
        flag=-2
    #
    return (infile,sp_counts,bg_counts,sp_netcts,sp_exp,flag,sp_snr)


def test_get_spectral_counts():
    output_dir = './test_data/test'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Set up logging
    log_file = os.path.join(output_dir, 'test_get_spectral_counts.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)

    def print_result(title, spec_tuple):
        print(f"\n🔹 {title}")
        print(f"""  → Output:
    Spectrum file:       {spec_tuple[0]}
    Source counts:       {spec_tuple[1]}
    Background counts:   {spec_tuple[2]}
    Net counts:          {spec_tuple[3]:.2f}
    Exposure time (s):   {spec_tuple[4]:.2f}
    Flag:                {spec_tuple[5]}
    SNR:                 {spec_tuple[6]:.2f}
""")

    # === Test 1: Nonexistent file
    infile = 'dummy_file_FJC.txt'
    spec_tuple = get_spectral_counts(infile, log_file)
    print_result("Test 1: Inexistent spectrum file", spec_tuple)
    assert spec_tuple[5] == -2

    # === Test 2: Valid PN spectrum — background via 'pps'
    infile = './test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ'
    spec_tuple = get_spectral_counts(infile, log_file, background_file='pps')
    print_result("Test 2: Valid PN spectrum — background via 'pps'", spec_tuple)
    assert spec_tuple[1] == 766
    assert spec_tuple[2] == 8474
    assert abs(spec_tuple[3] - 271.88) <= 0.01
    assert abs(spec_tuple[4] - 82181.94) <= 0.01
    assert spec_tuple[5] == 0

    # === Test 3: Valid PN spectrum — explicit background file
    background_file = './test_data/0760940101/pps/P0760940101PNS003BGSPEC0017.FTZ'
    spec_tuple = get_spectral_counts(infile, log_file, background_file=background_file)
    print_result("Test 3: Valid PN spectrum — background passed explicitly", spec_tuple)
    assert spec_tuple[1] == 766
    assert spec_tuple[2] == 8474
    assert abs(spec_tuple[3] - 271.88) <= 0.01
    assert abs(spec_tuple[4] - 82181.94) <= 0.01
    assert spec_tuple[5] == 0

    # === Test 4: Valid PN spectrum — background from header using 'path'
    spec_tuple = get_spectral_counts(infile, log_file, background_file='path')
    print_result("Test 4: Valid PN spectrum — background via header + path", spec_tuple)
    assert spec_tuple[1] == 766
    assert spec_tuple[2] == 8474
    assert abs(spec_tuple[3] - 271.88) <= 0.01
    assert abs(spec_tuple[4] - 82181.94) <= 0.01
    assert spec_tuple[5] == 0

    print("\n✅ All get_spectral_counts tests passed.")


    
if __name__ == "__main__":
    test_get_spectral_counts()

