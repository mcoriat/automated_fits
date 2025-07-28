import os
import numpy as np
from astropy.io import fits
import logging

from get_spectral_counts import get_spectral_counts


logger = logging.getLogger(__name__)



# Function to check which spectra are suitable for fitting
def check_spectra(list_spectra, responses_dir, output_dir, log_file):
    '''
    Check the spectra for a given list of spectra. The good spectra are ungrouped and copied to output_dir
       together with soft links to their corresponding background, arf and rmf files
    
    If no pn spectra are selected, but at least one background spectrum has >0 counts or one source
        spectrum has >0 counts, one element is added to pn_spectra with the highest flag that applies to all spectra
    Similarly for MOS
    
    In addition to the two lists of tuples, it generates in output_dir symbolic links
        for all "good spectra" (flag=0) and their background, response and ancillary files
    
    
    
    Parameters:
    - list_spectra (str): List of spectra in data_dir
    - responses_dir (str): Top directory with the response matrices (rmf) for pn and MOS
    - output_dir (str): Directory to copy the good spectra, background, arf and rmf files to
    - log_file (str): The log file to write the messages.
    Returns:
    - pn_spectra (list): A list of tuples containing the full path and name for the good pn spectra, total source counts, total background counts, total net counts, exposure time, a flag, the signal-to-noise ratio, and a string with the instrument
    - mos_spectra (list): A list of tuples containing the full path and name for the good mos spectra, total source counts, total background counts, total net counts, exposure time, a flag, the signal-to-noise ratio, and a string with the instrument

    '''

    pn_spectra = []
    mos_spectra = []
    pn_flag=-1
    mos_flag=-1
    
    
    # Making output_dir and responses_dir absolute paths, as needed by symbolic links
    responses_dir=os.path.abspath(responses_dir)
    output_dir=os.path.abspath(output_dir)
        
    
    for spectrum_file in list_spectra:
        message=f'\n\n Working on file {spectrum_file}'
        logger.info(message)
        # getting absolute paths
        spectrum_file=os.path.abspath(spectrum_file)
        # Extracting just the spectrum name for output messages
        name=spectrum_file.split('/')[-1]
        banner=name
        #
        # finding out whether it is a pn or MOS spectrum
        pn= (banner.find("PN")>0)
        if (pn):
            instrument='pn'
        else:
            instrument='MOS'
        #
        # Initializing values for output tuple
        sp_counts=np.nan
        bg_counts=np.nan
        sp_netcts=np.nan
        sp_exp=np.nan
        sp_snr=np.nan
        # Initializing output tuples, only used when no spectra suitable for fitting are found
        pn_tuple=('',sp_counts,bg_counts,sp_netcts,sp_exp,-1,sp_snr,instrument)
        mos_tuple=('',sp_counts,bg_counts,sp_netcts,sp_exp,-1,sp_snr,instrument)
        # background filename corresponding to the source
        background_file=spectrum_file.replace("SRSPEC","BGSPEC")
        #
        # reading the input and background spectra, and returning the counts and a flag
        spec_tuple=get_spectral_counts(spectrum_file,log_file,background_file=background_file)
        #
        sp_counts=spec_tuple[1]
        bg_counts=spec_tuple[2]
        sp_netcts=spec_tuple[3]
        sp_exp=spec_tuple[4]
        sp_flag=spec_tuple[5]
        sp_snr=spec_tuple[6]
        # Check if the source spectrum exists
        if (sp_flag==-2):
            #    going for next spectrum if it does not
            message = f"       spectrum:{banner} - Missing source spectrum"
            logger.info(message)
            continue
        else:
            # source spectrum exists, what about the background spectrum?
            # Check if the background spectrum exists
            #
            if (sp_flag==-1):
                #     going for the next spectrum if it does not
                message = f"       spectrum:{banner} - Missing background spectrum {background_file}"
                logger.info(message)
                continue
            else:
                #
                # bits below left almost as they were before using get_spectral_counts
                #      not using sp_flag, left for future versions
                #
                if bg_counts > 0:
                    # Only proceed if background has counts
                    if sp_counts > 0:
                        # Only proceed if source spectrum has counts
                        #
                        # verifying now that the arf file exists
                        arf_file=spectrum_file.replace('SRSPEC','SRCARF')
                        if not os.path.exists(arf_file):
                            message = f"       spectrum:{banner} - Missing arf file {arf_file}"
                            logger.info(message)
                            continue
                        else:
                            # verifying now that the rmf file exists
                            hdul=fits.open(spectrum_file)
                            sp_header=hdul[1].header
                            hdul.close()
                            del hdul
                            response=sp_header['RESPFILE']
                            # full path to and name of rmf file depend on whetherpn/MOS
                            if (pn):
                                # adding _v19.0 at the end of the root
                                response19=response.split('.')[0]+'_v19.0.rmf'
                                rmf_file=os.path.join(responses_dir+'/PN',response19)
                            else:
                                # finding out the resolution of the rmf
                                specdelt=sp_header['SPECDELT']
                                rmf_file=os.path.join(responses_dir+'/MOS','{:d}eV/'.format(specdelt)+response)
                            #
                            
                            if not os.path.exists(rmf_file):
                                message = f"       spectrum:{banner} - Missing rmf file {rmf_file}"
                                logger.info(message)
                                continue
                            else:
                                # all necessary files present
                                
                                if (sp_netcts >0):
                                    # This spectrum is suitable for fitting
                                    message = f"       spectrum:{banner} - Good for fitting"
                                    logger.info(message)
                                    #
                                    spectrum_fit=spectrum_file.split('/')[-1]
                                    spectrum_fit=os.path.join(output_dir,spectrum_fit)
                                    #
                                    # appending spectrum to be fitted to output list
                                    out_tuple=(spectrum_fit,sp_counts,bg_counts,sp_netcts,sp_exp,0,sp_snr,instrument)
                                    message=f"        (fit file, source counts, background counts, net counts, exposure time, flag, SNR, instrument) = {out_tuple} "
                                    logger.info(message)
                                    if (pn):
                                        pn_spectra.append(out_tuple)
                                    else:
                                        mos_spectra.append(out_tuple)
                                    #
                                    # linking spectrum, background and arf in output dir
                                    for file in [spectrum_file, background_file, arf_file]:
                                        fitfile=file.split('/')[-1]
                                        fitfile=os.path.join(output_dir,fitfile)
                                        if(os.path.islink(fitfile)):os.unlink(fitfile)
                                        os.symlink(file,fitfile)
                                    # 
                                    # rmf file has to be treated differently because
                                    #     its name for pn is harwired in header and is named differently
                                    rmf_fit=os.path.join(output_dir,response)
                                    if(os.path.islink(rmf_fit)):os.unlink(rmf_fit)
                                    os.symlink(rmf_file,rmf_fit)
                                else:
                                    message = f"       spectrum:{banner} - Source spectrum has <=0 net counts"
                                    logger.info(message)
                                    if (pn):
                                        # always updating pn_tuple
                                        # spectra with sp_counts>0 but sp_netcts<0 take precedence in output
                                        # and, among them, the latest to be processed
                                        pn_flag=2
                                        pn_tuple=(spectrum_file,sp_counts,bg_counts,sp_netcts,sp_exp,pn_flag,sp_snr,instrument)
                                    else:
                                        # always updating mos_tuple
                                        # spectra with sp_counts>0 but sp_netcts<0 take precedence in output
                                        # and, among them, the latest to be processed
                                        mos_flag=2
                                        mos_tuple=(spectrum_file,sp_counts,bg_counts,sp_netcts,sp_exp,mos_flag,sp_snr,instrument)
                                    #spectrum
                                    
                                #
                            #
                        #
                    else:
                        message = f"       spectrum:{banner} - Source spectrum has <=0 total counts"
                        logger.info(message)
                        if (pn):
                            if (pn_flag<=1 or pn_flag[3]<=0):
                                # only updating pn_tuple if only spectra with <=0 bgd counts or
                                #      only spectra with <=0 total source counds found yet
                                # spectra with sp_counts>0 but sp_netcts<0 take precedence in output
                                pn_flag=2
                                pn_tuple=(spectrum_file,sp_counts,bg_counts,sp_netcts,sp_exp,pn_flag,sp_snr,instrument)
                        else:
                            if (mos_flag<=1 or mos_flag[3]<=0):
                                # only updating mos_tuple if only spectra with <=0 bgd counts or
                                #      only spectra with <=0 total source counds found yet
                                # spectra with sp_counts>0 but sp_netcts<0 take precedence in output
                                mos_flag=2
                                mos_tuple=(spectrum_file,sp_counts,bg_counts,sp_netcts,sp_exp,mos_flag,sp_snr,instrument)
                        #
                else:
                    # no background counts for this spectrum: updating the corresponding maximum
                    #    flag and going for the next spectrum
                    message = f"       spectrum:{banner} - Background spectrum has <=0 total counts"
                    logger.info(message)
                    if (pn):
                        if (pn_flag<=1):
                            # only updating pn_tuple if no pn spectra with >0 counts found yet
                            pn_flag=1
                            pn_tuple=(spectrum_file,sp_counts,bg_counts,sp_netcts,sp_exp,pn_flag,sp_snr,instrument)
                    else:
                        if (mos_flag<=1):
                            # only updating mos_tuple if no MOS spectra with >0 counts found yet
                            mos_flag=1
                            mos_tuple=(spectrum_file,sp_counts,bg_counts,sp_netcts,sp_exp,mos_flag,sp_snr,instrument)

                    #
                    continue
                #
            #
        #
        
    # if no pn spectra suitable for fitting, formatting output
    if (len(pn_spectra)==0 and pn_flag>0):
        # only enters here if no pn spectrum suitable for fitting is found
        #    but at least one spectrum 
        pn_spectra.append(pn_tuple)
    # otherwise, returning empty list

    # if no MOS spectra suitable for fitting, formatting output
    if (len(mos_spectra)==0 and mos_flag>0):
        mos_spectra.append(mos_tuple)
    # otherwise, returning empty list

                    
    return pn_spectra,mos_spectra


def test_check_spectra():
    output_dir = './test_data/test'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Set up logging
    log_file = os.path.join(output_dir, 'test_process_log.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)

    responses_dir = './test_data/RESPONSES'

    # === Test 1: PN only ===
    print("\n🔹 Test 1: PN only")
    list_spectra = ['./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ']
    pn_list, mos_list = check_spectra(list_spectra, responses_dir, output_dir, log_file)

    print(f"PN spectra found: {len(pn_list)}")
    for s in pn_list:
        print("  →", s)
    print(f"MOS spectra found: {len(mos_list)}")

    assert len(pn_list) == 1
    assert pn_list[0][5] == 0
    assert len(mos_list) == 0

    net = 271.88
    tot = 766
    snr = net / np.sqrt(2 * tot - net)
    assert abs(pn_list[0][6] - snr) <= 0.01

    # === Test 2: MOS only ===
    print("\n🔹 Test 2: MOS only")
    list_spectra = [
        './test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ',
        './test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ'
    ]
    pn_list, mos_list = check_spectra(list_spectra, responses_dir, output_dir, log_file)

    print(f"PN spectra found: {len(pn_list)}")
    print(f"MOS spectra found: {len(mos_list)}")
    for s in mos_list:
        print("  →", s)

    assert len(pn_list) == 0
    assert len(mos_list) == 2
    assert mos_list[0][5] == 0
    assert mos_list[1][5] == 0

    # === Test 3: All instruments ===
    print("\n🔹 Test 3: PN + MOS")
    list_spectra = [
        './test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ',
        './test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ',
        './test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ'
    ]
    pn_list, mos_list = check_spectra(list_spectra, responses_dir, output_dir, log_file)

    print(f"PN spectra found: {len(pn_list)}")
    for s in pn_list:
        print("  →", s)

    print(f"MOS spectra found: {len(mos_list)}")
    for s in mos_list:
        print("  →", s)

    assert len(pn_list) == 1
    assert pn_list[0][5] == 0
    assert len(mos_list) == 2
    assert mos_list[0][5] == 0
    assert mos_list[1][5] == 0


if __name__ == "__main__":
    test_check_spectra()

  
    
