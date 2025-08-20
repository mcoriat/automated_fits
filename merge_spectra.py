import os
import shutil
import numpy as np
#from astropy.io import fits
import logging

from rebin_spectrum import rebin_spectrum


logger = logging.getLogger(__name__)



# Function to merge and rebin the spectra suitable for fitting
def merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1, test=False):
    """
    Given lists of pn and MOS spectra suitable for spectral fitting, the spectra are merged,
        and rebinned to >=1 count per bin.
    The merging is done per instrument, so at most one spectrum per each of pn and MOS is created,
        as well as some basic information (total, background and net counts, and exposure time)
    Not all spectra are merged, only the set that produces the highest cumulative signal-to-noise (SNR),
        going through the spectra in decreasing individual SNR

    The output is a list with two tuples, one per instrument, containing the spectrum name (with full path)
        and some basic information

    Parameters:
    - pn_spectra: list of tuples containing spectra suitable for fitting from pn
    - mos_spectra: list of tuples containing spectra suitable for fitting from MOS
    - srcid (long): SRCID of the source under study
    - output_dir (str): full path to the directory where the output files should be written
    - log_file (str): The log file to write the messages
    - mincts (int): minimum number per bin in the binned output merged spectrum. Optional, default is 1
    - test (boolean): if False, doing the actual merging of the selected spectra, if True, returning the spectrum with the highest
          individual SNR, for testing purposes. Optional, default is False
    Returns:
    - merged_spectra (list): A list of tuples containing the full path and name for the merged spectra, total source+background counts, total background counts, total net counts, exposure time, a flag, the signal-to-noise ratio, and a string stating whether it is a pn or MOS spectrum

    """

    message='\n\nMerging spectra'

    merged_spectra=[]

    for spec_list in [pn_spectra,mos_spectra]:
        nspec=len(spec_list)
        if (nspec==0):
            # no spectra suitable for fitting, no need to do anything
            continue
        elif (nspec==1 and spec_list[0][5]>0):
            # no spectra suitable for fitting, but flag=1 or 2 (negative or null spectrum or total or net counts)
            #   so just propagating the output basic information, and the instrument
            output_tuple=spec_list[0]
            merged_spectra.append(output_tuple)
            continue
        else:
            # there is at least one spectrum suitable for fitting
            # finding the instrument
            instrument=spec_list[0][7]
            message=f'\n\nWorking on instrument {instrument} with {nspec} spectra'
            logger.info(message)
             #
            # generating the merged spectrum filename
            outname='{}_{}.pha'.format(srcid,instrument)
            #
            merged_spectrum=os.path.join(output_dir,outname)
            #
            if(os.path.exists(merged_spectrum)):
                message=f' File {merged_spectrum} already exists in directory {output_dir}, it will be overwritten'
                logger.warning(message)
            #
            if (nspec>1):
                # merging the spectra here
                #
                # extracting from the list of tuples the individual snr and total and net counts for all spectra
                # in order to get the SNR we need to use tot and net, because bgd come from the background spectrum
                #     which is extracted from a different area
                specs=[val[0] for val in spec_list]
                snrs=np.array([val[6] for val in spec_list])
                tots=[val[1] for val in spec_list]
                nets=[val[3] for val in spec_list]
                # getting sorted indices
                sorted_indices = np.argsort(-snrs)
                print("   sorted_indices={} ".format(sorted_indices))
                #
                # calculating now the cumulative SNRs in decreasing order of individual SNR
                #
                # original index of highest individual SNR
                i0=sorted_indices[0]
                # cumulative values, initialized with the highest individual SNR
                cumsnrs=[snrs[i0]]
                cumtot=tots[i0]
                cumnet=nets[i0]
                # keeping the maximum cumulative SNR, so far
                max_cumsnr=snrs[i0]
                # sorted index of the maximum so far, the first in the list
                jmax=0
                # going through the full list of spectra to find the set that produces the maximum SNR
                j=0
                for i in sorted_indices[1:]:
                    j+=1
                    cumtot+=tots[i]
                    cumnet+=nets[i]
                    # the expresion below is equivalent to (tot-bgd)/sqrt(tot+bgd)
                    #     but, as discussed above, we cannot use directly bgd because it is normalised to a different area
                    cumsnr=(cumnet)/np.sqrt(2*cumtot-cumnet)
                    cumsnrs.append(cumsnr)
                    if(cumsnr>max_cumsnr):
                        # higher cumulative SNR, storing it
                        max_cumsnr=cumsnr
                        # its sorted index
                        jmax=j
                    #
                #
                # We now keep all spectra up to the maximum SNR found
                out_indices=sorted_indices[0:jmax+1]
                if (test):   
                   # just testing the algorithm above
                    # printing some results
                    message=f"Information about the merging of spectra for instrument {instrument}"
                    message+="   Original_index  Individual_SNR  Cumulative_SNR Spectrum"
                    for i in range(nspec):
                        j=sorted_indices[i]
                        line='\n  {:2d}  {:8.2f}  {:8.2f} {}'.format(j,snrs[j],cumsnrs[i],specs[j])
                        message+=line
                    logger.info(message)
                    print(message)
                    #
                    # generating a "merged" spectrum from the input spectrum with the highest SNR
                    #   just copying it
                    message=f"Test mode: Output spectrum is just the input spectrum with the highest SNR {spec_list[i0][0]}"
                    logger.warning(message)
                    shutil.copy2(spec_list[i0][0],merged_spectrum)
                    bkg_file = spec_list[i0][0].replace("SRSPEC", "BGSPEC")
                else:
                    # full merging needed, using a SAS script
                    # not implemented yet
                    continue
                #
            else:
                # only one input spectrum, no merging needed
                #
                # generating a "merged" spectrum from the only input spectrum
                #    by just copying it
                #
                shutil.copy2(spec_list[0][0],merged_spectrum)
                bkg_file = spec_list[0][0].replace("SRSPEC", "BGSPEC")
                #
            #
            # changing the extension of the merged spectrum to .grp for the binned spectrum
            binned_spectrum=os.path.splitext(merged_spectrum)[0]+'.grp'
            #
            # rebinning the spectrum, now with background file
            spec_tuple=rebin_spectrum(merged_spectrum,binned_spectrum,log_file,mincts=1, background_file=bkg_file)
            merged_spectra.append(spec_tuple)
        #
    #
    message='n\nFinished merging spectra. Output results:'
    for spec_tuple in merged_spectra:
        message=f"        (merged and binned file, source counts, background counts, net counts, exposure time, flag, signal-to-noise ratio, instrument) = {spec_tuple} "
        logger.info(message)                

    return merged_spectra



def test_merge_spectra():
    output_dir='./test_data/test'
    if not os.path.exists(output_dir) : os.mkdir(output_dir)
    # Set up logging
    log_file = os.path.join(output_dir, 'test_merge_spectra.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)
    # srcid
    srcid=3067718060100029
    #
    #
    # it should never get to call merge_spectra with two empty lists, but checking anyway
    pn_spectra=[]
    mos_spectra=[]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1)
    # output list should have no elements
    assert len(merged_list)==0

    #
    # One empty list and one with errors from before
    pn_spectra=[('dummyPN.fits',-1,1,-1,1000.0,2,-1,'pn')]
    mos_spectra=[]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1)
    # output list should have just one element, copied from pn_spectra
    assert len(merged_list)==1
    assert merged_list[0][5]==2

    #
    # only pn
    #  tuple just below from test_check_spectra.log
    pn_spectra=[('./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ', 766, 8474, 271.8810110703645, 82181.936317917, 0, 7.659018142527241,'pn')]
    mos_spectra=[]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1)
    # output list should have just one element, with the grouped version of the spectrum above
    assert len(merged_list)==1
    name=merged_list[0][0].split('/')[-1]
    assert name=='3067718060100029_pn.grp'

    #
    # only mos, test mode
    pn_spectra=[]
    mos_spectra=[('./test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ', 308, 14296, 63.03960341552707, 104469.411107063, 0, 2.6808126142724875,'MOS'),('./test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ', 236, 19138, 99.06503350999267, 105554.512163162, 0, 5.129840220900734,'MOS') ]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1, test=True)
    print("only MOS: merged_list ",merged_list)
    # output list should have just one element, with the grouped version of the second spectrum above
    assert len(merged_list)==1
    name=merged_list[0][0].split('/')[-1]
    assert name=='3067718060100029_MOS.grp'
    # signal-to-noise-ratio to check that chose the second MOS spectrum
    assert abs(mos_spectra[1][6]-merged_list[0][6])<=0.01

    #
    # everything, test mode
    pn_spectra=[('./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ', 766, 8474, 271.8810110703645, 82181.936317917, 0, 7.659018142527241,'pn')]
    mos_spectra=[('./test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ', 308, 14296, 63.03960341552707, 104469.411107063, 0, 2.6808126142724875,'MOS'),('./test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ', 236, 19138, 99.06503350999267, 105554.512163162, 0, 5.129840220900734,'MOS') ]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1, test=True)
    # output list should contain two elements, with the grouped versions of the first spectrum in the pn list and the second spectrum in the second list above
    assert len(merged_list)==2
    name=merged_list[0][0].split('/')[-1]
    assert name=='3067718060100029_pn.grp'
    name=merged_list[1][0].split('/')[-1]
    assert name=='3067718060100029_MOS.grp'
    # signal-to-noise-ratio to check that chose the second MOS spectrum
    assert abs(mos_spectra[1][6]-merged_list[1][6])<=0.01

    assert True
    
if __name__ == "__main__":
    test_merge_spectra()
    print("All merge and background tests completed successfully.")

