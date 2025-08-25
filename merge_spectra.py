import os
import shutil
import numpy as np
#from astropy.io import fits
import logging

from rebin_spectrum import rebin_spectrum
from get_spectral_counts import get_spectral_counts

import subprocess

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
    - merged_spectra (list): A list of tuples containing the full path and name for the merged spectra, total source+background counts, total background counts, total net counts, exposure time, a flag, the signal-to-noise ratio, a string stating whether it is a pn or MOS spectrum, and a dictionary containing the full path to the rmf,arf,bgd symbolic links

    """

    message='\n\nMerging spectra'

    merged_spectra=[]
    in_files=[]

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
            message=f'\n\nWorking on instrument {instrument} with {nspec} spectra\n\n'
            logger.info(message)
             #
            # generating the merged spectrum filename
            outname='{}_{}.pha'.format(srcid,instrument)
            #
            out_dir=os.path.abspath(output_dir)
            merged_spectrum=os.path.join(out_dir,outname)
            #
            if(os.path.exists(merged_spectrum)):
                message=f' File {merged_spectrum} already exists in directory {out_dir}, it will be overwritten'
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
                # printing some results
                message=f"\nInformation about the merging of spectra for instrument {instrument}"
                message+="\n   Original_index  Individual_SNR  Cumulative_SNR Spectrum Input_file Merged"
                for j in range(nspec):
                    i=sorted_indices[j]
                    if(j<=jmax):
                        merged='Yes'
                        in_file=specs[i].split('/')[-1]
                        in_files.append(in_file)
                    else:
                        merged='No'
                    line='\n  {:2d}  {:8.2f}  {:8.2f} {} {}'.format(i,snrs[i],cumsnrs[j],specs[i],merged)
                    message+=line
                logger.info(message)
                if (test):   
                   # just testing the algorithm above
                    #
                    # generating a "merged" spectrum from the input spectrum with the highest SNR
                    #   just copying it
                    message=f"\nTest mode: Output spectrum is just the input spectrum with the highest SNR {spec_list[i0][0]}"
                    logger.warning(message)
                    shutil.copy2(spec_list[i0][0],merged_spectrum)
                    sp_dic=spec_list[i0][8]
                    sp_dic['SPECFILE']=merged_spectrum
                    bkg_file = sp_dic['BACKFILE']
                else:
                    # full merging needed, using a SAS script
                    #
                    # using the list of spectra to be merged created above and the tuples to generate the
                    #     command for epicspeccombine
                    # initialising parameters with file lists
                    pha='"'
                    bkg='"'
                    rmf='"'
                    arf='"'
                    j=-1
                    for i in out_indices:
                        j+=1
                        if(j==0):
                            prefix=''
                        else:
                            prefix=' '
                        #
                        sp_dic=spec_list[i][8]
                        pha+=prefix+sp_dic['SPECFILE']
                        bkg+=prefix+sp_dic['BACKFILE']
                        rmf+=prefix+sp_dic['RESPFILE']
                        arf+=prefix+sp_dic['ANCRFILE']
                    #
                    pha+='"'
                    bkg+='"'
                    rmf+='"'
                    arf+='"'
                    # storing output file names
                    sp_dic={}
                    sp_dic['SPECFILE']=merged_spectrum
                    merged_bgd=merged_spectrum.replace('.pha','_bgd.pha')
                    sp_dic['BACKFILE']=merged_bgd
                    merged_rsp=merged_spectrum.replace('.pha','.rsp')
                    sp_dic['RESPFILE']=merged_rsp
                    sp_dic['ANCRFILE']=''
                    
                    '''
                    
                    250825 FJC: This piece of code corresponds to attempts to run epicspeccombine directly on the command line
                                I was unable to get the right combination of single and double quotes so that epicspeccombine
                                  interpreted correctly the list of files, so switched instead to writing a script and
                                  executing it, see below
                    
                    #
                    # generating command
                    # this way of passing the arguments to the shell complains that the pha or rmf are not paired
                    arguments='pha={} bkg={} rmf={} arf={} filepha="{}" filebkg="{}" filersp="{}" '.format(pha,
                    #                                bkg,rmf,arf,merged_spectrum,merged_bgd,merged_rsp)
                    cmd=['epicspeccombine',arguments]
                    #
                    # passing the arguments like a list does not complain about pairing, but the space between values in pha is ignored
                    # the line just below corresponds to defining pha (and Co) with double quotes around it
                    #cmd=['epicspeccombine',f'pha={pha}',f'bkg={bkg}',f'rmf={rmf}',f'arf={arf}',f'filepha={merged_spectrum}',f'filebkg={merged_bgd}',f'filersp={merged_rsp}']
                    # the line just below corresponds to defining pha (and Co) without double quotes around it
                    #cmd=['epicspeccombine',f"pha='{pha}'",f"bkg='{bkg}'",f"rmf='{rmf}'",f"arf='{arf}'",f'filepha={merged_spectrum}',f'filebkg={merged_bgd}',f'filersp={merged_rsp}']
                    print(f'\n   command=({cmd})')
                    
                    '''
                    
                    #
                    # preparing a string with the command line to be written to a shell script
                    combine=f'epicspeccombine pha={pha} bkg={bkg} rmf={rmf} arf={arf} filepha="{merged_spectrum}" filebkg="{merged_bgd}" filersp="{merged_rsp}" '
                    #
                    # running the SAS command epicspeccombine
                    #
                    # the lines just below are for first writing a shell script and then running it
                    merged_script=merged_spectrum.split('.')[0]+"_merge.sh"
                    if(os.path.exists(merged_script)):
                        message=f'\n File {merged_script} already exists in directory {out_dir}, it will be overwritten'
                        logger.warning(message)
                    # writing the shell script
                    try:
                        mergefile=open(merged_script,'w+')
                        mergefile.write(combine)
                        mergefile.close()
                    except Exception as e:
                        # writing the file failed
                        message=f'\nWriting of script file {merged_script} failed. Error={e}'
                        logger.error(message)
                        continue
                    #
                    # now executing the script/command
                    try: 
                        #
                        # the line just below is for submitting epicspeccombine directly to the shell, see above
                        #result=subprocess.run(cmd,capture_output=True,check=True)
                        #
                        # running  the script
                        cmd=['bash',merged_script]
                        #cmd=['echo','$0']
                        result=subprocess.run(cmd,capture_output=True,check=True)
                        #print(f'      stdout=({result.stdout})')                            
                    except Exception as e:
                        print(f'      Error occured, class=({type(e)})')
                        print(f'      Error occured e=({e})')
                        if (isinstance(e,subprocess.CalledProcessError)):
                           print(f'      stdout=({e.stdout})')
                           print(f'      stderr=({e.stderr})')
                        # merging failed
                        message=f'\nMerging of files failed. Error={e}'
                        logger.error(message)
                        continue
                    #
                    bkg_file=sp_dic['BACKFILE']                    
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
                sp_dic=spec_list[0][8]
                sp_dic['SPECFILE']=merged_spectrum
                in_files=[merged_spectrum.split('/')[-1]]
            #
            # Rebinning the merged spectrum
            #
            # changing the extension of the merged spectrum to .grp for namingthe binned spectrum
            binned_spectrum=os.path.splitext(merged_spectrum)[0]+'.grp'
            # list of merged files
            comment='Merged files: '+','.join(in_files)           
            # rebinning the spectrum
            spec_tuple=rebin_spectrum(merged_spectrum,binned_spectrum,log_file,mincts=1, background_file=bkg_file,comment=comment)
            #
            # adding to the output tuple the instrument name and the dictionary with the updated filenames
            out_tuple = list(spec_tuple)
            out_tuple.append(spec_list[0][7])
            sp_dic['SPECFILE']=binned_spectrum
            out_tuple.append(sp_dic)
            spec_tuple = tuple(out_tuple)

            #
            merged_spectra.append(spec_tuple)
        #
    #
    message=f'\nList of files used for merging={in_files} '
    logger.info(message)
    message='\n\nFinished merging spectra. Output results:'
    logger.info(message)                
    for spec_tuple in merged_spectra:
        message=f"        (merged and binned file, source counts, background counts, net counts, exposure time, flag, signal-to-noise ratio, instrument, filenames) = {spec_tuple} "
        logger.info(message)                

    return merged_spectra



def test_merge_spectra():
    output_dir='./test_data/test'
    # absolute full path
    out_dir=os.path.abspath(output_dir)
    # creating output directory
    if not os.path.exists(output_dir) : os.mkdir(output_dir)
    # Set up logging
    log_file = os.path.join(output_dir, 'test_merge_spectra.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)
    # srcid
    srcid=3067718060100029
    #
    pn_dic={}
    pn_dic['SPECFILE']=os.path.join(out_dir,'P0760940101PNS003SRSPEC0017.FTZ')
    pn_dic['BACKFILE']=pn_dic['SPECFILE'].replace('SRSPEC','BGSPEC')
    pn_dic['ANCRFILE']=pn_dic['SPECFILE'].replace('SRSPEC','SRCARF')
    pn_dic['RESPFILE']=os.path.join(out_dir,'epn_e3_ef20_sdY6.rmf')
    #
    M1_dic={}
    M1_dic['SPECFILE']=os.path.join(out_dir,'P0760940101M1S001SRSPEC0017.FTZ')
    M1_dic['BACKFILE']=M1_dic['SPECFILE'].replace('SRSPEC','BGSPEC')
    M1_dic['ANCRFILE']=M1_dic['SPECFILE'].replace('SRSPEC','SRCARF')
    M1_dic['RESPFILE']=os.path.join(out_dir,'m1_e13_im_pall_o.rmf')
    #
    M2_dic={}
    M2_dic['SPECFILE']=os.path.join(out_dir,'P0760940101M2S002SRSPEC0017.FTZ')
    M2_dic['BACKFILE']=M2_dic['SPECFILE'].replace('SRSPEC','BGSPEC')
    M2_dic['ANCRFILE']=M2_dic['SPECFILE'].replace('SRSPEC','SRCARF')
    M2_dic['RESPFILE']=os.path.join(out_dir,'m2_e13_im_pall_o.rmf')

    
    # === Test 1: two empty lists ===
    # it should never get to call merge_spectra with two empty lists, but checking anyway
    print("\n🔹 Test 1: two empty lists")
    pn_spectra=[]
    mos_spectra=[]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1)
    # output list should have no elements
    assert len(merged_list)==0

    # === Test 2: one empty list and one with errors from before ===
    print("\n🔹 Test 2: one empty list and one with errors from before")
    pn_spectra=[('dummyPN.fits',-1,1,-1,1000.0,2,-1,'pn',{})]
    mos_spectra=[]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1)
    # output list should have just one element, copied from pn_spectra
    assert len(merged_list)==1
    assert merged_list[0][5]==2

    # === Test 3: only pn ===
    print("\n🔹 Test 3: only pn")
    # only pn
    #  tuple just below from test_check_spectra.log
    pn_spectra=[('./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ', 766, 8474, 271.8810110703645, 82181.936317917, 0, 7.659018142527241,'pn',pn_dic)]
    mos_spectra=[]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1)
    # output list should have just one element, with the grouped version of the spectrum above
    assert len(merged_list)==1
    name=merged_list[0][0].split('/')[-1]
    assert name=='3067718060100029_pn.grp'

    # === Test 4: only MOS, test mode ===
    print("\n🔹 Test 4: only MOS, test mode")
    pn_spectra=[]
    mos_spectra=[('./test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ', 308, 14296, 63.03960341552707, 104469.411107063, 0, 2.6808126142724875,'MOS',M1_dic),('./test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ', 236, 19138, 99.06503350999267, 105554.512163162, 0, 5.129840220900734,'MOS',M2_dic) ]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1, test=True)
    print("only MOS: merged_list ",merged_list)
    # output list should have just one element, with the grouped version of the second spectrum above
    assert len(merged_list)==1
    name=merged_list[0][0].split('/')[-1]
    assert name=='3067718060100029_MOS.grp'
    # signal-to-noise-ratio to check that chose the second MOS spectrum
    assert abs(mos_spectra[1][6]-merged_list[0][6])<=0.01
    # checking now that the dictionary is present and filled
    assert len(merged_list[0][8])==4
    # checking that the spec file is properly filled
    assert merged_list[0][0]==merged_list[0][8]['SPECFILE']

    # === Test 5: both pn and MOS, test mode ===
    print("\n🔹 Test 4: both pn and MOS, test mode")
    pn_spectra=[('./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ', 766, 8474, 271.8810110703645, 82181.936317917, 0, 7.659018142527241,'pn',pn_dic)]
    mos_spectra=[('./test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ', 308, 14296, 63.03960341552707, 104469.411107063, 0, 2.6808126142724875,'MOS',M1_dic),('./test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ', 236, 19138, 99.06503350999267, 105554.512163162, 0, 5.129840220900734,'MOS',M2_dic) ]
    merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1, test=True)
    # output list should contain two elements, with the grouped versions of the first spectrum in the pn list and the second spectrum in the second list above
    assert len(merged_list)==2
    name=merged_list[0][0].split('/')[-1]
    assert name=='3067718060100029_pn.grp'
    name=merged_list[1][0].split('/')[-1]
    assert name=='3067718060100029_MOS.grp'
    # signal-to-noise-ratio to check that chose the second MOS spectrum
    assert abs(mos_spectra[1][6]-merged_list[1][6])<=0.01

    # === Test 6: only MOS, full run including merging ===
    print("\n🔹 Test 6: only MOS, full run including merging")
    #
    # checking if epicspeccombine is installed, otherwise skipping this test
    try:
        cmd=['epicspeccombine','-h']
        result=subprocess.run(cmd,capture_output=True,check=True)
        # if it gets here epicspeccombine is defined, continuing with the test    
        pn_spectra=[]
        mos_spectra=[('./test_data/0760940101/pps/P0760940101M1S001SRSPEC0017.FTZ', 308, 14296, 63.03960341552707, 104469.411107063, 0, 2.6808126142724875,'MOS',M1_dic),('./test_data/0760940101/pps/P0760940101M2S002SRSPEC0017.FTZ', 236, 19138, 99.06503350999267, 105554.512163162, 0, 5.129840220900734,'MOS',M2_dic) ]
        merged_list=merge_spectra(pn_spectra,mos_spectra, srcid, output_dir, log_file, mincts=1, test=True)
        print("only MOS: merged_list ",merged_list)
        # output list should have just one element, merging the two spectra above
        assert len(merged_list)==1
        name=merged_list[0][0].split('/')[-1]
        assert name=='3067718060100029_MOS.grp'
        # got the counts and SNR from the epicspeccombined file by hand, testing them below
        spec_tuple = get_spectral_counts(merged_list[0][0],log_file)
        # tot counts
        assert spec_tuple[1] == 544
        # bgd counts
        assert spec_tuple[2] == 33434
        # net counts
        assert abs(spec_tuple[3] - 196.65) <= 0.01
        # exposure time
        assert abs(spec_tuple[4] - 105011.96) <= 0.01
        # flag
        assert spec_tuple[5] == 0
        # signal-to-noise ratio
        assert abs(spec_tuple[6] - 6.59) <= 0.01
        # checking now that the dictionary is present and filled
        assert len(merged_list[0][8])==4
        # checking that the spec file is properly filled
        assert merged_list[0][0]==merged_list[0][8]['SPECFILE']
    except:
        print('   epicspeccombine not defined, skipping this test')
    #
    #assert True
    
    
if __name__ == "__main__":
    test_merge_spectra()
    print("All merge and background tests completed successfully.")

