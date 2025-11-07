#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import numpy as np
from astropy.io import fits
import logging

from get_spectral_counts import get_spectral_counts


logger = logging.getLogger(__name__)


# Function to rebin a spectrum, writing to disk a copy of it and returning some basic info
def rebin_spectrum(infile,outfile,log_file,mincts=1):
    """
    Reads in an input file, rebins it to have >=mincts counts per bin, taking into account
        if there are bad channels at the beginning, and flagging as "bad" empty channels at the end.
    It writes the rebinned file as outfile to disk, and returns a tuple with the
        full path to the output file, total counts in the spectrum, total background counts, total net counts, exposure time, 
        and a flag which is 1 if the background counts are <=0, 2 if the net counts are <=0, and 0 otherwise    
        
    Parameters:
    - infile (str): input spectrum file
    - outfile (str): ouput rebinned spectrum file
    - log_file (str): The log file to write the messages.
    - mincts (int): minimum number of counts per bin (default=1)
    Returns:
    - spec_tuple (tuple): A tuple containing the full path and name for the output spectrum, total source counts, total background counts, total net counts, exposure time, a flag, and the signal-to-noise ratio

    """

    message=f'\n\n Rebinning file {infile} with a minimum of {mincts} counts per bin'
    logger.info(message)


    
    #
    # full path to output spectrum
    outfile=os.path.abspath(outfile)
    #output_dir=os.path.dirname(fullname)


    #
    # getting the values for the output tuple, containing counts, exposure, a flag and the signal-to-noise ratio
    #   once generated, tuples cannot be changed, so going through a list to change the name
    #        of the spectral file
    #
    spec_tuple=list(get_spectral_counts(infile,log_file,background_file='path'))
    spec_tuple[0]=outfile
    spec_tuple=tuple(spec_tuple)    

    #
    # checking if input file can be opened
    if(spec_tuple[5]==-2):
        # cannot open input file, skipping the rest of the function
        message=f' Cannot open FITS spectrum file {infile}'
        logger.error(message)
        return spec_tuple
    
    
    ########
    #
    # Reading the input file header and data
    #
    with fits.open(infile) as hdul:
        data1=hdul[1].data
        header0=hdul[0].header.copy()
        header1=hdul[1].header.copy()
    
    #
    # reading count information
    counts=data1['COUNTS']
    nchan=len(counts)
      
    #
    # checking if quality info provided, otherwise fill with 0 (good quality)
    try:
        quality=data1['QUALITY']
        hasQuality=True
        nbad=0
        # counting how many channels at the beginning are marked as bad
        for i in range(nchan):
            if(quality[i]>0):
                nbad+=1
            else:
                break
            #
        #
        message=f'    {nbad} initial channels are marked as bad'
        logger.info(message)
        #
        badmask=quality>0
        nbadmask=sum(badmask)
        message=f'    {nbadmask} total channels are marked as bad'
        logger.info(message)
        if(nbadmask>nbad):
            message='    Additional bad channels not at the beginning'
            logger.warning(message)
    except:
        # no quality information, setting all channels as good
        hasQuality=False
        quality=np.full(nchan,0)
        nbad=0
    #
    
    #
    # checking if grouping info provided
    try:
        hasGrouping=True
        grouping=data1['GROUPING']
    except:
        # no grouping information 
        hasGrouping=False
    #
    # reinitializing grouping with one channel per group
    grouping=np.full(nchan,1)
    
    
    #
    # Redoing the grouping, so that each group has >=1 count
    # cumulative number of counts in the current group
    groupCounts=0
    # is this the first channel of the current group?
    first=True
    # index of the first channel of the current group
    first_i=0
    for i in range(nchan):
        if (i>=nbad):
            if(i== nbad-1):first_i=i
            #print('   i {} cts {} '.format(i,counts[i]))
            # not one of the bad channels at the beginning
            groupCounts+=counts[i]
            #print('      groupCounts {} '.format(groupCounts))
            if (groupCounts>=mincts):
                # minimum number of counts reached, this is the last channel of this group
                if (first):
                    # this is a single channel with =>1 counts, no need to change grouping
                    #     already set to 1
                    pass
                else:
                    # adding this channel to the existing group
                    grouping[i]=-1                
                #
                # initialising the cumulative counts to 0 for the next channel and group
                #    since the current group is closed here
                groupCounts=0
                # initialising next channel as first of its group
                first=True
                # initialising index of first channel of the next group
                first_i=i+1
                #print('         new Group grouping_i {} '.format(grouping[i]))
            elif (groupCounts<1):
                # not enough counts in this channel to close the group
                #
                if (first):
                    # this is the first channel of this new group
                    #     grouping already set to 1
                    pass
                else:
                    # still need to add this channel to the group
                    grouping[i]=-1
                # 
                # in any case, the next channel is not going to be the first one of the group
                first=False
                #print('         existing Group grouping_i {} '.format(grouping[i]))
            #
        else:
            # initial bad channel, going for the next
            continue
        #
    #
    # if the last group does not reach the minimum number of counts, setting all channels in it to bad
    #     and setting each channel as a group
    # this behaviour mimics the result of grppha
    if (groupCounts<mincts):
        grouping[first_i:]=1
        quality[first_i:]=2
        message=f'    Setting channels {first_i} to last as bad (quality=2)'
        logger.info(message)
    #
    
    
    
    #
    # Updating the grouping information
    if (hasGrouping):
        # just updating the grouping
        data1['GROUPING']=grouping
        # creating new extension with existing columns and updated one
        hdu=fits.BinTableHDU.from_columns(data1.columns, header=header1)
    else:
        # creating new column
        colGrouping=fits.colDefs([fits.Colum(name='GROUPING',format='I',array=grouping)])
        # creating new extension with existing columns + new one
        hdu=fits.BinTableHDU.from_columns(data1.columns+colGrouping, header=header1)
    #
    
    #
    # Updating the quality
    if (hasQuality):
        # just updating the grouping
        hdu.data['QUALITY']=quality
        # creating new extension with existing columns and updated one
        hdu=fits.BinTableHDU.from_columns(hdu.data.columns, header=hdu.header)
    else:
        # creating new column
        colQuality=fits.colDefs([fits.Colum(name='QUALITY',format='I',array=quality)])
        # creating new extension with existing columns + new one
        #    at this stage hdu has already been created in both options for grouping
        hdu=fits.BinTableHDU.from_columns(hdu.columns+colQuality,header=hdu.header)
    #
    
    # renaming hdu as SPECTRUM
    hdu.name='SPECTRUM'
    
    '''
    
    Disconnecting this for the time being, since the names are long and an '&' is added to the
       end of the filename, which makes xspec crash
    
    # following the symbolic links for the ancillary files
    for keyword in ['respfile','ancrfile','backfile']:
        filename=hdu.header[keyword]
        fullname=os.path.join(output_dir,filename)
        if (os.path.islink(fullname)):
            # need to check this for all files, because merged files will have a filename, not a link 
            fullname=os.readlink(fullname)
        print('   keyword ({})  fullname ({}) '.format(keyword,fullname))
        hdu.header[keyword]=fullname
    '''
    
    #creating a new primary HDU with the original primary header
    hdu0=fits.PrimaryHDU(header=header0)
    
    # concatenating the two hdu
    hdu_new = fits.HDUList([hdu0, hdu])
    #
    # writing the output file
    hdu_new.writeto(outfile,overwrite=True)
    message=f' {nchan} channels written out to file {outfile}'
    logging.info(message)
    
    

    #
    # cleaning up to free memory
    del hdul,hdu,hdu_new,hdu0,grouping,quality,counts,data1,header0,header1
    

    #      
    return spec_tuple


def test_rebin_spectrum():
    output_dir='./test_data/test'
    if not os.path.exists(output_dir) : os.mkdir(output_dir)
    # Set up logging
    log_file = os.path.join(output_dir, 'test_rebin_spectrum.txt')
    logging.basicConfig(filename=log_file, level=logging.INFO)
    #
    # inexistent spectrum file, testing if the function captures this case
    infile='dummy_file_FJC.txt'
    spec_tuple=rebin_spectrum(infile,'',log_file)
    assert spec_tuple[5]==-2
    #
    # existent spectrum file
    mincts=1
    infile='./test_data/0760940101/pps/P0760940101PNS003SRSPEC0017.FTZ'
    outfile=os.path.join(output_dir,'test_rebin_spectrum.grp1')
    spec_tuple=rebin_spectrum(infile,outfile,log_file,mincts=mincts)
    # no error
    assert spec_tuple[5]==0
    # ouput file exists
    assert os.path.exists(outfile)
    # properties of output file
    with fits.open(outfile) as sp_hdul:
        counts = sp_hdul[1].data['COUNTS']
        quality = sp_hdul[1].data['QUALITY']
        grouping = sp_hdul[1].data['GROUPING']
        #
        good = quality==0
        tot_good_cts=sum(counts[good])
        message=f'\n\ngrouping[good] {grouping[good]}'
        logger.info(message)
        grouping_good=grouping[good]
        tot_good_bins=sum( grouping_good==1 )
        #
        # comparing to values from same spectrum binned with grppha
        # total number of counts in good bins
        assert tot_good_cts==730
        # total number of good bins
        assert tot_good_bins==546
        #
        good_cts=counts[good]
        good_grouping=grouping[good]
        ngood=len(good_cts)
        tot_cts=good_cts[0]
        for i in range(1,ngood):
            if(good_grouping[i]==1):
                # new bin
                # previous bin has to have total counts >=mincts
                assert tot_cts>=mincts
                # initialising total counts in new bin
                tot_cts=good_cts[i]
            else:
                # still previous bin, increasing counts
                tot_cts+=good_cts[i]
            #
        #
        # checking last bin
        assert tot_cts>=mincts
    # correctly calculated signal to noise ratio
    # net number of counts after subtracting background, in all channels
    net=271.88
    # total number of counts in all channels
    tot=766
    # signal-to-noise ratio
    snr=net/np.sqrt(2*tot-net)
    assert abs(spec_tuple[6]-snr)<=0.01
            
                    
        
    #



