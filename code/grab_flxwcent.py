# -*- coding: utf-8 -*-
"""
The routine gets the flux weighted centroids for sectors
Supports multi sector
Will also be used to gather PDC goodness statistics when
I can talk with Jeff Smith.
The Time re-sampling of these is on an hour time scale

AUTHOR: Christopher J. Burke
"""

import numpy as np
import pickle
from gather_tce_fromdvxml import tce_seed
import os
import math
import h5py
from statsmodels import robust
import matplotlib.pyplot as plt
from astropy.io import fits
import argparse

def make_data_dirs(prefix, sector, epic):
    secDir = 'S{0:02d}'.format(sector)
    localDir = os.path.join(prefix,secDir)
    if not os.path.exists(localDir):
        os.mkdir(localDir)
    epcDir = '{0:04d}'.format(int(math.floor(epic/1000.0)))
    localDir = os.path.join(prefix,secDir,epcDir)
    if not os.path.exists(localDir):
        os.mkdir(localDir)
    return localDir

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list

# parse a page/sector list
# from https://stackoverflow.com/a/4248689
def parse_range(astr):
    result=set()
    for part in astr.split(','):
        x=part.split('-')
        result.update(range(int(x[0]),int(x[-1])+1))
    return sorted(result)


if __name__ == '__main__':
    # Parse the command line arguments for multisector sector specifications
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",nargs='?',\
                        help="String that specificies the sectors involved in the multisector run. It works like a printer page specification. For example '14-26,40-50' would include sectors 14-26 and sectors 40-50")
    

    args = parser.parse_args()
    secstr = args.s
    if not secstr is None:
        sectorswant = parse_range(secstr)
        maxsec = np.max(sectorswant)
        print('Requesting Sectors ',sectorswant)
    
    dirOutputs = '/pdo/users/cjburke/spocvet/sector54/'
    RESAMP = 31  ###  USE AN ODD NUMBER ###
    SECTOR = 54# =-1 if multi-sector

    if (not 'maxsec' in locals()) and SECTOR == -1:
        print('Code expects multi-sector (SECTOR = 54) but no sector argument was given. EXITING!')
        sys.exit()

    if SECTOR > 0:
        print('Single Sector')
        # Single sector block of file prefixes fill with fake values
        fileInputPrefixList = []
        for i in np.arange(1,SECTOR):
            fileInputPrefixList.append('/foo{0:d}'.format(i))
        fileInputPrefixList.append('/pdo/spoc-data/sector-054/light-curve/tess2022190063128-s0054-')
        fileInputSuffixList = []
        for i in np.arange(1,SECTOR):
            fileInputSuffixList.append('/foo{0:d}'.format(i))
        fileInputSuffixList.append('-0227-s_lc.fits.gz')
    else:
        print('Multi-Sector')
        # multisector read in table that has all the paths to data files
        dataBlock = np.genfromtxt('tec_datafile_paths.dat',dtype=['i8','U80','U80','i8'])
        sectorNumbers = dataBlock['f0']
        pathPrefixes = dataBlock['f1']
        filePrefixes = dataBlock['f2']
        fileSuffixNumbers = dataBlock['f3']
        fileInputPrefixList = []
        fileInputSuffixList = []
        for i in np.arange(1,maxsec+1):
            if i in sectorswant:
                ia = np.where(sectorNumbers == i)[0]
                if len(ia) > 0:
                    curStr = os.path.join(pathPrefixes[ia][0],'light-curve','{0}-s{1:04d}-'.format(filePrefixes[ia][0],i))
                    fileInputPrefixList.append(curStr)
                    curStr = '-{0:04d}-s_lc.fits.gz'.format(fileSuffixNumbers[ia][0])
                    fileInputSuffixList.append(curStr)
                else:
                    fileInputPrefixList.append('/foo{0:d}'.format(i))
                    fileInputSuffixList.append('/foo{0:d}'.format(i))
            else:
                fileInputPrefixList.append('/foo{0:d}'.format(i))
                fileInputSuffixList.append('/foo{0:d}'.format(i))

    nSector = len(fileInputPrefixList)    

    #fileOut = 'spoc_pdcstats_sector-54_20220907.txt'
    #fom = open(fileOut, 'w')
    vetFile = 'spoc_fluxtriage_sector-54_20220907.txt'
    #vetFile = 'junk.txt'
    tceSeedInFile = 'sector-54_20220907_tce.h5'

    # Load the tce data h5
    tceSeedInFile = 'sector-54_20220907_tce.h5'
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=int)
    allrp = np.array([x.at_rp for x in all_tces])
    allrstar = np.array([x.rstar for x in all_tces])
    alllogg = np.array([x.logg for x in all_tces])
    allper = np.array([x.at_period for x in all_tces])
    alltmags = np.array([x.tmag for x in all_tces])
    allmes = np.array([x.mes for x in all_tces])
    allsnr = np.array([x.at_snr for x in all_tces])
    alldur = np.array([x.at_dur for x in all_tces])
    allsolarflux = np.array([x.at_effflux for x in all_tces])
    allatdep = np.array([x.at_depth for x in all_tces])
    allatepoch = np.array([x.at_epochbtjd for x in all_tces])
    alltrpvalid = np.array([x.trp_valid for x in all_tces])
    allatrpdrstar = np.array([x.at_rpDrstar for x in all_tces])
    allatrpdrstare = np.array([x.at_rpDrstar_e for x in all_tces])
    allatadrstar = np.array([x.at_aDrstar for x in all_tces])

    # Load the  flux vetting
    dataBlock = np.genfromtxt(vetFile, dtype=[int,int,int,'S1'])
    fvtic = dataBlock['f0']
    fvpn = dataBlock['f1']
    fvvet = dataBlock['f2']
    
    allvet = np.zeros_like(allpn)
    for i in range(len(allvet)):
        idx = np.where((alltic[i] == fvtic) & (allpn[i] == fvpn))[0]
        if len(idx) > 0:
            allvet[i] = fvvet[idx]
    # only keep tces with both valid dv and trapezoid fits
    # and flux vetted pass
    idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
    alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
            allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
            allatrpdrstar, allatrpdrstare, allatadrstar = idx_filter(idx, \
            alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
            allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
            allatrpdrstar, allatrpdrstare, allatadrstar)
            
    # Get flux weighted centroids  and PDC stats over flux triage passing TCEs
    for i, curTic in enumerate(alltic):
        print('{:d} of {:d}'.format(i, len(alltic)))
        curPn = allpn[i]
        hasSector = np.zeros((nSector,), dtype=int)
        pdcTot = np.zeros((nSector,), dtype=float)
        pdcNoise = np.zeros((nSector,), dtype=float)
        pdcCor = np.zeros((nSector,), dtype=float)
        
        # Find files for each sector
        for k in range(nSector):
            fileInput = '{0}{1:016d}{2}'.format(fileInputPrefixList[k], curTic, fileInputSuffixList[k])
            fileOutput = os.path.join(make_data_dirs(dirOutputs, SECTOR, curTic), 'tess_flxwcent_{0:016d}_{1:02d}_{2:02d}.h5d'.format(curTic,curPn, k+1))
            if os.path.isfile(fileInput):
                hdulist = fits.open(fileInput)
                
                # Get some PDC statistics
                hasSector[k] = 1
                try:
                    pdcTot[k] = hdulist[1].header['PDC_TOT']
                except:
                    print('Cannot find PDC_TOT {0:d} sec: {1:d}'.format(curTic, k))
                pdcNoise[k] = hdulist[1].header['PDC_NOI']
                pdcCor[k] = hdulist[1].header['PDC_COR']
                pdcStats = np.array([pdcTot[k], pdcNoise[k], pdcCor[k]], dtype=float)
                flxw_centr1 = hdulist[1].data['MOM_CENTR1']
                flxw_centr2 = hdulist[1].data['MOM_CENTR2']
                time = hdulist[1].data['TIME']
                dqflgs = hdulist[1].data['QUALITY']
                nImage = len(time)
                newNImage = int(np.floor(nImage / RESAMP))
                oldNImage = newNImage*RESAMP
                # trim off the excess images not integral into resamp
                idx = np.arange(0,oldNImage)
                flxw_centr1, flxw_centr2, time, dqflgs = idx_filter(idx, \
                    flxw_centr1, flxw_centr2, time, dqflgs)
        
                # Do downsampling of data stream
                time = np.mean(np.reshape(time, (newNImage, RESAMP)), axis=1)
                flxw_centr1 = np.median(np.reshape(flxw_centr1, (newNImage, RESAMP)), axis=1)
                flxw_centr2 = np.median(np.reshape(flxw_centr2, (newNImage, RESAMP)), axis=1)
                dqflgs = np.bitwise_or.reduce(np.reshape(dqflgs, (newNImage, RESAMP)), axis=1)
    
                # Identify data that is missing or NaN
                idx = np.where((np.isfinite(time)) & (np.isfinite(flxw_centr1)) & (np.isfinite(flxw_centr2)))[0]
                valid_data_flag = np.zeros((newNImage,), dtype=np.bool_)
                valid_data_flag[idx] = True
                idx = np.where(dqflgs == 0)[0]
                valid_data_flag[idx] = True 
                # Now save data as h5py 
                f = h5py.File(fileOutput, 'w')
                tmp = f.create_dataset('time', data=time, compression='gzip') 
                tmp = f.create_dataset('flxw_centr1', data=flxw_centr1, compression='gzip')
                tmp = f.create_dataset('flxw_centr2', data=flxw_centr2, compression='gzip')
                tmp = f.create_dataset('dqflgs', data=dqflgs, compression='gzip')
                tmp = f.create_dataset('valid_data_flag', data=valid_data_flag, compression='gzip')
                tmp = f.create_dataset('pdc_stats', data=pdcStats)
                
                print(curTic, 'alpha')                


#        fo = open(fileOutput, 'w')
#        # Write out the best fit parameters as a header
#        #for j in range(len(ioblk.physvals)): 
#        #    fo.write('# {:s} {:f}\n'.format(ioblk.physval_names[j], ioblk.bestphysvals[j]))
#
#        for j in range(len(ioblk.normts)):
#            strout = '{:f} {:f} {:f}\n'.format(ioblk.normts[j], ioblk.normlc[j]-1.0, ioblk.modellc[j]-1.0)
#            fo.write(strout)
#        fo.close()
#        
#
#    fom.close()
