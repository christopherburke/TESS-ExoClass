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



if __name__ == '__main__':
    #  Directory list for Sector light curve files
    # This block is for the multi-sector case
# In the case of a single sector One needs dummy entries for
#  every sector
#    fileInputPrefixList = ['/pdo/spoc-data/sector-001-20210219/light-curve/tess2018206045859-s0001-', \
#                          '/pdo/spoc-data/sector-002-20210219/light-curve/tess2018234235059-s0002-', \
#                          '/pdo/spoc-data/sector-003-20210219/light-curve/tess2018263035959-s0003-', \
#                          '/pdo/spoc-data/sector-004-20210219/light-curve/tess2018292075959-s0004-', \
#                          '/pdo/spoc-data/sector-005-20210219/light-curve/tess2018319095959-s0005-', \
#                          '/pdo/spoc-data/sector-006-20210219/light-curve/tess2018349182459-s0006-', \
#                          '/pdo/spoc-data/sector-007/light-curve/tess2019006130736-s0007-', \
#                          '/pdo/spoc-data/sector-008/light-curve/tess2019032160000-s0008-', \
#                          '/pdo/spoc-data/sector-009/light-curve/tess2019058134432-s0009-', \
#                          '/pdo/spoc-data/sector-010/light-curve/tess2019085135100-s0010-', \
#                          '/pdo/spoc-data/sector-011/light-curve/tess2019112060037-s0011-', \
#                          '/pdo/spoc-data/sector-012/light-curve/tess2019140104343-s0012-', \
#                          '/pdo/spoc-data/sector-013/light-curve/tess2019169103026-s0013-',\
#                          '/foo14','/foo15','/foo16','/foo17','/foo18',\
#                          '/foo19','/foo20','/foo21','/foo22','/foo23',\
#                          '/foo24','/foo25','/foo26',\
#                          '/pdo/spoc-data/sector-027/light-curve/tess2020186164531-s0027-',\
#                          '/pdo/spoc-data/sector-028/light-curve/tess2020212050318-s0028-',\
#                          '/pdo/spoc-data/sector-029/light-curve/tess2020238165205-s0029-',\
#                          '/pdo/spoc-data/sector-030/light-curve/tess2020266004630-s0030-',\
#                          '/pdo/spoc-data/sector-031/light-curve/tess2020294194027-s0031-',\
#                          '/pdo/spoc-data/sector-032/light-curve/tess2020324010417-s0032-',\
#                          '/pdo/spoc-data/sector-033/light-curve/tess2020351194500-s0033-',\
#                          '/pdo/spoc-data/sector-034/light-curve/tess2021014023720-s0034-',\
#                          '/pdo/spoc-data/sector-035/light-curve/tess2021039152502-s0035-',\
#                          '/pdo/spoc-data/sector-036/light-curve/tess2021065132309-s0036-',\
#                          '/pdo/spoc-data/sector-037/light-curve/tess2021091135823-s0037-',\
#                          '/pdo/spoc-data/sector-038/light-curve/tess2021118034608-s0038-',\
#                          '/pdo/spoc-data/sector-039/light-curve/tess2021204101404-s0041-']

#    fileInputSuffixList = ['-0120-s_lc.fits.gz', \
#                           '-0121-s_lc.fits.gz', \
#                           '-0123-s_lc.fits.gz', \
#                           '-0124-s_lc.fits.gz', \
#                           '-0125-s_lc.fits.gz', \
#                           '-0126-s_lc.fits.gz', \
#                           '-0131-s_lc.fits.gz', \
#                           '-0136-s_lc.fits.gz', \
#                           '-0139-s_lc.fits.gz', \
#                           '-0140-s_lc.fits.gz', \
#                           '-0143-s_lc.fits.gz', \
#                           '-0144-s_lc.fits.gz', \
#                           '-0146-s_lc.fits.gz',\
#                           '/foo14','/foo15','/foo16','/foo17','/foo18',\
#                           '/foo19','/foo20','/foo21','/foo22','/foo23',\
#                           '/foo24','/foo25','/foo26',\
#                           '-0189-s_lc.fits.gz',\
#                           '-0190-s_lc.fits.gz',\
#                           '-0193-s_lc.fits.gz',\
#                           '-0195-s_lc.fits.gz',\
#                           '-0198-s_lc.fits.gz',\
#                           '-0200-s_lc.fits.gz',\
#                           '-0203-s_lc.fits.gz',\
#                           '-0204-s_lc.fits.gz',\
#                           '-0205-s_lc.fits.gz',\
#                           '-0207-s_lc.fits.gz',\
#                           '-0208-s_lc.fits.gz',\
#                           '-0209-s_lc.fits.gz',\
#                           '-0212-s_lc.fits.gz']

# Single sector block of file prefixes fill with fake values
    fileInputPrefixList = ['/foo1','/foo2','/foo3','/foo4','/foo5',\
                           '/foo6','/foo7','/foo8','/foo9','/foo10',\
                           '/foo11','/foo12','/foo13','/foo14','/foo15',
                           '/foo16','/foo17','/foo18','/foo19','/foo20',\
                           '/foo21', '/foo22','/foo23','/foo24','/foo25',\
                           '/foo26','/foo27','/foo28','/foo29','/foo30',\
                           '/foo31','/foo32','/foo33','/foo34','/foo35',\
                           '/foo36','/foo37','/foo38','/foo39','/foo40',\
                           '/pdo/spoc-data/sector-041/light-curve/tess2021204101404-s0041-']
    fileInputSuffixList = ['/foo1','/foo2','/foo3','/foo4','/foo5',\
                           '/foo6','/foo7','/foo8','/foo9','/foo10',\
                           '/foo11','/foo12','/foo13','/foo14','/foo15',
                           '/foo16','/foo17','/foo18','/foo19','/foo20',\
                           '/foo21','/foo22','/foo23','/foo24','/foo25',\
                           '/foo26','/foo27','/foo28','/foo29','/foo30',\
                           '/foo31','/foo32','/foo33','/foo34','/foo35',\
                           '/foo36','/foo37','/foo38','/foo39','/foo40',\
                           '-0212-s_lc.fits.gz']
    nSector = len(fileInputPrefixList)    
    dirOutputs = '/pdo/users/cjburke/spocvet/sector41/'
    RESAMP = 31  ###  USE AN ODD NUMBER ###
    SECTOR = 41# =-1 if multi-sector

    #fileOut = 'spoc_pdcstats_sector41_20210917.txt'
    #fom = open(fileOut, 'w')
    vetFile = 'spoc_fluxtriage_sector41_20210917.txt'
    #vetFile = 'junk.txt'
    tceSeedInFile = 'sector41_20210917_tce.h5'

    # Load the tce data h5
    tceSeedInFile = 'sector41_20210917_tce.h5'
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=np.int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=np.int)
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
        hasSector = np.zeros((nSector,), dtype=np.int)
        pdcTot = np.zeros((nSector,), dtype=np.float)
        pdcNoise = np.zeros((nSector,), dtype=np.float)
        pdcCor = np.zeros((nSector,), dtype=np.float)
        
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
                pdcStats = np.array([pdcTot[k], pdcNoise[k], pdcCor[k]], dtype=np.float)
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