#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 19:05:23 2018

@author: cjburke
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import h5py
from statsmodels import robust
import glob
import os
import math
from gather_tce_fromdvxml import tce_seed
import cjb_utils as cjb
import argparse
import sys

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

def tpf_resamp(file, fileOut, RESAMP, lcFile):
    """ Resample TESS target pixel file and save as h5d format
        resamp - Resample factor just make it odd okay"""    
    hdulist = fits.open(file)
    arr = hdulist[1].data[0]['FLUX']
    nImage = len(hdulist[1].data[:]['CADENCENO'])
    shp = arr.shape
    nx = shp[0]
    ny = shp[1]
    saturate_pixel = np.zeros((nx, ny), dtype=int)
    median_image = np.zeros((nx, ny))

    # Get header information that we should keep
    keepprihdr = ['TICID','SECTOR','CAMERA','CCD','PXTABLE','RA_OBJ', \
                  'DEC_OBJ','PMRA','PMDEC','PMTOTAL','TESSMAG','TEFF', \
                  'LOGG','RADIUS']
    formatprihdr = [np.uint32, int, int, int,int, \
                    float,float,float,float,float, \
                    float,float,float,float,float]

    keep1hdr = ['1CRV4P','2CRV4P','1CRPX4','2CRPX4']
    format1hdr = [int,int,float,float]
    
    cadenceNo = hdulist[1].data[:]['CADENCENO']
    timetbjd = hdulist[1].data[:]['TIME']
    flux_array = hdulist[1].data[:]['FLUX']
    flux_bkg_array = hdulist[1].data[:]['FLUX_BKG']
    dq_flag = hdulist[1].data[:]['QUALITY']
    
    f = h5py.File(lcFile,'r')
    cadNo = np.array(f['cadenceNo'])
    cadNoBeg = np.array(f['cadenceNoBeg'])
    cadNoEnd = np.array(f['cadenceNoEnd'])
    ia, ib = cjb.intersect(cadNoBeg, cadenceNo)
    # In rare instances the light curve data doesnt exist for this sector and
    #  ib will be empty causeing error if so
    #  just return and do nothing for this sector
    try:
        frstIdx = ib[0]
    except:
        return
    ia, ib = cjb.intersect(cadNoEnd, cadenceNo)    
    endIdx = ib[-1]    

    #Make a fix for Sector 3 where not all cadences were used
    # in the backend DV
    #idx = np.where((cadenceNo>=114115) & (cadenceNo<=128706))[0]
    #nImage = len(idx)
    #cadenceNo, timetbjd, dq_flag = idx_filter(idx, cadenceNo, timetbjd, dq_flag)
    #flux_array = flux_array[idx, :, :]
    #flux_bkg_array = flux_bkg_array[idx, :, :]

    # trim off the excess images not integral into resamp
    cadenceNo = cadenceNo[frstIdx:endIdx+1]
    timetbjd = timetbjd[frstIdx:endIdx+1]
    flux_array = flux_array[frstIdx:endIdx+1, :, :]
    flux_bkg_array = flux_bkg_array[frstIdx:endIdx+1, :,:]
    dq_flag = dq_flag[frstIdx:endIdx+1]
    newNImage = len(cadenceNo) // RESAMP    
    # Do downsampling of data stream
    cadenceNo = np.mean(np.reshape(cadenceNo, (newNImage, RESAMP)), axis=1, dtype=int)
    timetbjd = np.mean(np.reshape(timetbjd, (newNImage, RESAMP)), axis=1)
    flux_array = np.sum(np.reshape(flux_array, (newNImage, RESAMP, nx, ny)), axis=1)
    flux_bkg_array = np.sum(np.reshape(flux_bkg_array, (newNImage, RESAMP, nx, ny)), axis=1)
    dq_flag = np.sum(np.reshape(dq_flag, (newNImage, RESAMP)), axis=1, dtype=int)

    # Identify data that is missing or NaN
    idx = np.where((np.isfinite(timetbjd)) & (np.isfinite(np.squeeze(flux_array[:,0,0]))) & (np.isfinite(np.squeeze(flux_bkg_array[:,0,0]))))[0]
    valid_data_flag = np.zeros((newNImage,), dtype=np.bool_)
    valid_data_flag[idx] = True
    


    # Identify saturated pixels
    for i in range(nx):
        for j in range(ny):
            curflux = flux_array[:,i,j]
            diff_flux = np.diff(curflux[valid_data_flag])
            robmad = robust.mad(diff_flux)
            medval = np.median(curflux[valid_data_flag])
            median_image[i,j] = medval
            if medval > 1000.0 and np.log10(robmad/medval) < -3.5:
                saturate_pixel[i,j] = 1
#                print("Saturated Pixel detected x: {0:d} y: {1:d}".format(i, j))

 
    # Now save data as h5py
    epic = hdulist[0].header['TICID']
    sec = hdulist[0].header['SECTOR']
#    fileoutput = os.path.join(make_data_dirs(dirOut,sec,epic), 'tess_tpf_{0:016d}.h5d'.format(epic))
    f = h5py.File(fileOut, 'w')
    tmp = f.create_dataset('cadenceNo', data=cadenceNo, compression='gzip') 
    tmp = f.create_dataset('timetbjd', data=timetbjd, compression='gzip')
    tmp = f.create_dataset('flux_array', data=flux_array, compression='gzip')
    tmp = f.create_dataset('flux_bkg_array', data=flux_bkg_array, compression='gzip')
    tmp = f.create_dataset('dq_flag', data=dq_flag, compression='gzip') 
    tmp = f.create_dataset('valid_data_flag', data=valid_data_flag, compression='gzip')
    tmp = f.create_dataset('saturate_pixel', data=saturate_pixel, compression='gzip')
    tmp = f.create_dataset('median_image', data=median_image, compression='gzip')
    # Now make many datasets from the header parameters
    for i in range(len(keepprihdr)):
        curval = hdulist[0].header[keepprihdr[i]]
        if np.isscalar(curval):
            tmp = f.create_dataset(keepprihdr[i], data=np.array([hdulist[0].header[keepprihdr[i]]], dtype=formatprihdr[i]))
        else:
            tmp = f.create_dataset(keepprihdr[i], data=np.array([-1], dtype=formatprihdr[i]))
            
    for i in range(len(keep1hdr)):
        curval = hdulist[1].header[keep1hdr[i]]
        if np.isscalar(curval):
            tmp = f.create_dataset(keep1hdr[i], data=np.array([hdulist[1].header[keep1hdr[i]]], dtype=format1hdr[i]))
        else:
            tmp = f.create_dataset(keep1hdr[i], data=np.array([-1], dtype=format1hdr[i]))
            
    f.close()

# parse a page/sector list
# from https://stackoverflow.com/a/4248689
def parse_range(astr):
    result=set()
    for part in astr.split(','):
        x=part.split('-')
        result.update(range(int(x[0]),int(x[-1])+1))
    return sorted(result)

if __name__ == "__main__":
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
    
    dirOutputs = '/pdo/users/cjburke/spocvet/sector56/'
    SECTOR = 56# =-1 if multi-sector
    RESAMP = 5  ###  USE AN ODD NUMBER HELPS WITH CADENCE NO ###
    overwrite = False

    if (not 'maxsec' in locals()) and SECTOR == -1:
        print('Code expects multi-sector (SECTOR = 56) but no sector argument was given. EXITING!')
        sys.exit()

    if SECTOR > 0:
        print('Single Sector')
        # In the case of a single sector One needs dummy entries for
        #  every sector
        fileInputPrefixList = []
        for i in np.arange(1,SECTOR):
            fileInputPrefixList.append('/foo{0:d}'.format(i))
        fileInputPrefixList.append('/pdo/spoc-data/sector-056/target-pixel/tess2022244194134-s0056-')
        fileInputSuffixList = []
        for i in np.arange(1,SECTOR):
            fileInputSuffixList.append('/foo{0:d}'.format(i))
        fileInputSuffixList.append('-0243-s_tp.fits.gz')
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
                    curStr = os.path.join(pathPrefixes[ia][0],'target-pixel','{0}-s{1:04d}-'.format(filePrefixes[ia][0],i))
                    fileInputPrefixList.append(curStr)
                    curStr = '-{0:04d}-s_tp.fits.gz'.format(fileSuffixNumbers[ia][0])
                    fileInputSuffixList.append(curStr)
                else:
                    fileInputPrefixList.append('/foo{0:d}'.format(i))
                    fileInputSuffixList.append('/foo{0:d}'.format(i))
            else:
                fileInputPrefixList.append('/foo{0:d}'.format(i))
                fileInputSuffixList.append('/foo{0:d}'.format(i))
    #print(fileInputPrefixList)
    #print(fileInputSuffixList)
        
    nSector = len(fileInputPrefixList)    

    # Only do tpfs for the targets with TCEs
    #  You can specify a multisector tce seed file because
    #   al that it uses is TIC.  If it exists it is made
    # Load the tce data h5
    tceSeedInFile = 'sector-56_20221023_tce.h5'
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)

    alltic = np.unique(np.array([x.epicId for x in all_tces], dtype=np.int64))

    cnt=0

    for curTic in alltic:
        print(cnt)
        fileLCInputList = glob.glob(os.path.join(make_data_dirs(dirOutputs, SECTOR, curTic), 'tess_dvts_{0:016d}_*.h5d'.format(curTic)))
        if len(fileLCInputList)>0:    
            # Sort list so that the first PC is used to determine data quality
            fileLCInputList.sort()
            for k in range(nSector):
                fileInput = '{0}{1:016d}{2}'.format(fileInputPrefixList[k], curTic, fileInputSuffixList[k])
    
                fileOutput = os.path.join(make_data_dirs(dirOutputs, SECTOR, curTic), 'tess_tpf_{0:016d}_{1:02d}.h5d'.format(curTic, k+1))
                if os.path.isfile(fileInput):
                    fileExists=os.path.isfile(fileOutput)
                    if (not fileExists) or overwrite:
                        tpf_resamp(fileInput, fileOutput, RESAMP, fileLCInputList[0])
                    else:
                        print('Skipping ', fileOutput)
        cnt = cnt + 1
