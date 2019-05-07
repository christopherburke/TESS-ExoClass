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
import pickle
from gather_tce_fromdvxml import tce_seed
import cjb_utils as cjb


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
    saturate_pixel = np.zeros((nx, ny), dtype=np.int)
    median_image = np.zeros((nx, ny))

    # Get header information that we should keep
    keepprihdr = ['TICID','SECTOR','CAMERA','CCD','PXTABLE','RA_OBJ', \
                  'DEC_OBJ','PMRA','PMDEC','PMTOTAL','TESSMAG','TEFF', \
                  'LOGG','RADIUS']
    formatprihdr = [np.uint32, np.int, np.int, np.int,np.int, \
                    np.float,np.float,np.float,np.float,np.float, \
                    np.float,np.float,np.float,np.float,np.float]

    keep1hdr = ['1CRV4P','2CRV4P','1CRPX4','2CRPX4']
    format1hdr = [np.int,np.int,np.float,np.float]
    
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
    frstIdx = ib[0]
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
    cadenceNo = np.mean(np.reshape(cadenceNo, (newNImage, RESAMP)), axis=1, dtype=np.int)
    timetbjd = np.mean(np.reshape(timetbjd, (newNImage, RESAMP)), axis=1)
    flux_array = np.sum(np.reshape(flux_array, (newNImage, RESAMP, nx, ny)), axis=1)
    flux_bkg_array = np.sum(np.reshape(flux_bkg_array, (newNImage, RESAMP, nx, ny)), axis=1)
    dq_flag = np.sum(np.reshape(dq_flag, (newNImage, RESAMP)), axis=1, dtype=np.int)

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
    

if __name__ == "__main__":
    #  Directory list for Sector light curve files
    fileInputPrefixList = ['/pdo/spoc-data/sector-01/sector_early_look/target-pixel/tess2018206045859-s0001-', \
                          '/pdo/spoc-data/sector-02/target-pixel/tess2018234235059-s0002-', \
                          '/pdo/spoc-data/sector-03/target-pixel/tess2018263035959-s0003-', \
                          '/pdo/spoc-data/sector-04/target-pixel/tess2018292075959-s0004-', \
                          '/pdo/spoc-data/sector-05/target-pixel/tess2018319095959-s0005-', \
                          '/pdo/spoc-data/sector-06/target-pixel/tess2018349182459-s0006-']
    fileInputSuffixList = ['-0120-s_tp.fits.gz', \
                           '-0121-s_tp.fits.gz', \
                           '-0123-s_tp.fits.gz', \
                           '-0124-s_tp.fits.gz', \
                           '-0125-s_tp.fits.gz', \
                           '-0126-s_tp.fits.gz']
# In the case of a single sector One needs dummy entries for
#  every sector
#    fileInputPrefixList = ['/foo1','/foo2','/foo3','/foo4','/foo5',\
#                           '/f006',\
#                           '/pdo/spoc-data/sector-01-06/target-pixel/tess2019032160000-s0008-']
#    fileInputSuffixList = ['-foo1','-foo2','-foo3','-foo4','-foo5',\
#                           '/f006',\
#                           '-0136-s_tp.fits.gz']

    nSector = len(fileInputPrefixList)    
    dirOutputs = '/pdo/users/cjburke/spocvet/sector1-6/'
    SECTOR = -1# =-1 if multi-sector
    RESAMP = 5  ###  USE AN ODD NUMBER HELPS WITH CADENCE NO ###

    # Only do tpfs for the targets with TCEs
    #  You can specify a multisector tce seed file because
    #   al that it uses is TIC.  If it exists it is made
    tceSeedInFile = 'sector1-6_20190428_tce.pkl'
    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()
    alltic = np.unique(np.array([x.epicId for x in all_tces], dtype=np.int64))


    cnt=0
    for curTic in alltic:
        print(cnt)
        fileLCInputList = glob.glob(os.path.join(make_data_dirs(dirOutputs, SECTOR, curTic), 'tess_dvts_{0:016d}_*.h5d'.format(curTic)))
        if len(fileLCInputList)>0:    
            for k in range(nSector):
                fileInput = '{0}{1:016d}{2}'.format(fileInputPrefixList[k], curTic, fileInputSuffixList[k])
    
                fileOutput = os.path.join(make_data_dirs(dirOutputs, SECTOR, curTic), 'tess_tpf_{0:016d}_{1:02d}.h5d'.format(curTic, k+1))
                if os.path.isfile(fileInput):
                    tpf_resamp(fileInput, fileOutput, RESAMP, fileLCInputList[0])
        cnt = cnt + 1
