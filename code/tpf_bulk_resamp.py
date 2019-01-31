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

def tpf_resamp(file, dirOut, RESAMP):
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
    #Make a fix for Sector 3 where not all cadences were used
    # in the backend DV
    #idx = np.where((cadenceNo>=114115) & (cadenceNo<=128706))[0]
    #nImage = len(idx)
    #cadenceNo, timetbjd, dq_flag = idx_filter(idx, cadenceNo, timetbjd, dq_flag)
    #flux_array = flux_array[idx, :, :]
    #flux_bkg_array = flux_bkg_array[idx, :, :]

    newNImage = int(np.floor(nImage / RESAMP))
    oldNImage = newNImage*RESAMP
    # trim off the excess images not integral into resamp
    cadenceNo = cadenceNo[0:oldNImage]
    timetbjd = timetbjd[0:oldNImage]
    flux_array = flux_array[0:oldNImage, :, :]
    flux_bkg_array = flux_bkg_array[0:oldNImage, :,:]
    dq_flag = dq_flag[0:oldNImage]
    
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
    fileoutput = os.path.join(make_data_dirs(dirOut,sec,epic), 'tess_tpf_{0:016d}.h5d'.format(epic))
    f = h5py.File(fileoutput, 'w')
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

    dirInputs = '/pdo/spoc-data/sector-04/target-pixel'
    dirOutputs = '/pdo/users/cjburke/spocvet/sector4/'
    RESAMP = 5  ###  USE AN ODD NUMBER HELPS WITH CADENCE NO ###

    # Only do tpfs for the targets with TCEs
    #  You can specify a multisector tce seed file because
    #   al that it uses is TIC.  If it exists it is made
    tceSeedInFile = 'sector4_20190129_tce.pkl'
    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()
    alltic = np.unique(np.array([x.epicId for x in all_tces], dtype=np.int64))


    fileList = glob.glob(os.path.join(dirInputs, '*tp.fits.gz'))

    cnt=0
    for fil in fileList:
        dirName, fileName = os.path.split(fil)
        tokens = fileName.split('-')
        curTic = np.int64(tokens[2])
        idx = np.where(alltic == curTic)[0]
        if len(idx) > 0: 
            print(cnt)
            tpf_resamp(fil, dirOutputs, RESAMP)
            cnt = cnt + 1
