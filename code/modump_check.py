# -*- coding: utf-8 -*-
"""
Check for events falling near momentum dumps

AUTHOR: Christopher J. Burke
"""

import numpy as np
import pickle
from gather_tce_fromdvxml import tce_seed
import os
from subprocess import Popen, PIPE
import math
import h5py
from statsmodels import robust
from pgmcmc import pgmcmc_ioblk, pgmcmc_setup
from pgmcmc import pgmcmc_run_mcmc, pgmcmc_run_minimizer
import matplotlib.pyplot as plt

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


def phaseData(t, per, to):
    """Phase the data at period per and centered at to
       INPUT:
         t - time of data
         per - period to phase time data period and t should
               be in same units
         to - epoch of phase zero
       OUTPUT:
         phi - data phased running from -0.5<phi<=0.5
     """
    phi = np.mod(t - to, per) / per
    phi = np.where(phi > 0.5, phi - 1.0, phi)
    return phi

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list



if __name__ == '__main__':
    #  Directory storing the ses mes time series
    sesMesDir = '/pdo/users/cjburke/spocvet/sector51'
    SECTOR = 51

    
    fileOut = 'spoc_modump_sector-51_20220624.txt'    
    fom = open(fileOut, 'w')
    vetFile = 'spoc_fluxtriage_sector-51_20220624.txt'
    #vetFile = 'junk.txt'
    tceSeedInFile = 'sector-51_20220624_tce.h5'

    # cadence number time mapping has momentum dump flags in it
    # It is generated in dvts_bulk_resamp.py
    dataBlock = np.genfromtxt('cadnoVtimemap.txt', dtype=['i4','f8','i4','i4','i4'])
    cadmap = dataBlock['f0']
    timemap = dataBlock['f1']
    momdump = dataBlock['f3']
    idx = np.where(momdump == 1)[0]
    bdTime = timemap[idx]
    
    # Load the tce data h5
    tceSeedInFile = 'sector-51_20220624_tce.h5'
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
    # These lines can be used for debugging
    #idx = np.where((alltic == 101955023))[0]
    #alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
    #        allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
    #        allatrpdrstar, allatrpdrstare, allatadrstar = idx_filter(idx, \
    #        alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
    #        allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
    #        allatrpdrstar, allatrpdrstare, allatadrstar)
            
    # Recalculate MES after removing momentum dump impacted events
    for i, curTic in enumerate(alltic):
        print('{:d} of {:d}'.format(i, len(alltic)))
        curPn = allpn[i]
        curDur = alldur[i]
        curDurDay = curDur/24.0
    
        fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(curTic,curPn))
        f = h5py.File(fileInput,'r')
        allCorr = np.array(f['allCorr'])
        allNorm = np.array(f['allNorm'])
        allTime = np.array(f['allTime'])
        allCadNo = np.array(f['allCadNo'])
        nBd = 0
        for j in range(len(allCorr)):
            curTime = allTime[j]
            tDiff = np.abs(curTime - bdTime)
            idx = np.where(tDiff < curDurDay)[0]
            if len(idx)>0:
                nBd = nBd + 1
        fracBd = float(nBd)/float(len(allCorr))
        #print('{0:d} {1:f}'.format(curTic, fracBd))

        
        fom.write('{:016d} {:02d} {:f}\n'.format( \
                      curTic, curPn, fracBd))

    fom.close()