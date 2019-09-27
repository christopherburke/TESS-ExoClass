#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:51:04 2018
Calculate the skyline histogram of cadences that contribute to all
TCEs.  This is used to avoid cadences contribute strongly to detections
implying that something is wrong with this data to induce significantly
more detections than other cadences.

@author: Christopher J. Burke
"""

import numpy as np
import cjb_utils as cjb
import scipy.special as spec
import toidb_federate as fed
from gather_tce_fromdvxml import tce_seed
import pickle
import matplotlib.pyplot as plt
from statsmodels import robust
import os


def genericFed(per, epc, tryper, tryepc, trydur, trypn, trytic, tStart, tEnd):
    ts = fed.timeseries(tStart, tEnd)
    federateResult = fed.federateFunction(per, epc, ts, \
                    tryper, tryepc, trydur)
    bstpn = trypn[federateResult[0]]
    bsttic = trytic[federateResult[0]]
    bstMatch = int(federateResult[1])
    bstStat = federateResult[2]
    bstPeriodRatio = federateResult[3]
    bstPeriodRatioFlag = federateResult[4]
    bstFederateFlag = federateResult[5]
    
    return bstpn, bsttic, bstMatch, bstStat, bstPeriodRatio, bstPeriodRatioFlag, bstFederateFlag


if __name__ == '__main__':
    fout = open('skyline_data_sector14_20190918.txt', 'w')
    
    # Load the tce data pickle    
    tceSeedInFile = 'sector14_20190918_tce.pkl'
    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()
    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=np.int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=np.int)
    allrp = np.array([x.at_rp for x in all_tces])
    allper = np.array([x.at_period for x in all_tces])
    alldur = np.array([x.at_dur for x in all_tces])
    allepc = np.array([x.at_epochbtjd for x in all_tces])
    alltrpvalid = np.array([x.trp_valid for x in all_tces], dtype=np.int)
    alltrpdur = np.array([x.trp_dur for x in all_tces])
    alltrpepc = np.array([x.trp_epochbtjd for x in all_tces])
    alltcedur = np.array([x.pulsedur for x in all_tces])
    alltceepc = np.array([x.tce_epoch for x in all_tces])
    alltceper = np.array([x.tce_period for x in all_tces])
    # Go through each tce and use valid fits from dv, trpzd, tce in that order for matching
    useper = np.zeros_like(allper)
    useepc = np.zeros_like(allper)
    usedur = np.zeros_like(allper)
    idx = np.where((allatvalid == 0) & (alltrpvalid == 0) )[0]
    useper[idx] = alltceper[idx]
    useepc[idx] = alltceepc[idx]
    usedur[idx] = alltcedur[idx]
    idx = np.where((alltrpvalid == 1))[0]
    useper[idx] = alltceper[idx]
    useepc[idx] = alltrpepc[idx]
    usedur[idx] = alltrpdur[idx]
    idx = np.where((allatvalid == 1))[0]
    useper[idx] = allper[idx]
    useepc[idx] = allepc[idx]
    usedur[idx] = alldur[idx]

    uowStart = np.min(useepc) - 1.0
    uowEnd = np.max(useepc) + 13.0
    # Check if cadnoVtimemap.txt  cadence no to time mapping exists
    # created by dvts_bulk_resamp.py use it if so
    if os.path.exists('cadnoVtimemap.txt'):
        dataBlock = np.genfromtxt('cadnoVtimemap.txt', dtype=['i4','f8','i4','i4'])
        uowStart = np.min(dataBlock['f1'])
        uowEnd = np.max(dataBlock['f1'])
    print('Data Start: {0} Data End: {1}'.format(uowStart, uowEnd))
    # put an upper limit on duration
    idx = np.where(usedur>5.0)[0]
    usedur[idx] = 5.0

    ts = fed.timeseries(uowStart, uowEnd)
    skylineData = np.zeros_like(ts.ephemCentral, dtype=np.int64)
    
    print('alpha')

    for i in range(len(alltic)):
        curTic = alltic[i]
        curper = useper[i]
        curepc = useepc[i]
        curdur = usedur[i]
        ts = ts.makeEphemVector(curper,curepc,curdur)
        skylineData = skylineData + np.copy(ts.ephemFull)

    medSkyline = np.median(skylineData)
    madSkyline = robust.mad(skylineData)
    tmp = np.arange(len(skylineData))
    plt.plot(tmp, skylineData, '.')
    idxBad = np.where((skylineData-medSkyline)/madSkyline > 2.75)[0]
    plt.plot(tmp[idxBad], skylineData[idxBad], '.r')

    plt.show()

    for j in idxBad:
        fout.write('{:11.5f}\n'.format(ts.ts[j]))
    fout.close()

    
    print('hello world')