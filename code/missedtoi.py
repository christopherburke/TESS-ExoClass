#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:51:04 2018

@author: cjburke
"""

import pickle
import numpy as np
from gather_tce_fromdvxml import tce_seed
import csv
import sys
import time
import json
import cjb_utils as cjb




if __name__ == '__main__':
    
    # load the TOI data
    qlpfile = 'hlsp_tess-data-alerts_tess_phot_alert-summary-s01+s02+s03+s04_tess_v9_spoc.csv'
    dtypeseq = ['i4','f8','U2']
    dtypeseq.extend(['f8']*10)
    dataBlock = np.genfromtxt(qlpfile, \
                              dtype=dtypeseq, delimiter=',',skip_header=1)
    gtTIC = dataBlock['f0']
    gtTOI = dataBlock['f1']
    gtDisp = dataBlock['f2']
    gtPer = dataBlock['f9']
    gtEpc = dataBlock['f7']
    gtDur = dataBlock['f11']
    # Need to fix single transit TOIs with nan for period
    idx = np.where(np.logical_not(np.isfinite(gtPer)))[0]
    gtPer[idx] = 1000.0

    # Load the tce data pickle    
    tceSeedInFile = 'sector5_20190223_tce.pkl'
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
#    ia = np.argsort(alltic)
#    alltic, allpn, useper, useepc, usedur = cjb.idx_filter(ia, alltic, \
#                                allpn, useper, useepc, usedur)

    vetFile = 'spoc_sector5_fluxtriage_20190223.txt'
    # Load the  flux vetting
    dataBlock = np.genfromtxt(vetFile, dtype=[int,int,int,'S1'])
    fvtic = dataBlock['f0']
    fvpn = dataBlock['f1']
    fvvet = dataBlock['f2']

    for i in range(len(gtTIC)):
        curTic = gtTIC[i]
        curper = gtPer[i]
        curepc = gtEpc[i]
        curToi = gtTOI[i]
        idx = np.where(fvtic == curTic)[0]
        if len(idx)>0:
            prat = np.array([], dtype=np.float)
            passTriage = False
            for j in idx:
                curPn = fvpn[j]
                if fvvet[j] == 1:
                    passTriage = True
                    idx2 = np.where((alltic == curTic) & (allpn == curPn))[0]
                    tcePer = useper[idx2]
                    prat = np.append(prat, tcePer/curper)
            if not passTriage:
                print('{0:5.2f} not Passing triage {1:010d} {2}'.format(curToi, curTic, curper))
                