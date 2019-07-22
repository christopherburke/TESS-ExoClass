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
    
    # load the TEV TCE  data

    qlpfile = 'csv-file-2019-05-07.csv'
    dtypeseq = ['i4','i4','U30','U30','U30','U30']
    dataBlock = np.genfromtxt(qlpfile, \
                              dtype=dtypeseq, delimiter=',',skip_header=1)
    gtTIC = dataBlock['f0']
    gtPN = dataBlock['f1']

    # Load the tce data pickle    
    tceSeedInFile = 'sector12_20190712_tce.pkl'
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

    print('Total TCEs {0:d}'.format(len(alltic)))
    uniqTic = np.unique(alltic)
    print('Total Targets {0:d}'.format(len(uniqTic)))
    
    vetFile = 'spoc_fluxtriage_sector12_20190712.txt'
    # Load the  flux vetting
    dataBlock = np.genfromtxt(vetFile, dtype=[int,int,int,'S1'])
    fvtic = dataBlock['f0']
    fvpn = dataBlock['f1']
    fvvet = dataBlock['f2']

    nTot = 0
    nMiss = 0
    for i in range(len(gtTIC)):
        curTic = gtTIC[i]
        curPN = gtPN[i]
        idx = np.where(fvtic == curTic)[0]
        if len(idx)>0:
            prat = np.array([], dtype=np.float)
            passTriage = False
            nTot += 1
            for j in idx:
                curPn = fvpn[j]
                if fvvet[j] == 1:
                    passTriage = True
            if not passTriage:
                print('Not Passing triage {0:010d} {1:d}'.format( curTic, curPN))
                nMiss += 1
    print('{0:d} TOIs had TICS with TCEs'.format(nTot))
    print('Triage has {0:f} Recall Rate'.format((nTot-nMiss)/float(nTot)))
    idx = np.where(fvvet == 1)[0]
    triTic = fvtic[idx]
    triPn = fvpn[idx]
    print('\nTCEs pass triage: {0:d}'.format(len(triTic)))
    triUniqTic = np.unique(triTic)
    print('On {0:d} targets'.format(len(triUniqTic)))

