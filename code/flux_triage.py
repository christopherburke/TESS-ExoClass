#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform the Flux Triage stage of TEC


@author: Christopher J. Burke (MIT)
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import pickle
import math
from gather_tce_fromdvxml import tce_seed
import scipy.special as spec

def coughlin_sigmap(p1,p2):
    up1 = p1
    up2 = p2
    if up1 > up2:
        tmp = up2
        up2 = up1
        up1 = tmp
    delP = (up1 - up2)/up1
    delPp = abs(delP - round(delP))
    return np.sqrt(2.0)*spec.erfcinv(delPp)
    

def get_useable_ephems(all_tces):
    allepics = np.zeros((len(all_tces),), dtype=np.int64)
    allpns = np.zeros((len(all_tces),), dtype=np.int)
    allper = np.zeros((len(all_tces),))
    allepoch = np.zeros_like(allper)
    allduration = np.zeros_like(allper)
    for i,td in enumerate(all_tces):
        allepics[i] = td.epicId
        allpns[i] = td.planetNum
        if td.at_valid == 1:
            allper[i] = td.at_period
            allepoch[i] = td.at_epochbtjd
            allduration[i] = td.at_dur
        elif td.trp_valid ==1:
            allper[i] = td.tce_period
            allepoch[i] = td.trp_epochbtjd
            allduration[i] = td.trp_dur
        else:
            allper[i] = td.tce_period
            allepoch[i] = td.tce_epoch
            allduration[i] = td.pulsedur
    return allepics, allpns, allper, allepoch, allduration

if __name__ == '__main__':
    
    # Load the pickle file that contains TCE seed information
    # The pickle file is created by gather_tce_fromdvxml.py
    tceSeedInFile = 'sector1-13_20190812_tce.pkl'

    #  Directory storing the ses mes data
    sesDataDir = '/pdo/users/cjburke/spocvet/sector1-13'
    SECTOR = -1
    fluxVetOut = 'spoc_fluxtriage_sector1-13_20190812.txt'
#    fluxVetOut = 'junk.txt'

    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()
    
    allepics, allpns, allpers, allepochs, alldurations = get_useable_ephems(all_tces)
    
    fout = open(fluxVetOut, 'w')
    debug = False
    # Loop over tces and perform flux vetting
    # For a particular id
    #debug=True
    #alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    #idxdebug = np.where(alltic == 123702439)[0]
    cnt = 0
    #for td in [all_tces[idxdebug[0]],all_tces[idxdebug[0]]]:
    for td in all_tces:
        print(cnt)
        cnt = cnt+1
        epicid = td.epicId
        pn = td.planetNum
        period = 0.0
        epoch = 0.0
        duration = 0.0
        if td.at_valid == 1:
            period = td.at_period
            epoch = td.at_epochbtjd
            duration = td.at_dur
        elif td.trp_valid == 1:
            period = td.tce_period
            epoch = td.trp_epochbtjd
            duration = td.trp_dur
        else:
            period = td.tce_period
            epoch = td.tce_epoch
            duration = td.pulsedur
        searchDurationHours = td.pulsedur
        origMes = td.mes
        print('Orig ',origMes, epicid, pn)
        epcDir = '{0:04d}'.format(int(math.floor(epicid/1000.0)))
        localDir = os.path.join(sesDataDir,'S{0:02d}'.format(SECTOR),epcDir)
        fileInput = os.path.join(localDir, 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(epicid,pn))
        f = h5py.File(fileInput,'r')

        if debug:
            phi = np.array(f['phi'])
            eve = np.array(f['events'])
            tphi = phi+eve
            vd = np.array(f['validData'])
            flx = np.array(f['initFlux'])
            altflx = np.array(f['altDetrend'])
            normTS = np.array(f['normTS'])
            corrTS = np.array(f['corrTS'])
            sesTS = corrTS/normTS
            cdppTS = 1.0e6/normTS
            
            vd_ext = np.array(f['valid_data_flag_ext'])
            plt.plot(tphi, flx[vd], '.')
            plt.show()
            plt.plot(tphi, altflx[vd], '.')
            plt.show()
            plt.plot(tphi, sesTS[vd_ext], '.' )
            plt.show()
            plt.plot(tphi, cdppTS[vd_ext], '.')
            plt.show()
            
        # First check that the detrending worked
        flux_pass = True
        flux_str = ['']
        validAltDet = np.array(f['validAltDet'])[0]
        if validAltDet == 0:
            flux_pass = False
            flux_str.append('AltDetFail')
        # Check for ses stats worked
        if flux_pass:
            validSes = np.array(f['validSes'])[0]
            if validSes == 0:
                flux_pass = False
                flux_str.append('SesStatFail')
        # Check for new mes below thresh
        if flux_pass:
            newMes = np.array(f['newMes'])[0]
            if newMes < 4.5:
                flux_pass = False
                flux_str.append('newMesBelowThresh')
        # Check for new mes shrink fractionally too much
        if flux_pass:
            lgOrigMes = np.log10(origMes)
            delMesCutLine = (lgOrigMes-1.0)*(0.76-0.5)/(2.0-1.0) + 0.53
            delMes = (origMes - newMes) / origMes
            if delMes > delMesCutLine:
                flux_pass = False
                flux_str.append('newMesShrink')
        # Check for ses2Mes
        if flux_pass:
            mnSes2mnMes = np.array(f['mnSes2mnMes'])[0]
            mnSes2mnMes_r = np.array(f['mnSes2mnMes_r'])[0]
            nTran = np.array(f['newNTran'])[0]
            #print(origMes, mnSes2mnMes, mnSes2mnMes_r, nTran)
            if (origMes<=13.0 and mnSes2mnMes>0.97) or (origMes>13.0 and mnSes2mnMes>0.85):
                if (mnSes2mnMes>mnSes2mnMes_r*1.1):
                    flux_pass = False
                    flux_str.append('ses2MesFail')
        # Check for Chases
        if flux_pass:
            chasesSumry = np.array(f['chasesSumry'])[0]
            singSes = newMes / np.sqrt(nTran)
            print(singSes, chasesSumry, newMes, origMes, nTran)
            if singSes > 2.8 and chasesSumry < 0.3:
                flux_pass = False
                # Do an override for very strong signals that the shoulders throw a false chases fail
                if singSes > 20.0 and nTran >= 4 and chasesSumry > 0.1:
                    flux_pass = True
                else:
                    flux_str.append('ChasesFail')
                
        # Check for TCE being a secondary of previous TCE
        if flux_pass:
            if pn >= 2:
                idx = np.where(allepics == epicid)[0]
                tmppns = allpns[idx]
                tmppers = allpers[idx]
                tmpepochs = allepochs[idx]
                tmpdurations = alldurations[idx]
                haveMatch = False
                for jj in range(len(tmppns)):
                    if (tmppns[jj] < pn) and (not haveMatch):
                        sigp = coughlin_sigmap(period, tmppers[jj])
                        if sigp > 2.9:
                            haveMatch = True
                            flux_pass = False
                            flux_str.append('SecondaryOfPN_{:02d}'.format(tmppns[jj]))
        if flux_pass:            
            str = '{0:9d} {1:2d} {2:1d} PASS\n'.format(epicid, pn, int(flux_pass))
        else:
            str = '{0:9d} {1:2d} {2:1d} {3}\n'.format(epicid, pn, int(flux_pass), flux_str[1])
        if debug:
            print(str)
        else:
            fout.write(str)
        
    fout.close()