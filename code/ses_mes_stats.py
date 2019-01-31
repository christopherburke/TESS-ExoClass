#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 09:25:25 2018
This routine performs detrending of the flux time series (see smoothn.py)
The calculates the single event statistics (SES) using a python port
of the Kepler pipeline TPS wavelet search methodology.  The SES timeseries
allows recalculating the transit SNR (or MES in TPS parlance) under alternative
detrending and allows metrics to be calculated to determine whether the
signal is consistent a transit and not a systematic.
In particular from the Kepler robovetter the SES/MES test is implemented
here and the CHASES test.  


@author: Christopher J. Burke
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import pickle
import os
import math
import fluxts_conditioning as flux_cond
import kep_wavelets as kw
from gather_tce_fromdvxml import tce_seed
import scipy.stats as st
from statsmodels import robust

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

def assignEvents(t, epc, phi, per, phiDur):
    """Phase the data at period per and centered at to
       INPUT:
         phi - data phase running from -0.5<phi<=0.5
       OUTPUT:
         events - transit event number first appearing in data is event 1
     """
    events = np.zeros_like(phi)
    idx = np.where(np.abs(phi) <= phiDur)[0]
    if len(idx)>0:
        runcad = np.arange(len(t))
        minidx = np.min(idx)
        minto = t[minidx]
        idx2 = np.where(np.abs(t-minto) < 2.0*phiDur*per)[0]
        tmpcad = runcad[idx2]
        tmpphi = phi[idx2]
        ia = np.argmin(np.abs(tmpphi))
        cntrto = t[tmpcad[ia]]
        oldeve = np.round((cntrto-epc)/per)
        newcntrto = epc + oldeve*per
        events = np.round((t-newcntrto)/per) + 1.0
    return events

def get_ses_stats(corr, norm, phi, phiDur, events, debug=True):
    chasesWindowFac = 6
    chasesPeakFac = 0.7
    all_ses = np.array([])
    all_chases = np.array([])
    all_corr = np.array([])
    all_norm = np.array([])
    runcad = np.arange(len(corr))
    idx1 = np.where(np.abs(phi) <= 1.0*phiDur)[0]
    wideLimit = np.min([chasesWindowFac*phiDur, 0.2])

    tmpcad = runcad[idx1]
    idxGaps = np.where(np.diff(tmpcad)>1)[0]
    if len(idxGaps) > 0:
        idxStrts = tmpcad[idxGaps+1]
        idxStrts = np.insert(idxStrts, 0, tmpcad[0])
        idxEnds = tmpcad[idxGaps]+1
        idxEnds = np.append(idxEnds, tmpcad[-1] + 1)
        for i in range(len(idxStrts)):
            js = idxStrts[i]
            je = idxEnds[i]
            curcorr = corr[js:je]
            curnorm = norm[js:je]
            curcad = runcad[js:je]
            curevents = events[js:je]
            ia = np.argmax(curcorr/curnorm)
            all_corr = np.append(all_corr, curcorr[ia])
            all_norm = np.append(all_norm, curnorm[ia])
            maxSes = curcorr[ia]/curnorm[ia]
            all_ses = np.append(all_ses, curcorr[ia]/curnorm[ia])
            if debug:
                tmpx = np.arange(len(curcorr))
                plt.plot(tmpx, curcorr/curnorm, '.')
                plt.plot(tmpx[ia], curcorr[ia]/curnorm[ia], '.')
                plt.show()
            # Do Chases
            cureve = curevents[ia]
            idx2 = np.where((events == cureve) & (np.abs(phi)<wideLimit))[0]
            usePhi = phi[idx2]
            useSes = corr[idx2]/norm[idx2]
            xxx = np.arange(len(usePhi))
            oUseSes = np.copy(useSes)
            idx3 = np.where((np.abs(usePhi) > phiDur) | ((np.abs(usePhi)<phiDur) & (useSes<0.0)))[0]
            if len(idx3) > 0:
                usePhi = usePhi[idx3]
                useSes = useSes[idx3]
                closePhi = np.min(np.abs(usePhi))
                if debug:
                    plt.plot(xxx, oUseSes, '.')
                    plt.plot(xxx[idx3], useSes, '.')
                    plt.show()
                idx4 = np.where(np.abs(useSes) > chasesPeakFac*maxSes)[0]
                chase_val = 1.0
                if len(idx4) > 0:
                    nexPhi = np.min(np.abs(usePhi[idx4]))
                    chase_val = (nexPhi-closePhi)/wideLimit
            else:
                chase_val = 0.0
            all_chases = np.append(all_chases, chase_val)
            
            
    else:
        # only a single transit event in data set this data bad
        all_ses = np.append(all_ses, 0.0)
        all_chases = np.append(all_chases, 0.0)
        all_corr = np.append(all_corr, 0.0)
        all_norm = np.append(all_norm, 1.0)
        
    # Recalculate the mes
    sumNumer = np.sum(all_corr)
    sumDenum = np.sqrt(np.sum(all_norm*all_norm))
    newMes = sumNumer/sumDenum
    # Calculate ses to mes  and noise average ses to mes
    newNtran = len(all_ses)
    if newNtran > 1:
        ses2Mes = np.max(all_ses) / newMes
        mnNorm = np.mean(norm)
        depthEst = all_ses/all_norm
        all_mnses = depthEst * mnNorm
        all_mncorr = all_mnses * mnNorm
        sumNumer = np.sum(all_mncorr)
        sumDenum = np.sqrt(mnNorm*mnNorm*len(all_mncorr))
        newMnMes = sumNumer/sumDenum
        mnSes2mnMes = np.max(all_mnses) / newMnMes
    else:
        ses2Mes = 3.0
        newMnMes = newMes
        mnSes2mnMes = 3.0
    # Calculate chases
    if newNtran>2 & newNtran<5:
        chases_sum = np.median(all_chases)
    elif newNtran>=5:
        chases_sum = np.mean(all_chases)
    elif newNtran==2:
        chases_sum = np.min(all_chases)
    else:
        chases_sum = 0.0

    print('new ', newMes, ses2Mes, newNtran, chases_sum, newMes/np.sqrt(newNtran), newMnMes, mnSes2mnMes)

    return newMes, ses2Mes, newNtran, chases_sum, all_ses, all_chases, newMnMes, mnSes2mnMes
    
if __name__ == "__main__":
    # These are for parallel procoessing
    wID = 5
    nWrk = 6
    # Load the pickle file that contains TCE seed information
    # The pickle file is created by gather_tce_fromdvxml.py
    tceSeedInFile = 'sector4_20190129_tce.pkl'
    #  Directory storing the resampled dv time series data
    dvDataDir = '/pdo/users/cjburke/spocvet/sector4'
    # Directory of output hd5 files
    outputDir = dvDataDir
    SECTOR=4
    # What fraction of data can be missing and still calculat ses_mes
    # In Sector 1 due to the 2 days of missing stuff it was 0.68
    validFrac = 0.52
    overWrite = False

    # Skyline data excises loud cadecnes
    dataBlock = np.genfromtxt('skyline_data_sector4_20190129.txt', dtype=['f8'])
    badTimes = dataBlock['f0']

    # Search and filter parameters
    cadPerHr = 6
    firstFilterScaleFac = 10 # low frequency median filter will be
                            # firstFilterScaleFac*searchDurationHours medfilt window

    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()

    # These next few lines can be used to examine a single target    
    #all_epics = np.array([x.epicId for x in all_tces], dtype=np.int64)
    #all_pns = np.array([x.planetNum for x in all_tces], dtype=np.int)
    #ia = np.where((all_epics == 31850842) & (all_pns == 3))[0]
    #doDebug = True
    # Loop over tces and perform various ses, mes, chases tests
    cnt = 0
    doDebug = False
    # This for loop can be used for debugging
    #for td in all_tces[ia[0]:ia[0]+1]:
    # Normal for loop
    for td in all_tces:
        epicid = td.epicId
        print(cnt, epicid)
        cnt = cnt+1
        if np.mod(cnt, nWrk) == wID:
            if np.mod(cnt,10) == 0:
                print(cnt)
    
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
            fileOutput = os.path.join(make_data_dirs(outputDir, SECTOR, epicid), 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(epicid,pn))
            fileExists = os.path.isfile(fileOutput)
            if (not fileExists) or overWrite:
                print('pulse: {:f} fitdur: {:f}'.format(td.pulsedur, duration))
                searchDurationHours = np.max([td.pulsedur, duration])
                # Cap duration at 15 hours
                searchDurationHours = np.min([searchDurationHours, 15.0])
        
                origMes = td.mes
                print('Orig ',origMes, epicid, pn)
                localDir = make_data_dirs(dvDataDir,SECTOR,epicid)
                fileInput = os.path.join(localDir, 'tess_dvts_{0:016d}_{1:02d}.h5d'.format(epicid,pn))
                f = h5py.File(fileInput,'r')
                # Decide which flux time series to use pdc or lc_init
                # ***NOTE harmonic filter is not being used so pdf flux should be used
                #useFlux = np.array(f['lc_init'])
                #  ***If flux is norm subtracted at 1.0 to it
                #useFlux = useFlux + 1.0
                useFlux = np.array(f['pdc_flux'])
                #useFlux3 = np.array(f['lc_white'])
                #useFlux3 = useFlux3 + 1.0
                vd = np.array(f['valid_data_flag'])
                time = np.array(f['timetbjd'])
                for jj, curTime in enumerate(time):
                    diff = np.abs(curTime - badTimes)
                    if np.min(diff) < 5.0/60.0/24.0:
                        vd[jj] = False
                #if doDebug:
                #    plt.plot(useFlux, '.')
                #    plt.show()
                
                # Mark data in transit as deweighted during detrending
                ootvd = np.full_like(vd, True)
                phi = phaseData(time, period, epoch)
                idx = np.where(np.abs(phi) < duration/2.0/24.0/period)[0]
                ootvd[idx] = False
                
                # For ETE6 only there was unexpectedly extra data after 
                #  indice 1348. Mark those aftr 1348 as invalid now
                #vd[1349:] = False
                
                # For Sector 1 data there was a crappy part 
                #vd[2465:2750] = False
                # detrend data fixe edges
                # Verify that there is enough data to do analysis
                # Also do a period cut
         
                vdfrac = float(len(np.where(vd)[0]))/len(vd)
                if len(np.where(vd)[0]) > 200 and period > 0.3 and duration/24.0/period < 0.2 and vdfrac>validFrac:
                    print('Valid fraction: {:f}'.format(vdfrac))
                    # Look for large differences due to harmonic filter being applied.
                    # ***NOTE harmonic filter is not being used for TESS***
                    #tmpcad = np.arange(len(vd))
                    #idx = np.where(vd)[0]
                    #tmpcad = tmpcad[idx]
                    #tmpflux = useFlux[idx]
                    #linfit = st.linregress(tmpcad, tmpflux)
                    #tmp_linflat = tmpflux / (linfit[0]*tmpcad + linfit[1])
                    #tmpflux = useFlux2[idx]
                    #linfit = st.linregress(tmpcad, tmpflux)
                    #tmp_linflat2 = tmpflux / (linfit[0]*tmpcad + linfit[1])
                    #tmpflux = useFlux3[idx]
                    #linfit = st.linregress(tmpcad, tmpflux)
                    #tmp_linflat3 = tmpflux / (linfit[0]*tmpcad + linfit[1])
                    
                    #plt.plot(tmpcad, tmp_linflat, '.')
                    #plt.show()
                    #plt.plot(tmpcad, tmp_linflat2, '.')
                    #plt.show()
                    #plt.plot(tmpcad, tmp_linflat3, '.')
                    #plt.show()
                    #tmpdiff = tmp_linflat - tmp_linflat2
                    #print(robust.mad(tmpdiff))
        
                    #plt.plot(tmpcad, tmpdiff, '.')
                    #plt.show()
                    #tmpdiff = tmp_linflat - tmp_linflat3
                    #print(robust.mad(tmpdiff))
                    #plt.plot(tmpcad, tmpdiff, '.')
                    #plt.show()
        
                    final_smooth_flux, bad_edge_flag = flux_cond.detrend_with_smoothn_edgefix(\
                                        useFlux, vd, ootvd, int(np.ceil(cadPerHr*searchDurationHours)), fixEdge=True, \
                                         medfiltScaleFac=10, gapThreshold=5, edgeExamWindow=8, \
                                         edgeSig=6.0, edgeMinCad=50, debug=doDebug)
                    # fill gaps and extend to power of two
                    fillWindow = int(np.ceil(cadPerHr*searchDurationHours*firstFilterScaleFac*3))
                    final_smooth_flux_ext, vd_ext = flux_cond.fill_extend_fluxts(final_smooth_flux, \
                                                                vd, fillWindow, doExtend=True, debug=False)
                    #plt.cla()
                    #plt.plot(final_smooth_flux_ext, '.')
                    #if doDebug:
                    #    plt.show()
                    #else:
                    #    plt.pause(0.002)
                    waveletLen = 12
                    searchLen = np.int(np.round(cadPerHr * searchDurationHours))
                    varianceFilterFactor = 15
                    varianceFilterWindow = searchLen * varianceFilterFactor
                    wavObj = kw.waveletObject(waveletLen, final_smooth_flux_ext, varianceFilterWindow)
                    trial_pulse = kw.set_trial_transit_pulse(searchLen)
            
                    normTS, corrTS = kw.compute_statistic_time_series(wavObj, searchLen, trial_pulse)
                    #plt.cla()
                    #plt.plot(corrTS/normTS, '.')
                    #if doDebug:
                    #    plt.show()
                    #else:
                    #    plt.pause(0.002)
                    #plt.cla()
                    #plt.plot(1.0e6/normTS, '.')
                    #if doDebug:
                    #    plt.show()
                    #else:
                    #    plt.pause(0.002)
                    useNormTS = normTS[vd_ext]
                    useCorrTS = corrTS[vd_ext]
                    usePhase = phaseData(time[vd], period, epoch)
                    phaseDur = duration / 24.0 / period
                    useEvents = assignEvents(time[vd], epoch, usePhase, period, phaseDur)
                    newMes, ses2Mes, newNTran, Chases_sumry, allSes, allChases, \
                            newMnMes, mnSes2mnMes = get_ses_stats(useCorrTS, useNormTS,\
                                                usePhase, phaseDur, useEvents, debug=False)
                    if newNTran > 1:
                        validSes = 1
                    else:
                        validSes = 0
                    validAltDet = 1
                else:
                    print('Too little data for analysis epic: {0:d} pn: {1:d}'.format(epicid, pn))
                    validAltDet = 0
                    validSes = 0
                    newMes = 0.0
                    ses2Mes = 3.0
                    newNTran = 1
                    Chases_sumry = 0.0
                    allSes = np.array([0.0])
                    allChases = np.array([0.0])
                    final_smooth_flux = np.array([0.0])
                    bad_edge_flag = np.array([0.0])
                    vd_ext = np.array([0.0])
                    final_smooth_flux_ext = np.array([0.0])
                    normTS = np.array([0.0])
                    corrTS = np.array([0.0])
                    newMnMes = 0.0
                    mnSes2mnMes = 3.0
                    usePhase = np.array([0.0])
                    useEvents = np.array([0.0])
                    phaseDur = 0.0
                    
                #print(fileOutput)
                f = h5py.File(fileOutput,'w')
                tmp = f.create_dataset('altDetrend', data=final_smooth_flux, compression='gzip')
                tmp = f.create_dataset('validData', data=vd, compression='gzip')
                tmp = f.create_dataset('time', data=time, compression='gzip')
                tmp = f.create_dataset('phi', data=usePhase, compression='gzip')
                tmp = f.create_dataset('events', data=useEvents, compression='gzip')
                tmp = f.create_dataset('initFlux', data=useFlux, compression='gzip')
                tmp = f.create_dataset('bad_edge_flag', data=bad_edge_flag, compression='gzip')
                tmp = f.create_dataset('altDetrend_ext', data=final_smooth_flux_ext, compression='gzip')
                tmp = f.create_dataset('valid_data_flag_ext', data=vd_ext, compression='gzip')
                tmp = f.create_dataset('normTS', data=normTS, compression='gzip')
                tmp = f.create_dataset('corrTS', data=corrTS, compression='gzip')
                tmp = f.create_dataset('validSes', data=np.array([validSes], dtype=np.int))
                tmp = f.create_dataset('newMes', data=np.array([newMes], dtype=np.float))
                tmp = f.create_dataset('validAltDet', data=np.array([validAltDet], dtype=np.int))
                tmp = f.create_dataset('ses2Mes', data=np.array([ses2Mes], dtype=np.float))
                tmp = f.create_dataset('newNTran', data=np.array([newNTran], dtype=np.int))
                tmp = f.create_dataset('chasesSumry', data=np.array([Chases_sumry], dtype=np.float))
                tmp = f.create_dataset('allSes', data=allSes)
                tmp = f.create_dataset('allChases', data=allChases)
                tmp = f.create_dataset('newMnMes', data=np.array([newMnMes], dtype=np.float))
                tmp = f.create_dataset('mnSes2mnMes', data=np.array([mnSes2mnMes], dtype=np.float))
                tmp = f.create_dataset('phaseDur', data=np.array([phaseDur], dtype=np.float))
                f.close()
                