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
import os
import math
import fluxts_conditioning as flux_cond
import kep_wavelets as kw
from gather_tce_fromdvxml import tce_seed
import scipy.stats as st
from statsmodels import robust
import argparse

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

def trapezoid(t, depth, bigT, littleT):
    """Trapezoid shape for model
       INPUT:
       t -  [float] vector of independent values to evaluate
                    trapezoid model
       depth - [float] depth of trapezoid
       bigT - [float] full trapezoid duration
       littleT - [float] 'ingress/egress' duration
       OUTPUT:
       output - [float] vector of trapezoid model values
    """
    output = np.full_like(t, 1.0)
    t = np.abs(t)
#    output = np.where(t <= bigT/2.0 - littleT/2.0, 1.0 - depth, output)
#    output = np.where(np.logical_and(t > bigT/2.0 - littleT/2.0, \
#                      t < bigT/2.0 + littleT/2.0),  \
#                      1.0 - depth + ((depth/littleT)* \
#                      (t-bigT/2.0 + littleT/2.0)), output)
    output = np.where(t <= littleT/2.0, 1.0 - depth, output)
    output = np.where(np.logical_and(t > littleT/2.0, t<=bigT/2.0), \
                      1.0-depth + ((depth/((bigT-littleT)/2.0))*(t-littleT/2.0)), output)
    return output

def get_ses_stats(corr, norm, corr_r, norm_r, phi, phiDur, events, time, oCadNo, origMes, debug=True):
    chasesWindowFac = 6
    chasesPeakFac = 0.7
    all_ses = np.array([])
    all_chases = np.array([])
    all_corr = np.array([])
    all_norm = np.array([])
    all_ses_r = np.array([])
    all_chases_r = np.array([])
    all_corr_r = np.array([])
    all_norm_r = np.array([])
    all_time = np.array([])
    all_cadNo = np.array([], dtype=np.int)
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
        chsThrSumOld = 0.0
        chsThrSumNew = 0.0
        
        for i in range(len(idxStrts)):
            js = idxStrts[i]
            je = idxEnds[i]
            curcorr = corr[js:je]
            curnorm = norm[js:je]
            curcadNo = oCadNo[js:je]
            curTime = time[js:je]
            curevents = events[js:je]
            curcorr_r = corr_r[js:je]
            curnorm_r = norm_r[js:je]
            ia = np.argmax(curcorr/curnorm)
            all_corr = np.append(all_corr, curcorr[ia])
            all_norm = np.append(all_norm, curnorm[ia])
            all_time = np.append(all_time, curTime[ia])
            all_cadNo = np.append(all_cadNo, curcadNo[ia])
            maxSes = curcorr[ia]/curnorm[ia]
            all_ses = np.append(all_ses, curcorr[ia]/curnorm[ia])
            ib = np.argmax(curcorr_r/curnorm_r)
            all_corr_r = np.append(all_corr_r, curcorr_r[ib])
            all_norm_r = np.append(all_norm_r, curnorm_r[ib])
            maxSes_r = curcorr_r[ib]/curnorm_r[ib]
            all_ses_r = np.append(all_ses_r, curcorr_r[ib]/curnorm_r[ib])
            if debug:
                tmpx = np.arange(len(curcorr))
                plt.plot(tmpx, curcorr/curnorm, '.')
                plt.plot(tmpx[ia], curcorr[ia]/curnorm[ia], '.')
                plt.plot(tmpx, curcorr_r/curnorm_r, '.')
                plt.plot(tmpx[ib], curcorr_r[ib]/curnorm_r[ib], '.')
                plt.show()
                
            # Do Chases
            cureve = curevents[ia]
            idx2 = np.where((events == cureve) & (np.abs(phi)<wideLimit))[0]
            usePhi = phi[idx2]
            useSes = corr[idx2]/norm[idx2]
            useSes_r = corr_r[idx2]/norm_r[idx2]
            xxx = np.arange(len(usePhi))
            oUseSes = np.copy(useSes)
            idx3 = np.where((np.abs(usePhi) > phiDur) | ((np.abs(usePhi)<phiDur) & (useSes<0.0)))[0]
            if len(idx3) > 0:
                usePhi = usePhi[idx3]
                useSes = useSes[idx3]
                useSes_r = useSes_r[idx3]
                closePhi = np.min(np.abs(usePhi))
                # Set threshold based upon reference sig
                threshOld = chasesPeakFac*maxSes
                idxtmp = np.where(useSes_r<0.0)[0]
                if len(idxtmp) > 0:
                    threshNew = np.max(np.abs(useSes_r[idxtmp])*1.1)
                else:
                    threshNew = 0.0
                useThresh = np.max([threshOld, threshNew])
                # For low MES things just use original threshold
                if (origMes<10.0):
                    useThresh = threshOld
                chsThrSumOld = chsThrSumOld + threshOld
                chsThrSumNew = chsThrSumNew + threshNew
#                idx4 = np.where(np.abs(useSes) > chasesPeakFac*maxSes)[0]
                idx4 = np.where(np.abs(useSes) > useThresh)[0]
                chase_val = 1.0
                if len(idx4) > 0:
                    nexPhi = np.min(np.abs(usePhi[idx4]))
                    chase_val = (nexPhi-closePhi)/wideLimit
                if debug:
                    print('Chases Thresh Ref: {0:f} Old: {1:f}'.format(threshNew, threshOld))
                    print('Individual Chases Value: {0:f}'.format(chase_val))

                    plt.plot(xxx, oUseSes, '.')
                    plt.plot(xxx[idx3], useSes, '.')
                    plt.plot(xxx[idx3], useSes_r, '.')
                    plt.show()
            else:
                chase_val = 0.0
            all_chases = np.append(all_chases, chase_val)
        nChs = len(all_chases)
        print('Chases Threshold New: {0:f} Old: {1:f}'.format(chsThrSumNew/nChs, chsThrSumOld/nChs))    
    else:
        # only a single transit event in data set this data bad
        all_ses = np.append(all_ses, 0.0)
        all_chases = np.append(all_chases, 0.0)
        all_corr = np.append(all_corr, 0.0)
        all_norm = np.append(all_norm, 1.0)
        all_time = np.append(all_time, 0.0)
        all_cadNo = np.append(all_cadNo, 0)
        
    # Recalculate the mes
    sumNumer = np.sum(all_corr)
    sumDenum = np.sqrt(np.sum(all_norm*all_norm))
    newMes = sumNumer/sumDenum
    sumNumer_r = np.sum(all_corr_r)
    sumDenum_r = np.sqrt(np.sum(all_norm_r*all_norm_r))
    newMes_r = sumNumer_r/sumDenum_r
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
        # same calc for ref sig
        ses2Mes_r = np.max(all_ses_r) / newMes_r
        mnNorm_r = np.mean(norm_r)
        depthEst_r = all_ses_r/all_norm_r
        all_mnses_r = depthEst_r * mnNorm_r
        all_mncorr_r = all_mnses_r * mnNorm_r
        sumNumer_r = np.sum(all_mncorr_r)
        sumDenum_r = np.sqrt(mnNorm_r*mnNorm_r*len(all_mncorr_r))
        newMnMes_r = sumNumer_r/sumDenum_r
        mnSes2mnMes_r = np.max(all_mnses_r) / newMnMes_r
    else:
        ses2Mes = 3.0
        newMnMes = newMes
        mnSes2mnMes = 3.0
        ses2Mes_r = 3.0
        newMnMes_r = newMes_r
        mnSes2mnMes_r = 3.0
    # Calculate chases
    if newNtran>2 & newNtran<5:
        chases_sum = np.median(all_chases)
    elif newNtran>=5:
        chases_sum = np.mean(all_chases)
    elif newNtran==2:
        chases_sum = np.min(all_chases)
    else:
        chases_sum = 0.0

    print('new ', newMes, newMes_r, ses2Mes, ses2Mes_r, newNtran, chases_sum, newMes/np.sqrt(newNtran), newMnMes, newMnMes_r, mnSes2mnMes, mnSes2mnMes_r)

    return newMes, newMes_r, ses2Mes, ses2Mes_r, newNtran, chases_sum, all_ses, \
        all_chases, newMnMes, newMnMes_r, mnSes2mnMes, mnSes2mnMes_r, all_corr, \
        all_norm, all_time, all_cadNo
    
if __name__ == "__main__":
    # Parse the command line arguments for multiprocessing
    # With Gnu parallel with 13 cores
    # seq 0 12 | parallel --results ses_mes_results python ses_mes_stats.py -w {} -n 13
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", type=int,\
                        default = 0, \
                        help="Worker ID Number 0 through nWrk-1")
    parser.add_argument("-n", type=int,\
                        default = 1, \
                        help="Number of Workers")
    

    args = parser.parse_args() 
    # These are for parallel procoessing
    wID = int(args.w)
    nWrk = int(args.n)
    # Load the h5 file that contains TCE seed information
    # The h5 file is created by gather_tce_fromdvxml.py
    tceSeedInFile = 'sector32_20200125_tce.h5'
    #  Directory storing the resampled dv time series data
    dvDataDir = '/pdo/users/cjburke/spocvet/sector32'
    # Directory of output hd5 files
    outputDir = dvDataDir
    SECTOR = 32
    # What fraction of data can be missing and still calculat ses_mes
    # In Sector 1 due to the 2 days of missing stuff it was 0.68
    validFrac = 0.52
    overWrite = False

    # Skyline data excises loud cadecnes
    skyline_file = 'skyline_data_sector32_20200125.txt'
    if os.path.isfile(skyline_file):
        dataBlock = np.genfromtxt('skyline_data_sector32_20200125.txt', dtype=['f8'])
        badTimes = dataBlock['f0']
        if len(badTimes) < 2:
            badTimes = np.array([0.0])
    else:
        print('No skyline data found')
        badTimes = np.array([0.0])

    # Search and filter parameters
    cadPerHr = 6
    firstFilterScaleFac = 10 # low frequency median filter will be
                            # firstFilterScaleFac*searchDurationHours medfilt window

    # Load the tce data h5
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)

    # These next few lines can be used to examine a single target    
    #all_epics = np.array([x.epicId for x in all_tces], dtype=np.int64)
    #all_pns = np.array([x.planetNum for x in all_tces], dtype=np.int)
    #ia = np.where((all_epics == 24843062) & (all_pns == 1))[0]
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
            depth = 0.0
            if td.at_valid == 1:
                period = td.at_period
                epoch = td.at_epochbtjd
                duration = td.at_dur
                depth = td.at_depth
            elif td.trp_valid == 1:
                period = td.tce_period
                epoch = td.trp_epochbtjd
                duration = td.trp_dur
                depth = td.trp_depth
            else:
                period = td.tce_period
                epoch = td.tce_epoch
                duration = td.pulsedur
                depth = 1000.0
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
                cadNo = np.array(f['cadenceNo'])
                for jj, curTime in enumerate(time):
                    diff = np.abs(curTime - badTimes)
                    if np.min(diff) < 5.0/60.0/24.0:
                        vd[jj] = False
                if doDebug:
                    tmpx = np.arange(len(useFlux))
                    plt.plot(useFlux, '.')
                    plt.plot(tmpx[vd], useFlux[vd], '.')
                    
                    plt.show()
                    
                # assign sector numbers to data
                secnum = np.ones_like(cadNo, dtype=np.int)
                idxSec = np.where(td.all_sectors>=0)[0]
                td.all_sectors = td.all_sectors[idxSec]
                td.all_cadstart = td.all_cadstart[idxSec]
                td.all_cadend = td.all_cadend[idxSec]
                nSec = len(td.all_sectors)
                if len(td.all_sectors)>1:
                    for kk in range(nSec):
                        idx = np.where((cadNo >= td.all_cadstart[kk]) & (cadNo <= td.all_cadend[kk])  )[0]
                        secnum[idx] = td.all_sectors[kk]
                
                # Mark data in transit as deweighted during detrending
                ootvd = np.full_like(vd, True)
                phi = phaseData(time, period, epoch)
                idx = np.where(np.abs(phi) < duration/2.0/24.0/period)[0]
                ootvd[idx] = False
                
                # Make a reference transit signal shape
                ztmp = phi * period
                ref_sig = trapezoid(ztmp, depth/1.0e6, duration/24.0, duration/24.0/1.1)
#                if doDebug:
#                    plt.plot(ref_sig, '.')
#                    plt.show()
                # Determine empirical noise from flux time series
                vdootIdx = np.where(vd & ootvd)[0]
                flux_level = np.median(useFlux[vdootIdx])
                flux_diff = np.diff(useFlux[vdootIdx]) / flux_level
                emp_noise = robust.mad(flux_diff)/np.sqrt(2.0)
                print('Empirical Noise [ppm]: {0:.1f}'.format(emp_noise*1.0e6))
#                if doDebug:
#                    plt.plot(flux_diff, '.')
#                    plt.show()
                # Add noise to reference transit signal
                ref_noise = np.random.randn(len(useFlux))*emp_noise + 1.0
                ref_sig = ref_sig*ref_noise
#                if doDebug:
#                    plt.plot(ref_sig, '.')
#                    plt.show()
                
                # For ETE6 only there was unexpectedly extra data after 
                #  indice 1348. Mark those aftr 1348 as invalid now
                #vd[1349:] = False
                
                # For Sector 1 data there was a crappy part 
                #vd[2465:2750] = False
                # detrend data fixe edges
                # Verify that there is enough data to do analysis
                # Also do a period cut
                
                # Need to update the valid data fraction for multi-sector
                # For single sector this would catch data that was missing
                # because previous planets detected removed too much data
                #    avoid the 'swiss' cheese light curves
                #  but for multi-sector the light curve can have long gaps
                #  of missing sectors this is not the intent of the check
                # Need to look for long gaps and remove these from the calculation
                # Protect against no valid data
                idxGd = np.where(vd)[0]
                if not len(idxGd) == 0:
                    # Time differences
                    gdTime = time[idxGd]
                    difftime = np.diff(gdTime)
                    jumpidx = np.where(np.abs(difftime)>7.0)[0]
                    vdFracUse = np.copy(vd)
                    timeFracUse = np.copy(time)
                    for jj in jumpidx:
                        jumptimebeg = gdTime[jj]
                        jumptimeend = gdTime[jj+1]
                        itmp = np.where((timeFracUse >=jumptimebeg) & (timeFracUse<=jumptimeend))[0]
                        ibeg = np.min(itmp)
                        iend = np.max(itmp)
                        vdFracUse = np.delete(vdFracUse, np.arange(ibeg,iend))
                        timeFracUse = np.delete(timeFracUse, np.arange(ibeg,iend))
                    # Also need to remove from last valid time
                    tmpx = np.arange(len(vdFracUse))
                    maxGdidx = np.max(tmpx[vdFracUse])
                    itmp = np.where(timeFracUse>=timeFracUse[maxGdidx])[0]
                    vdFracUse = vdFracUse[0:itmp[0]]
                    if len(vdFracUse)>0:
                        vdfrac = float(len(np.where(vdFracUse)[0]))/len(vdFracUse)
                    else:
                        vdfrac = 0.0
                else:
                    vdfrac = 0.0
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
                    #doDebug=False
                    final_smooth_flux, bad_edge_flag = flux_cond.detrend_with_smoothn_edgefix(\
                                        useFlux, vd, ootvd, int(np.ceil(cadPerHr*searchDurationHours)), fixEdge=True, \
                                         medfiltScaleFac=10, gapThreshold=5, edgeExamWindow=8, \
                                         edgeSig=6.0, edgeMinCad=50, debug=doDebug, secNum=secnum)
                    #doDebug=True
                    # fill gaps and extend to power of two
                    fillWindow = int(np.ceil(cadPerHr*searchDurationHours*firstFilterScaleFac*3))
                    final_smooth_flux_ext, vd_ext = flux_cond.fill_extend_fluxts(final_smooth_flux, \
                                                                vd, fillWindow, doExtend=True, debug=False)
                    final_ref_sig_ext, vd_ext = flux_cond.fill_extend_fluxts(ref_sig, \
                                                                vd, fillWindow, doExtend=True, debug=False)
                    #plt.cla()
                    #plt.plot(final_smooth_flux_ext, '.')
                    #plt.plot(final_ref_sig_ext, '.')
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

                    wavObj_ref = kw.waveletObject(waveletLen, final_ref_sig_ext, varianceFilterWindow)
                    trial_pulse = kw.set_trial_transit_pulse(searchLen)
                    # Replace the whitening coefficients from the original time series
                    wavObj_ref.whiteningCoeffs = wavObj.whiteningCoeffs
                    normTS_ref, corrTS_ref = kw.compute_statistic_time_series(wavObj_ref, searchLen, trial_pulse)
                    #plt.cla()
                    #plt.plot(corrTS/normTS, '.')
                    #plt.plot(corrTS_ref/normTS_ref, '.')
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
                    useNormTS_ref = normTS_ref[vd_ext]
                    useCorrTS_ref = corrTS_ref[vd_ext]
                    usePhase = phaseData(time[vd], period, epoch)
                    phaseDur = duration / 24.0 / period
                    useEvents = assignEvents(time[vd], epoch, usePhase, period, phaseDur)
                    useTime = time[vd]
                    useCadNo = cadNo[vd]
                    #doDebug=False
                    newMes, newMes_r, ses2Mes, ses2Mes_r, newNTran, Chases_sumry, \
                    allSes, allChases, \
                            newMnMes,newMnMes_r, \
                            mnSes2mnMes, mnSes2mnMes_r, \
                            allCorr, allNorm, allTime, allCadNo = get_ses_stats(useCorrTS, useNormTS,
                                                useCorrTS_ref, useNormTS_ref, \
                                                usePhase, phaseDur, useEvents, useTime, useCadNo, origMes, debug=doDebug)
                    #doDebug = True
                    if newNTran > 1:
                        validSes = 1
                    else:
                        validSes = 0
                    validAltDet = 1
                else:
                    print('Too little data for analysis epic: {0:d} pn: {1:d}'.format(epicid, pn))
                    print('NDat: {0:d} P: {1:f} phaseCoverage: {2:f} validFraction {3:f}'.format(len(np.where(vd)[0]),period, duration/24.0/period,vdfrac))
                    validAltDet = 0
                    validSes = 0
                    newMes = 0.0
                    ses2Mes = 3.0
                    newNTran = 1
                    Chases_sumry = 0.0
                    allSes = np.array([0.0])
                    allChases = np.array([0.0])
                    allCorr = np.array([0.0])
                    allNorm = np.array([0.0])
                    allTime = np.array([0.0])
                    allCadNo = np.array([0], dtype=np.int)
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
                    newMes_r = 0.0
                    ses2Mes_r = 3.0
                    newMnMes_r = 0.0
                    mnSes2mnMes_r = 3.0
                    
                #print(fileOutput)
                f = h5py.File(fileOutput,'w')
                tmp = f.create_dataset('altDetrend', data=final_smooth_flux, compression='gzip')
                tmp = f.create_dataset('validData', data=vd, compression='gzip')
                tmp = f.create_dataset('time', data=time, compression='gzip')
                tmp = f.create_dataset('cadNo', data=cadNo, compression='gzip')
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
                tmp = f.create_dataset('allSes', data=allSes, compression='gzip')
                tmp = f.create_dataset('allChases', data=allChases, compression='gzip')
                tmp = f.create_dataset('allCorr', data=allCorr, compression='gzip')
                tmp = f.create_dataset('allNorm', data=allNorm, compression='gzip')
                tmp = f.create_dataset('allTime', data=allTime, compression='gzip')
                tmp = f.create_dataset('allCadNo', data=allCadNo, compression='gzip')
                
                tmp = f.create_dataset('newMnMes', data=np.array([newMnMes], dtype=np.float))
                tmp = f.create_dataset('mnSes2mnMes', data=np.array([mnSes2mnMes], dtype=np.float))
                tmp = f.create_dataset('phaseDur', data=np.array([phaseDur], dtype=np.float))
                tmp = f.create_dataset('newMes_r', data=np.array([newMes_r], dtype=np.float))
                tmp = f.create_dataset('ses2Mes_r', data=np.array([ses2Mes_r], dtype=np.float))
                tmp = f.create_dataset('newMnMes_r', data=np.array([newMnMes_r], dtype=np.float))
                tmp = f.create_dataset('mnSes2mnMes_r', data=np.array([mnSes2mnMes_r], dtype=np.float))

                f.close()
                