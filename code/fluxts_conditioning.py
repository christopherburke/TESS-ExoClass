#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 12:14:57 2018
This is a flux time series conditioner that detrends (see smoothn.py),
fills gaps, fixes edges, and extends data to a power of 2
to prepare the time series for the Kepler wavelet machinary.
It also allows for protecting/deweighting in transit data
from the detrending.
There is so much going on in here just pray you
don't need to fiddle with the parameters.  

@author: Christopher J. Burke
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.stats as st
import scipy.signal as sig
import smoothn as smth
from statsmodels import robust
from spectrum import arburg
import kep_wavelets as kw
import math

def fill_noise(nW, preY, postY, SIG=3.0, orderAR=4, debug=False):
    # SIG - Before calculating fill remove SIG outliers
    # orderAR - order of the AR model to fit to data

    # For misbehaving light curves
    # too much data is clipped out such that <orderAR data is left
    #  and arburg errors out
    # detect when too much data has been removed and increase SIG and 
    #  dont use robust mad
    oLen = len(preY)
    preYStd = robust.mad(preY)
    preYMn = np.mean(preY)
    preY = preY - preYMn
    idx = np.where((np.abs(preY/preYStd)<SIG))[0]
    sigmults = [1.0, 1.0, 1.5, 2.0]
    if len(idx) < math.floor(0.8*oLen):
        # Too much data removed try to increase SIG and less robust std
        cnt = 0
        while cnt < 2 and len(idx) < math.floor(0.8*oLen):
            cnt = cnt + 1
            preYStd = np.std(preY)
            newSig = SIG*sigmults[cnt]
            idx = np.where((np.abs(preY/preYStd)<newSig))[0]
    preY = preY[idx]
    oLen = len(postY)
    postYStd = robust.mad(postY)
    postYMn = np.mean(postY)
    postY = postY - postYMn
    idx = np.where((np.abs(postY/postYStd)<SIG))[0]
    if len(idx) < math.floor(0.8*oLen):
        # Too much data removed try to increase SIG and less robust std
        cnt = 0
        while cnt < 2 and len(idx) < math.floor(0.8*oLen):
            cnt = cnt + 1
            postYStd = np.std(postY)
            newSig = SIG*sigmults[cnt]
            idx = np.where((np.abs(postY/postYStd)<newSig))[0]
    postY = postY[idx]
    if len(preY) >= orderAR:
        # Get the autocorrelation coeffs from preData
        preAR, preVar, k = arburg(preY, orderAR)
        # Now make the preFillData
        # Construt initial conditions
        lfil_ic = sig.lfiltic([1], np.insert(preAR,0,1.0), preY[-orderAR:])
        lfilter_result = sig.lfilter([1], np.insert(preAR,0,1.0), \
                np.random.normal(scale=np.sqrt(preVar), size=(nW,)), \
                zi=lfil_ic)
        preFill = np.real(lfilter_result[0])
    else:
        # as backup if the length of the pre data is too short do random
        preFill = np.random.normal(scale=preYStd, size=(nW,))
        preVar = np.power(preYStd, 2)
    if debug:
        print(np.var(preFill), preVar)

    # Now do postData
    if len(postY) >= orderAR:
        postAR, postVar, k = arburg(postY, orderAR)
        lfil_ic = sig.lfiltic([1], np.insert(postAR,0,1.0), np.flip(postY[0:orderAR], axis=0))
        lfilter_result = sig.lfilter([1], np.insert(postAR,0,1.0), \
                np.random.normal(scale=np.sqrt(postVar), size=(nW,)), \
                zi=lfil_ic)
        postFill = np.real(lfilter_result[0])
    else:
        postFill = np.random.normal(scale=postYStd, size=(nW,))
        postVar = np.power(postYStd, 2)
    # final fill is linear weighted combination
    runN = np.arange(nW)
    preRamp = runN/np.float(nW)
    postRamp = 1.0 - preRamp
    # This scaling factor is necessary to linearly transition from the variance
    #  of the pre portion to the post portion  
    scl = np.sqrt((preRamp*(postVar-preVar)+preVar)/(preRamp*preRamp*preVar+np.power(1.0-preRamp,2)*postVar))

    finalFill = (preFill*preRamp + postFill*postRamp)*scl
    if debug:
        print(np.var(postFill), postVar)

        allN = len(preY) + len(finalFill) + len(postY)
        runN = np.arange(allN)
        showFill = preY
        showClr = np.zeros_like(preY)
        showFill = np.append(showFill, finalFill)
        showClr = np.append(showClr, np.ones_like(finalFill))
        showFill = np.append(showFill, postY)
        showClr = np.append(showClr, np.zeros_like(postY))
        plt.scatter(runN, showFill, c=showClr)
        plt.show()
    return finalFill
    
def detrend_with_smoothn_edgefix(flux, vd, ootvd, durat, fixEdge=True, \
                                 medfiltScaleFac=10, gapThreshold=5, edgeExamWindow=8, \
                                 edgeSig=6.0, edgeMinCad=50, debug=False, secNum=None):
    # Detrend light curve in pieces using the smoothn routine
    # vd is valid data True==valid ; False==bad
    # ootvd is out of transit data==True ; in transit data(protected from detrending) == False
    # durat - [cadences] transit duration to protect during smoothing
    # fixEdge - identify and fix edge
    # gapThreshold - [cadences] gaps of this size will be fit as separate pieces
    # medfiltScaleFac - the initial medfilt window is durat*medfiltScaleFac
    #       First fit linear trend across sector
    # edgeExamWindow - [cadences] Cadences at beg and end of data block are examine
        # for strong positive outliers indicating 'wild' data at edges
    # edgeSig - significance of the outliers in the edge to trigger correction
    # edgeMinCad - how many cadences in data block required to attempt bad edge detection
    # secNum if specified set each sector to have the same mean flux as all data combined before
    #   starting detrending
    filterCircularShift = 20 
    bad_edge_flag = np.zeros_like(flux)
    runcad = np.arange(len(flux))
    useflux = np.copy(flux)
    final_smooth_flux = np.zeros_like(useflux)
    # combvd is the valid data and oot data combined
    combvd = np.logical_and(vd, ootvd)
    tmpcad = runcad[combvd]
    tmpflux = useflux[combvd]
    tmpbadedge = bad_edge_flag[combvd]
    tmpcadwt = runcad[vd]
    tmpfluxwt = useflux[vd]
    tmpbadedgewt = bad_edge_flag[vd]
    tmpootvd = ootvd[vd]
    
    # pre step to make all sectors have same flux level
    if not secNum is None:
        if debug:
            plt.plot(tmpcadwt, tmpfluxwt, '.')
        fluxmedian = np.median(tmpflux)
        tmpsecnum = secNum[combvd]
        tmpsecnumwt = secNum[vd]
        uniqsec = np.unique(tmpsecnum)
        # Make sure all cadences have a sector number assigned in the valid data
        idxany = np.where(tmpsecnum == 0)[0]
        if len(idxany)>0:
            print('Valid data with no sector number!')
            exit()
        for cursec in uniqsec:
            idxsec = np.where(tmpsecnum == cursec)[0]
            tmpflxmedian = np.median(tmpflux[idxsec])
            tmpflux[idxsec] = tmpflux[idxsec] / tmpflxmedian * fluxmedian
            idxsecwt = np.where(tmpsecnumwt == cursec)[0]
            tmpfluxwt[idxsecwt] = tmpfluxwt[idxsecwt] / tmpflxmedian * fluxmedian
        if debug:
            plt.plot(tmpcadwt, tmpfluxwt, '.')
            plt.show()
    # First step is to linear fit to all data to make it somewhat
    #  more same and normalized
    linfit = st.linregress(tmpcad, tmpflux)
    tmp_linflat = tmpflux / (linfit[0]*tmpcad + linfit[1])
    tmp_linflat_wt = tmpfluxwt / (linfit[0]*tmpcadwt + linfit[1])
    if debug:
        plt.plot(tmpcadwt, tmpfluxwt, '.')
        plt.plot(tmpcadwt, linfit[0]*tmpcadwt + linfit[1], '-')
        plt.show()
        plt.plot(tmpcadwt, tmp_linflat_wt, '.')
        plt.show()
    
    # Step 2 make a median filtered lightcurve that only keeps very long
    #  time scale variability
    medfilterlen = np.int(durat * medfiltScaleFac)
    # Make sure medfilterlen is odd
    if np.mod(medfilterlen,2) == 0:
        medfilterlen = medfilterlen+1
    tmp_medfilt = sig.medfilt2d(tmp_linflat.reshape(1,-1), (1, medfilterlen))[0]
    # medfilt2d is much faster than medfilt
    #tmp_medfilt = sig.medfilt(tmp_linflat, medfilterlen)
    # We only need this median filter to get the alpha parameter from smoothn
    # Trim away the begining and end to remove edge effects
    # Protect agains tmp_medfilt being shorter than twice medfilterlen
    if len(tmp_medfilt) > 3*medfilterlen:
        tmp2_medfilt = tmp_medfilt[medfilterlen:-medfilterlen]
    else:
        medfilterlen = int(0.1*len(tmp_medfilt))
        tmp2_medfilt = tmp_medfilt[medfilterlen:-medfilterlen]
    smthresult = smth.smoothn(tmp2_medfilt, np.ones_like(tmp2_medfilt))
    medsmthparm = smthresult[2]
    # Scale the smoothing parameter by the variance in the median filtered and
    #  the native time series.  This results in smoothing of similar scale
    #  when smoothn is used with the native time scale.  Not too sure why the roll
    #   function is used, but doing it because this was done in matlab version
    fluxAmp = (np.max(smthresult[0]) - np.min(smthresult[0])) / np.std(tmp_linflat[medfilterlen:-medfilterlen]-smthresult[0])
    smthparm  = np.var(np.diff(np.roll(tmp_linflat, filterCircularShift) ) ) \
                    / np.var( np.diff(tmp2_medfilt-smthresult[0]) )  * medsmthparm * fluxAmp
    print(medsmthparm, smthparm, smthparm/medsmthparm, fluxAmp)
    if not np.isfinite(smthparm):
        print("Bad smoothing parameter!")
#    else:
#        print("Smoothing parameter: {0:f}".format(smthparm))
    if debug:
        plt.plot(tmpcad, tmp_linflat, '.')
        plt.plot(tmpcad, tmp_medfilt, '-')
        plt.show()
        plt.plot(tmp2_medfilt, '.')
        plt.plot(smthresult[0], '-')
        plt.show()
    # Look for gaps bigger than gapThreshold
    idxGaps = np.where(np.diff(tmpcadwt) > gapThreshold)[0]
    idxGapStart = np.array([0])
    idxGapEnd = np.array([len(flux)])
    if len(idxGaps)>0:
        idxGapStart = np.append(idxGapStart, idxGaps+1)
        idxGapEnd = np.append(idxGaps+1, len(flux))
    tmp_smooth_flux = np.zeros_like(tmp_linflat_wt)   
    for j in range(len(idxGaps)+1):
        ist = idxGapStart[j]
        ien = idxGapEnd[j]
        curflux = tmp_linflat_wt[ist:ien]
        wght = np.ones_like(curflux)
        curootvd = tmpootvd[ist:ien]
        wght = np.where(np.logical_not(curootvd), 0.0, wght)
        # Set flux to nan for in transit data so that the smoothn
        #  function will know to linearly interpolate over these 
        # cadences for its initial guess and ignore them for the
        # smoothing solution.  At least that is what is supposed
        #  to happen in principle.
        ucurflux = np.where(np.logical_not(curootvd), np.nan, curflux)
        # Check to make sure there is at least three valid data points in current chunk
        if len(np.where(np.isfinite(ucurflux))[0])>3:
            smthresult = smth.smoothn(ucurflux, wght, smthparm)
            cursmoothflux = smthresult[0]
        else:
            cursmoothflux = np.ones_like(curflux)
        # Check for bad smoothing
        if not np.all(np.isfinite(cursmoothflux)):
            print("non finite smoothing detected")
            print("{:d} {:d}".format(len(curflux), len(np.where(curootvd)[0])))
            cursmoothflux = np.ones_like(curflux)
        if debug:
            plt.plot(curflux, '.')
            plt.plot(cursmoothflux, '-')
            plt.show()
        tmp_smooth_flux[ist:ien] = curflux / cursmoothflux
        if debug:
            plt.plot(tmp_smooth_flux[ist:ien], '.')
            plt.show()
        # Check for edge issues
        if len(curflux) > edgeMinCad and fixEdge:
            # Calculate noise in data before edge at end
            tmpsmth = tmp_smooth_flux[ist:ien]
            preEdge = tmpsmth[len(curflux)-4*edgeExamWindow:len(curflux)-edgeExamWindow]
            preEdgeMad = robust.mad(preEdge)
            edge = tmpsmth[len(curflux)-edgeExamWindow-1:]
            edgePosOutlier  = np.max((edge-1.0)/preEdgeMad)
            if edgePosOutlier > edgeSig:
                #print('Bad end edge detected')
                fullEdge = np.copy(tmpsmth[len(curflux)-4*edgeExamWindow:])
                subsmthresult = smth.smoothn(fullEdge, np.ones_like(fullEdge)*preEdgeMad, smthparm/100.0, robust=False)
                subsmoothflux = subsmthresult[0]
                tmp_smooth_flux[ist:ien][len(curflux)-4*edgeExamWindow:] = (fullEdge/subsmoothflux)
                tmpbadedgewt[ist:ien][len(curflux)-edgeExamWindow-1:] = 1.0
                if debug:
                    plt.plot(fullEdge, '.')
                    plt.plot(subsmoothflux,'-')
                    plt.show()
                    plt.plot(tmp_smooth_flux[ist:ien], '.')
                    plt.show()
            # Calculat noise in data after edge at beginning
            tmpsmth = tmp_smooth_flux[ist:ien]
            postEdge = tmpsmth[edgeExamWindow+1:4*edgeExamWindow]
            postEdgeMad = robust.mad(postEdge)
            edge = tmpsmth[0:edgeExamWindow+1]
            edgePosOutlier = np.max((edge-1.0)/postEdgeMad)
            if edgePosOutlier > edgeSig:
                #print('Band starting edge detected')
                fullEdge = np.copy(tmpsmth[0:4*edgeExamWindow])
                subsmthresult = smth.smoothn(fullEdge, np.ones_like(fullEdge)*postEdgeMad, smthparm/100.0, robust=False)
                subsmoothflux = subsmthresult[0]
                tmp_smooth_flux[ist:ien][0:4*edgeExamWindow] = (fullEdge/subsmoothflux)
                tmpbadedgewt[ist:ien][0:edgeExamWindow+1] = 1.0
                if debug:
                    plt.plot(fullEdge, '.')
                    plt.plot(subsmoothflux)
                    plt.show()
                    plt.plot(tmp_smooth_flux[ist:ien], '.')
                    plt.show()
                
    final_smooth_flux[vd] = tmp_smooth_flux
    bad_edge_flag[vd] = tmpbadedgewt
    if debug:
        print("done detrending")
        plt.plot(runcad[vd], final_smooth_flux[vd], '.')
        plt.show()
    return final_smooth_flux, bad_edge_flag

def fill_extend_fluxts(flux, vd, noiseWindow, doExtend=True, debug=False):
    #  fill all gaps and extend to the next power of two
    # Work on extending flux, detrending, filling gaps and identify outliers
    nBand = np.int(np.ceil(np.log2(len(flux))))
    wantNExtend = np.int(np.power(2,nBand))
    origN = len(flux)

    # Deal with case where first cadence is invalid
    frstCadOK = True
    if not vd[0]:
        frstCadOK = False
        # temporarily make the first cadence valid
        vd[0] = True
        print('First cadence invalid')
   # Next is extend to the power of two
    valid_data_flag_ext = np.zeros((wantNExtend,), dtype=np.bool_)
    valid_data_flag_ext[0:origN] = np.copy(vd)
    final_smooth_flux_ext = np.zeros((wantNExtend,))
    final_smooth_flux_ext[0:origN] = np.copy(flux)

    runcad_ext = np.arange(wantNExtend)
    vd_ext = np.where(valid_data_flag_ext)[0]
    tmpcad_ext = runcad_ext[vd_ext]
    
    # Look for even single gaps
    idxGaps = np.where(np.diff(tmpcad_ext) > 1)[0]
    #plt.plot(flux, '.')
    #plt.show()
     # want the index from the original length
    idxGapStart = tmpcad_ext[idxGaps]
    # Add in the last valid cadence
    idxGapStart = np.append(idxGapStart, tmpcad_ext[-1])+1
    idxGapEnd = tmpcad_ext[idxGaps+1]
    idxGapEnd = np.append(idxGapEnd, wantNExtend)
    # Now that we have determined gaps unfix the situation where the first cadence was invalid
    if not frstCadOK:
        # if there are additional false after first this is taken care of
        if vd[1] == False:
            vd[0] = False
            valid_data_flag_ext[0] = False
            idxGapStart[0] = 0
        else: # first and only first was bad
            vd[0] = False
            valid_data_flag_ext[0] = False
            idxGapStart = np.insert(idxGapStart, 0, 0)
            idxGapEnd = np.insert(idxGapEnd, 0, 1)
    # Do the fills backward order in order to
    # First do the extension fill
    # Allow previous filled data to now be valid
    # Determine the median gap size
    medGapSz = np.int(np.median(idxGapEnd-idxGapStart))
    # set noiseWindow to bigger of noiseWindow or twice medGapSz
    oNoiseWindow = noiseWindow
    noiseWindow = np.max([noiseWindow, 2*medGapSz])
    #print("Fill Window Used: {0:d} request: {1:d}".format(noiseWindow, oNoiseWindow))
    tmp_valid_data_flag_ext = np.copy(valid_data_flag_ext)
    for j in reversed(range(len(idxGapStart))):
        nFill = idxGapEnd[j] - idxGapStart[j]
        js = idxGapStart[j]
        je = idxGapEnd[j]
        noiseExtWindow = noiseWindow
        prejs = js-noiseExtWindow
        preje = js
        if (prejs < 0):
            prevalid = tmp_valid_data_flag_ext[prejs:]
            prevalid = np.append(prevalid, tmp_valid_data_flag_ext[0:preje])
            preData = final_smooth_flux_ext[prejs:]
            preData = np.append(preData, final_smooth_flux_ext[0:preje])
        else:
            prevalid = tmp_valid_data_flag_ext[prejs:preje]
            preData = final_smooth_flux_ext[prejs:preje]

        postjs = je
        postje = je+noiseExtWindow
        if (postje > wantNExtend):
            postvalid = tmp_valid_data_flag_ext[postjs:]
            postvalid = np.append(postvalid, tmp_valid_data_flag_ext[0:np.mod(postje,wantNExtend)])
            postData = final_smooth_flux_ext[postjs:]
            postData = np.append(postData, final_smooth_flux_ext[0:np.mod(postje, wantNExtend)])
        else:
            postvalid = tmp_valid_data_flag_ext[postjs:postje]
            postData = final_smooth_flux_ext[postjs:postje]
            
        if debug:
            print("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}".format(js, je, prejs, preje, \
                  postjs, postje))
        idx = np.where(prevalid)[0]
        preData = preData[idx]
        idx = np.where(postvalid)[0]
        postData = postData[idx]

        fillData = fill_noise(nFill, preData, postData, debug=debug)
        final_smooth_flux_ext[js:je] = fillData + 1.0
        tmp_valid_data_flag_ext[js:je] = True

    if debug:
        plt.plot(final_smooth_flux_ext, '.')
        plt.show()
    return final_smooth_flux_ext, valid_data_flag_ext

if __name__ == '__main__':
#    fileinput = 'tess_dvts_0000000229684206_01.h5d'
    fileinput = 'tess_dvts_0000000180349353_01.h5d'
   
    f = h5py.File(fileinput, 'r')
    cadenceNo = np.array(f['cadenceNo'])
    timetbjd = np.array(f['timetbjd'])
    lc_init = np.array(f['lc_init'])
    lc_init_err = np.array(f['lc_init_err'])
    lc_white = np.array(f['lc_white'])
    pdc_flux = np.array(f['pdc_flux'])
    pdc_flux_err = np.array(f['pdc_flux_err'])
    deweights = np.array(f['deweights'])
    valid_data_flag = np.array(f['valid_data_flag'])
    # For ETE6 only there was unexpectedly extra data after 
    #  indice 1348. Mark those aftr 1348 as invalid now
    valid_data_flag[1349:] = False
    
    vd = valid_data_flag
    # Search and filter parameters
    cadPerHr = 2
    searchDurationHours = 6.0
    firstFilterScaleFac = 10 # low frequency median filter will be
                            # firstFilterScaleFac*searchDurationHours medfilt window
    filterCircularShift = 20  
    gapThreshold = 5  # Gaps bigger than gapThreshold cadences are fit separately                      
    edgeExamWindow = 8 # 8 Cadences at end of data block are examine
        # for issues with the fit
    edgeSig = 6.0 # significance of the outliers in the edge to trigger correction
    edgeMinCad = 50 # how many cadences in data block required to attempt bad edge detection

    ootvd = np.full_like(vd, True)
    # detrend data fixe edges
    final_smooth_flux, bad_edge_flag = detrend_with_smoothn_edgefix(\
                                pdc_flux, vd, ootvd, int(np.ceil(cadPerHr*searchDurationHours)), fixEdge=True, \
                                 medfiltScaleFac=10, edgeExamWindow=8, \
                                 edgeSig=6.0, edgeMinCad=50, debug=True)
    # fill gaps and extend to power of two
    fillWindow = int(np.ceil(cadPerHr*searchDurationHours*firstFilterScaleFac*3))
    final_smooth_flux_ext, valid_data_flag_ext = fill_extend_fluxts(final_smooth_flux, \
                                                        vd, fillWindow, doExtend=True, debug=True)

    plt.plot(final_smooth_flux_ext, '.')
    plt.show()
    
    waveletLen = 12
    searchLen = np.int(np.round(cadPerHr * searchDurationHours))
    varianceFilterFactor = 15
    varianceFilterWindow = searchLen * varianceFilterFactor
    wavObj = kw.waveletObject(waveletLen, final_smooth_flux_ext, varianceFilterWindow)
    trial_pulse = kw.set_trial_transit_pulse(searchLen)
    
    normTS, corrTS = kw.compute_statistic_time_series(wavObj, searchLen, trial_pulse)

    plt.plot(corrTS/normTS, '.')
    plt.show()
    plt.plot(1.0e6/normTS, '.')
    plt.show()
    print("hello world")