#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a port of the Kepler wavelet search TPS machinary from matlab to python.
The Kepler pipeline is open source and available here
https://github.com/nasa/kepler-pipeline

The original matlab code that this code is a modification of was released under 
NASA Open Source Agreement v1.3 (NOSA v1.3)
NOSA require the following statement for modifications to the original software
to appear prominently 

Copyright  2017 United States Government as represented by the Administrator 
of the National Aeronautics and Space Administration.  All Rights Reserved.

NASA acknowledges the SETI Institute’s primary role in authoring and 
producing the Kepler Data Processing Pipeline under Cooperative Agreement
 Nos. NNA04CC63A, NNX07AD96A, NNX07AD98A, NNX11AI13A, NNX11AI14A, NNX13AD01A 
 & NNX13AD16A.

Portions of the Kepler Data Processing Pipeline software constitute modified
 Matlab scripts and may be utilized only in a manner consistent with the
 terms and conditions of the Mathworks Limited License, rather than the
 terms of this Agreement.  The Mathworks Limited License may be found 
 in the file MATHWORKS-LIMITED-LICENSE.docx.

Further NOSA license details for the original source code are available 
in the included document
kep_wavelets_py-NASA-OPEN-SOURCE-AGREEMENT.doc

kep_wavelets.py is the only code as part of the TESS-ExoClass project
that NOSA is applicable.

Modification from original matlab source to Python by
@author: Christopher J. Burke (MIT)
"""

import numpy as np
import matplotlib.pyplot as plt
from statsmodels import robust
import scipy.signal as sig
import scipy.io as sio

class waveletObject:
    h0 = np.array([])
    H = np.array([], dtype=np.cdouble)
    G = np.array([], dtype=np.cdouble)
    nBands = 0
    def __init__(self, waveletlen, fluxTS, varwindow):
        """ fluxTS must be a power of 2 """
        self.fluxTS = np.array(fluxTS).flatten()
        self.whiteningCoeffs = np.array([])
        self.waveletCoeffs = np.array([])
        # Build the wavelet and whitening coeffs   
        self.h0, tmp = daubcqf(waveletlen)
        self.nBands = calcNBands(waveletlen, len(fluxTS))
        self.H, self.G = self.set_filter_bank()
        self.waveletCoeffs = self.overcomplete_wavelet_transform()
        self.whiteningCoeffs = self.set_whitening_coefficients(varwindow)
        
        
    # Set the filter banks
    def set_filter_bank(self):
        wavObj = self
        nSamples = len(wavObj.fluxTS)
        filterLength = len(wavObj.h0)
        nBands = wavObj.nBands
        
        # construct the 4 basis vectors
        # matlab assumed h0 is 1,filterLength
        h0 = np.reshape(wavObj.h0, (filterLength,1))
        h1 = np.flipud(h0) * np.reshape(np.power(-1, np.arange(0,filterLength)), (filterLength,1))
        g0 = np.flipud(h0)
        g1 = np.flipud(h1)

        # construct the FFT of each of the vectors, with appropriate padding -- note that here we
        # explicitly show which are low-pass and which are high-pass
        HL = np.fft.fft(h0.flatten(), nSamples)
        HH = np.fft.fft(h1.flatten(), nSamples)
        GL = np.fft.fft(g0.flatten(), nSamples)
        GH = np.fft.fft(g1.flatten(), nSamples)
        #np.save('sfb_HLR', np.real(HL))
        #np.save('sfb_HLI', np.imag(HL))
        #np.save('sfb_HHR', np.real(HH))
        #np.save('sfb_HHI', np.imag(HH))
        #np.save('sfb_GLR', np.real(GL))
        #np.save('sfb_GLI', np.imag(GL))
        #np.save('sfb_GHR', np.real(GL))
        #np.save('sfb_GHI', np.imag(GL))
        # define the filters
        wavObj.G = np.zeros((nSamples, nBands), dtype=np.cdouble)
        wavObj.H = np.zeros((nSamples, nBands), dtype=np.cdouble)
        
        # define 2 vectors which will hold product of low-pass filters
        GLProduct = np.ones((nSamples,), dtype=np.cdouble)
        HLProduct = np.ones((nSamples,), dtype=np.cdouble)

        # Loop over bands
        for iBand in range(0,nBands):
            #on the last band, the GH and HH vectors have to be set to one, since the lowest band
            #sees only low-pass filters all the way down
            if iBand == nBands -1:
                HH = np.ones((nSamples,))
                GH = np.ones((nSamples,))
            wavObj.G[:,iBand] = GH * GLProduct
            wavObj.H[:,iBand] = HH * HLProduct
            # Increment the products of the low-pass filters
            GLProduct = GLProduct * GL
            HLProduct = HLProduct * HL
            # convert the elemental filters to the next band down in freq
            tmp = GL[0::2]
            GL = np.append(tmp, tmp)
            tmp = HL[0::2]
            HL = np.append(tmp, tmp)
            tmp = GH[0::2]
            GH = np.append(tmp, tmp)
            tmp = HH[0::2]
            HH = np.append(tmp, tmp)
            
#        print("hello world")
        #np.save('sfb_HR',np.real(wavObj.H))
        #np.save('sfb_HI',np.imag(wavObj.H))
        #np.save('sfb_GR',np.real(wavObj.G))
        #np.save('sfb_GI',np.imag(wavObj.G))
        return wavObj.H, wavObj.G

    def overcomplete_wavelet_transform(self, usets=None):
        wavObj = self
        nBands = wavObj.nBands
        filterLength = len(wavObj.h0)
        nSamples = len(wavObj.fluxTS)
        default = False
        if usets is None:
            default = True
            usets = wavObj.fluxTS
        if not len(usets) == nSamples:
            print("Warning the input time series to owt is not the same as was used to setup wavelet object!!!")
        waveletCoefficients = -1.0 * np.ones((nSamples, nBands))
        #% construct the FFT of the initial vector and repmat it to the # of bands
        Xoneband = np.reshape(np.fft.fft(usets, axis=0), (nSamples,1))
        #if default:
            #np.save('owt_usets', usets)
            #np.save('owt_XonebandR',np.real(Xoneband))
            #np.save('owt_XonebandI',np.imag(Xoneband))
        X = np.tile(Xoneband, (1,nBands))
        #% the wavelet expansion is ALMOST just the IFFT of X multiplied by H ...

        waveletCoefficients = np.real(np.fft.ifft(X * wavObj.H, axis=0))
        #if default:
        #    np.save('owt_wc1',waveletCoefficients)
        # Except for some circshifts
        for iBand in range(nBands):
            shiftIndex = np.min([iBand+1, nBands-1])
            nShift = filterLength*int(np.power(2, shiftIndex-1)) - int(np.power(2, shiftIndex-1))
            waveletCoefficients[:,iBand] = np.roll(waveletCoefficients[:,iBand], -nShift)

#        print("hello world")
        #if default:
        #    np.save('owt_wc2', waveletCoefficients)
        return waveletCoefficients

    def set_whitening_coefficients(self, varWindow, usewavc=None):
        wavObj = self
        nBands = wavObj.nBands
        nSamples = len(wavObj.fluxTS)
        
        if usewavc is None:
            usewavc = wavObj.waveletCoeffs
        whitec = np.zeros_like(usewavc)
        for iBand in range(nBands):
            if iBand == nBands-1:
                subtractMedianFlag = True
            else:
                subtractMedianFlag = False
            decimationFactor = int(np.power(2, iBand))
            whitec[:,iBand] = np.power(self.moving_circular_mad(usewavc[:,iBand], \
                              varWindow*decimationFactor, subtractMedianFlag), -2.0)
        #np.save('swc_whitec1', whitec)
        #% Look for bands that have excessively large whitening coefficients and
        #% set them to something reasonable if they do
        waveletSupportBuffer = 50
        outlierSigmaMultiplier = 6.0 #% 6.0 may need to become a module parameter and tuned
 
        #% an impulse has support of 2*2^iBand so multiply by buffer to be safe
        waveletSupportInCadences = (waveletSupportBuffer * 2* np.power(2, np.arange(1,nBands+1))).astype(int)
        suspectBandIndicator = waveletSupportInCadences >= nSamples

        meanWhiteningCoefficients = np.mean(whitec,axis=0) 
        overallMeanWhiteningCoefficients = np.median(meanWhiteningCoefficients[np.logical_not(suspectBandIndicator)].flatten())
        stdWhiteningCoefficients = robust.mad(meanWhiteningCoefficients[np.logical_not(suspectBandIndicator)])
        badBands = (meanWhiteningCoefficients-overallMeanWhiteningCoefficients) > outlierSigmaMultiplier*stdWhiteningCoefficients
        idxBad = np.where((badBands) & (suspectBandIndicator))[0]
        if len(idxBad)>0:
            for i in idxBad:
                whitec[:,i] = overallMeanWhiteningCoefficients
        #np.save('swc_whitec2', whitec)
        return whitec
    
    def moving_circular_mad(self, vec, window, subMedian=True):
        madValues = np.zeros_like(vec)
        # do a circular extension
        nSamples = len(vec)
        if window < nSamples-2:        
            # window should be odd for sig.medfilt
            if np.mod(window,2)==0:
                window = window+1
            vecCirc = np.insert(vec, 0, vec[nSamples-window:])
            vecCirc = np.append(vecCirc, vec[0:window])
            nSamplesCirc = len(vecCirc)
            ##%     if median subtracted is desired, compute the median; otherwise, set equal to zero
            if subMedian:
                medianValue = sig.medfilt(vecCirc, window)
            else:
                medianValue = np.zeros((nSamplesCirc,))
            tmp = np.abs(vecCirc-medianValue)
            # medfilt2d is much faster than medfilt based upon
            #https://gist.github.com/f0k/2f8402e4dfb6974bfcf1
            madValuesCirc = sig.medfilt2d(tmp.reshape(1,-1), (1, window))[0]
#            madValuesCirc = sig.medfilt(np.abs(vecCirc-medianValue), window)
            # How about convolve for moving mean
            #  This is much faster than medfilt2, but the signal is supprresed
            #  because mean is less robust to the signal, thus the CDPP
            # spikes around the signal
            #madValuesCirc = np.convolve(tmp, np.ones((window,))/window, mode='full')
            madValues = madValuesCirc[window:nSamplesCirc-window]
            madValues = madValues / 0.6745
        else:
            if subMedian:
                medianValue = np.median(vec)
            else:
                medianValue = 0.0
            madValues = madValues + np.median(np.abs(vec-medianValue))
            madValues = madValues / 0.6745
        return madValues

        
  
def calcNBands(wN, fN):
    return int(np.log2(fN) - int(np.floor(np.log2(wN))) + 1)

def daubcqf(N):
    """literal translation of daubcqf.m of Kepler pipeline"""

    if np.mod(N,2) == 1:
        print("No Daubechies filter exists for ODD length {0:d}".format(N))
    K = int(N/2)
    a = 1
    p = 1
    q = 1
    h0 = np.array([1.0, 1.0])
    
    for j in range(1,K):
        a = -a*0.25*(j+K-1)/j
        h0 = np.insert(h0, 0, 0.0) + np.append(h0, 0.0)
        negp = -p
        p = np.insert(negp, 0, 0) + np.append(p, 0)
        negp = -p
        p = np.insert(negp, 0, 0) + np.append(p, 0)
        zqz = np.insert(q, 0, 0.0)
        zqz = np.append(zqz, 0.0)
        q = zqz + a*p         

    q = np.sort(np.roots(q))
    qt = q[0:K-1]
    h0 = np.convolve(h0, np.real(np.poly(qt)))
    h0 = np.sqrt(2.0)*h0/np.sum(h0) # normalize to sqrt(2)
    if (np.abs(np.sum(np.power(h0,2))) - 1.0) > 1.0e-4:
        print("Numerically unstable Daubechies for this value of N {0:d}".format(N))
    h1 = np.rot90(np.reshape(h0, (len(h0),1)), 2).flatten()
    h1[0:2:N] = -h1[0:2:N]
    

    return h0, h1

def set_trial_transit_pulse(duration):
    """ duration is transit search width in cadences integer
         do a box signal"""
    trial_pulse = np.zeros((duration+1),)
    trial_pulse[0:duration] = -0.5
    trial_pulse[-1] = 0.0
    return trial_pulse

def compute_statistic_time_series(wavObj, searchLen, trial_pulse):
    whtC = wavObj.whiteningCoeffs    
    x = wavObj.waveletCoeffs
    nSamples = x.shape[0]
    nBands = x.shape[1]
    shiftLength =  int(np.fix(searchLen/2.0)) + 1
    #% zero pad the pulse so its the same length as x
    full_trial_pulse = np.zeros((nSamples,))
    full_trial_pulse[0:len(trial_pulse)] = trial_pulse
    s = wavObj.overcomplete_wavelet_transform(full_trial_pulse)
    #np.save('csts_s', s)
    corrTS = np.zeros((nSamples,))
    normTS = np.zeros((nSamples,))    
    for iBand in range(0,nBands-1):
        factorOfTwo = np.power(2.0, -np.min([iBand+1, nBands-1]))
        SNRi = circfilt(np.flip(s[:,iBand]*s[:,iBand], 0), whtC[:,iBand])
        Li = circfilt(np.flip(s[:,iBand], 0), x[:,iBand]*whtC[:,iBand])
        SNRi = np.roll(SNRi, shiftLength)
        Li = np.roll(Li, shiftLength)
        normTS = normTS + SNRi*factorOfTwo
        corrTS = corrTS + Li*factorOfTwo
        if iBand == nBands-2:
            normTS = np.sqrt(normTS)
    #np.save('csts_normTS', normTS)
    #np.save('csts_corrTS', corrTS)        
    return normTS, corrTS

def circfilt(vec1, vec2):
    nLength = len(vec2)
    X = np.fft.fft(vec2)
    H = np.fft.fft(vec1, nLength)
    y = np.real(np.fft.ifft(H*X))
    return y
    
if __name__ == '__main__':
    waveletLen = 12
    nTrials = 200
    depth = 6.0
    durat = 3
    sesMax = np.zeros((nTrials,))
    for i in range(nTrials):
        if np.mod(i,10) == 0:
            print("{0:d}".format(i))
        fluxTS = np.random.randn(2048,1)
 #       fluxTS = np.load('test_ts.npy')
        oFluxTS = fluxTS
        fluxTS[1024:1024+durat] = fluxTS[1024:1024+durat] - depth
#        plt.plot(fluxTS, '.')
#        plt.show()
        searchLen = durat
        varianceFilterFactor = 10
        varianceFilterWindow = searchLen * varianceFilterFactor
        wavObj = waveletObject(waveletLen, fluxTS, varianceFilterWindow)
        trial_pulse = set_trial_transit_pulse(searchLen)
    
        normTS, corrTS = compute_statistic_time_series(wavObj, searchLen, trial_pulse)
        sesMax[i] = np.max(corrTS/normTS)
#        plt.plot(corrTS/normTS, '.')
#        plt.show()
#        plt.plot(1.0/normTS, '.')
#        plt.show()

    print("Mean: {0:f} std: {1:f}".format(np.mean(sesMax), np.std(sesMax)))
    print("Expect Mean: {0:f}".format(depth*np.sqrt(durat)))    
#    matin = sio.loadmat('test_tps_matlabout.mat')
#    matCorrTS = matin['corrTS'].flatten()
#    matNormTS = matin['normTS'].flatten()
#    matSES = matCorrTS/matNormTS.flatten()
#    ses = corrTS/normTS
#    plt.plot(matSES, '.')
#    plt.plot(ses, '.')
#    plt.show()
#    plt.plot(matSES-ses,'.')
#    plt.show()
#    
    
    
