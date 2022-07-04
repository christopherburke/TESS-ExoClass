# -*- coding: utf-8 -*-
"""
This routine implements the SWEET test from Kepler Robovetter
Uses the PGMCMC package for fitting of sinusoid.

AUTHOR: Christopher J. Burke
"""

import numpy as np
import pickle
from gather_tce_fromdvxml import tce_seed
import os
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


def pgmcmc_model(ioblk):
    """Generate a subsampled model at the current parameters
       INPUT:
       ioblk - [class] trp_ioblk class structure
       OUTPUT:
       ioblk - [class] modified ioblk
       err - [0 ok; 1 not ok] Error flag
    """
    err = 0
    per = ioblk.physvals[0]
    to = ioblk.physvals[1]
    amp = ioblk.physvals[2]
    zpt = ioblk.physvals[3]

    ts = ioblk.normts
    lc = np.zeros_like(ioblk.normts) + zpt
    lc = lc + amp*np.sin((ts-to)/per*2.0*np.pi)

    ioblk.modellc = lc
    if np.sum(np.isfinite(lc)) != lc.size:
        err = 1
    return ioblk, err

def pgmcmc_likehood(ioblk):
    """Return a residual time series of data minus model
       trp_setup(ioblk) should be called before this function is called
       INPUT:
       ioblk - [class] trp_ioblk class structure
       OUTPUT:
       ioblk - [class] modified ioblk
    """
    err = 0 # Error flag
    ioblk.likecount += 1
    #ioblk.likecount += 1
    #if ioblk.likecount > 75:
    #    exit()
    # Users model
    ioblk, err = pgmcmc_model(ioblk)
    
    # Model is valid calculate likelihood
    if (err == 0):
    # Calculate residuals
        residuals = (ioblk.normlc - ioblk.modellc)/(ioblk.errData)
        # Return scalar summed residuals
        ioblk.mcmc.chi2 = np.sum(residuals**2)
        ioblk.mcmc.like = np.log10(ioblk.mcmc.chi2)
        if (ioblk.mcmc.chi2 < ioblk.chi2min):
            if ioblk.parm.debugLevel>0:
                print("New Best Like: {0:f} Old: {1:f} Expected: {2:f}".format( \
                      ioblk.mcmc.chi2, ioblk.chi2min, ioblk.expchi2))
            ioblk.chi2min = ioblk.mcmc.chi2
            ioblk.bestphysvals = ioblk.physvals
            ioblk.bestcalcvals = ioblk.calcvals
    
    #if (ioblk.parm.debugLevel > 2 and \
    #    ioblk.mcmc.curpar == 0 and \
    #    ioblk.pt.curtempidx == ioblk.pt.ntemp - 1 and \
    #    np.mod(ioblk.mcmc.pos,ioblk.parm.likehoodmoddisplay) == 0):
        
        #pgmcmc_showmodel(ioblk)
        
    return ioblk, err

 
    
    


def pgmcmc_calcvals(ioblk):
    """Convert physical parameters to calculation variables
       This is where you scale physical variables or more generic
         variable transformation
       In this example we are not doing any variable transformations
       INPUT:
         ioblk - [class] pgmcmc_ioblk class
       OUTPUT: 
         ioblk - [class]
         err - [0 ok ; 1 not ok]
    """
    err = 0 # Error flag
    per = ioblk.physvals[0]
    to = ioblk.physvals[1]
    amp = ioblk.physvals[2]
    zpt = ioblk.physvals[3]
    
    ioblk.calcvals[0] = per
    ioblk.calcvals[1] = to
    ioblk.calcvals[2] = amp
    ioblk.calcvals[3] = zpt
    
    #ioblk.calcvals = np.copy(ioblk.physvals)
    if ~np.isfinite(ioblk.calcvals).all():
        print("Calculation Values Bad")
        print(ioblk.calcvals)
        print(ioblk.physvals)                 
        err = 1
    return ioblk, err

def pgmcmc_physvals(ioblk):
    """Convert calculation variables to physical variables
    This is where you undo the variable transformation performed in
    pgmcmc_calcvals.
    In this example there is not any variables transformations
       INPUT:
         ioblk - [class] pgmcmc_ioblk class
       OUTPUT: 
         ioblk - [class]
         err - [0 ok ; 1 not ok]
    """
    err = 0 # Error flag
    per = ioblk.calcvals[0]
    to = ioblk.calcvals[1]
    amp = ioblk.calcvals[2]
    zpt = ioblk.calcvals[3]
    
    ioblk.physvals[0] = per
    ioblk.physvals[1] = to
    ioblk.physvals[2] = amp
    ioblk.physvals[3] = zpt
    #ioblk.physvals = np.copy(ioblk.calcvals)
    if ~np.isfinite(ioblk.physvals).all():                 
        print("Physical Values Bad")
        print(ioblk.calcvals)
        print(ioblk.physvals)                 
        err = 1
    return ioblk, err

def pgmcmc_prior(ioblk):
    """Calculate prior and also check variable bounds
       This function must catch all variable values that would
       cause the model or likelihood to fail, nan, imaginary, or infs
    """
    err = 0 # Error flag
    prior = 0.0
    ioblk.calcvals[ioblk.mcmc.paridx] = ioblk.mcmc.pars
    # First the calcvals must be in bounds
    anybadcheck = np.any(np.logical_or(ioblk.calcvals < ioblk.calcval_mins, \
                         ioblk.calcvals > ioblk.calcval_maxs))
    if (anybadcheck):
        err = 1
    else:
        ioblk, err = pgmcmc_physvals(ioblk)
        # physvals must be in bounds
        anybadcheck2 = np.any(np.logical_or(ioblk.physvals < ioblk.physval_mins, \
                              ioblk.physvals > ioblk.physval_maxs))
        if (anybadcheck2):
            err = 1

    # If variables are within bounds then proceed to calculate prior
    if (err == 0):
        priorvals = np.log(1.0 / (ioblk.physval_maxs - ioblk.physval_mins))
        prior = np.sum(priorvals)
        # Warn on nonfinite and nonreal prior 
        if ~np.isfinite(prior):
            print("Non Finite Prior")
            print(ioblk.physvals)
        if ~np.isreal(prior):
            print("Non Real Prior")
            print(ioblk.physvals)
    ioblk.mcmc.prior = prior
    return ioblk, err    


if __name__ == '__main__':
    # Maximum period for test
    MAXPER = 5.0
    
    #  Directory storing the ses mes time series
    sesMesDir = '/pdo/users/cjburke/spocvet/sector52'
    SECTOR = 52

    fileOut = 'spoc_sweet_sector-52.txt'
    fom = open(fileOut, 'w')
    vetFile = 'spoc_fluxtriage_sector-52.txt'
    tceSeedInFile = 'sector-52_tce.h5'

    # Load the tce data h5
    tceSeedInFile = 'sector-52_tce.h5'
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
    # and flux vetted pass and period < PERMAX
    idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1) & (allper<MAXPER))[0]
    alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
            allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
            allatrpdrstar, allatrpdrstare, allatadrstar = idx_filter(idx, \
            alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
            allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
            allatrpdrstar, allatrpdrstare, allatadrstar)
            
    # Calculate modshift over flux triage passing TCEs
    for i, curTic in enumerate(alltic):
        print('{:d} of {:d}'.format(i, len(alltic)))
        curPn = allpn[i]
        fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(curTic,curPn))
        f = h5py.File(fileInput,'r')
        flx = np.array(f['initFlux'])
        validData = np.array(f['validData'])
        time = np.array(f['time'])
        f.close()
            
        # Prepare for sin curve fitting
        # Instantiate pgmcmc_ioblk class and fill in values
        ioblk = pgmcmc_ioblk()
        ioblk.parm.samplen = 15
        idxGd = np.where(validData)[0]
        ioblk.parm.cadlen = np.abs(np.median(np.diff(time[idxGd])))
        print('Cadence Length: {:7.4f}'.format(ioblk.parm.cadlen))
        ioblk.parm.fitregion = 4.0
        ioblk.parm.debugLevel = 3
        tmpflx = flx[idxGd]
        tmpt = time[idxGd]
        tmpttwo = tmpt*tmpt
        tmptfour = tmpttwo*tmpttwo
        # Remove polynomial fit
        pvals  = np.polyfit(tmpt,tmpflx,4,full=False)
        tmpy = pvals[4] + pvals[3]*tmpt + pvals[2]*tmpttwo + pvals[1]*tmpttwo*tmpt + pvals[0]*tmptfour
        #plt.plot(tmpt, tmpflx, '.')
        #plt.plot(tmpt, tmpy, '-')
        #plt.show()
        ioblk.normlc = tmpflx/tmpy - 1.0
        
        ioblk.normes = robust.mad(ioblk.normlc)
        origstd = np.copy(ioblk.normes)
        ioblk.normts = time[idxGd]
        ioblk.modellc = np.copy(ioblk.normlc)
        ioblk.yData = np.copy(ioblk.normlc)
        ioblk.errData = np.full_like(ioblk.normlc, ioblk.normes)

        ioblk.timezpt = np.median(ioblk.normts)
        ioblk.normts = ioblk.normts - ioblk.timezpt

        
        ioblk.physval_names = ['Per', 'To', 'Amp', 'Zpt']
        ioblk.calcval_names = ['Per_c', 'To_c', 'Amp_c', 'Zpt_c']
        # Give seed starting values for the minimization
        ioblk.origests = np.array([allper[i], 0.0, ioblk.normes*3.0, 0.0])
        # Give integer array for variables you want fixed during fit
        # 0 - not fixed (solved for) ; 1 - fixed (not solved for)
        ioblk.fixed = np.array([1, 0, 0, 0])
        # Give upper and lower limits of variables
        #  Should comfortably constrain values
        #  and in this case will be used as hard constraints in the prior
        #  allowed range 
        ioblk.physval_mins = np.array([allper[i]-0.1, -allper[i], 0.0, -ioblk.normes*2.0])
        ioblk.physval_maxs = np.array([allper[i]+0.1, allper[i], ioblk.normes*20.0, ioblk.normes*2.0])
        ioblk.calcval_mins = np.copy(ioblk.physval_mins)
        ioblk.calcval_maxs = np.copy(ioblk.physval_maxs)

        # Assign the users problem specific functions
        ioblk.func_physvals = pgmcmc_physvals
        ioblk.func_calcvals = pgmcmc_calcvals
        ioblk.func_prior = pgmcmc_prior
        ioblk.func_likehood = pgmcmc_likehood
        ioblk.func_showmodel = []
        ioblk.parm.debugLevel = 0
        ioblk.parm.likehoodmoddisplay = 100
        ioblk.parm.outPrefix = 'sweet'        
        ioblk.likecount = 0

        # setup some more variables
        ioblk = pgmcmc_setup(ioblk)

        ioblk = pgmcmc_run_minimizer(ioblk, 6)
        
        #print(ioblk.bestphysvals)
        #print(ioblk.bestcalcvals)
        ioblk.physvals = np.copy(ioblk.bestphysvals)
        ioblk.calcvals = np.copy(ioblk.bestcalcvals)
        ioblk, err = pgmcmc_model(ioblk)
        onestd = robust.mad(ioblk.normlc-ioblk.modellc)
        #print(origstd, onestd, onestd/origstd)
        #plt.plot(ioblk.normts, ioblk.normlc, '.')
        #plt.plot(ioblk.normts, ioblk.modellc, '.')
        #plt.show()
        
        pervalues = np.copy(ioblk.bestphysvals)
        
        # Do half period
        ioblk.origests = np.array([allper[i]/2.0, 0.0, ioblk.normes*3.0, 0.0])
        ioblk.physval_mins = np.array([allper[i]/2.0-0.1, -allper[i], 0.0, -ioblk.normes*2.0])
        ioblk.physval_maxs = np.array([allper[i]/2.0+0.1, allper[i], ioblk.normes*20.0, ioblk.normes*2.0])
        ioblk.calcval_mins = np.copy(ioblk.physval_mins)
        ioblk.calcval_maxs = np.copy(ioblk.physval_maxs)
        ioblk.likecount = 0
        # setup some more variables
        ioblk = pgmcmc_setup(ioblk)

        ioblk = pgmcmc_run_minimizer(ioblk, 6)
        
        #print(ioblk.bestphysvals)
        #print(ioblk.bestcalcvals)
        ioblk.physvals = np.copy(ioblk.bestphysvals)
        ioblk.calcvals = np.copy(ioblk.bestcalcvals)
        ioblk, err = pgmcmc_model(ioblk)
        #print('half')
        halfstd = robust.mad(ioblk.normlc-ioblk.modellc)
        #print(origstd, halfstd, halfstd/origstd)

        #plt.plot(ioblk.normts, ioblk.normlc, '.')
        #plt.plot(ioblk.normts, ioblk.modellc, '.')
        #plt.show()

        halfvalues = np.copy(ioblk.bestphysvals)

        # Do double period
        ioblk.origests = np.array([allper[i]*2.0, 0.0, ioblk.normes*3.0, 0.0])
        ioblk.physval_mins = np.array([allper[i]*2.0-0.1, -allper[i], 0.0, -ioblk.normes*2.0])
        ioblk.physval_maxs = np.array([allper[i]*2.0+0.1, allper[i], ioblk.normes*20.0, ioblk.normes*2.0])
        ioblk.calcval_mins = np.copy(ioblk.physval_mins)
        ioblk.calcval_maxs = np.copy(ioblk.physval_maxs)
        ioblk.likecount = 0
        # setup some more variables
        ioblk = pgmcmc_setup(ioblk)

        ioblk = pgmcmc_run_minimizer(ioblk, 6)
        
        #print(ioblk.bestphysvals)
        #print(ioblk.bestcalcvals)
        ioblk.physvals = np.copy(ioblk.bestphysvals)
        ioblk.calcvals = np.copy(ioblk.bestcalcvals)
        ioblk, err = pgmcmc_model(ioblk)
        #print('double')
        doublestd = robust.mad(ioblk.normlc-ioblk.modellc)
        #print(origstd, doublestd, doublestd/origstd)

        #plt.plot(ioblk.normts, ioblk.normlc, '.')
        #plt.plot(ioblk.normts, ioblk.modellc, '.')
        #plt.show()

        doublevalues = np.copy(ioblk.bestphysvals)
        
        
        fom.write('{:016d} {:02d} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f} {:f}\n'.format( \
                      curTic, curPn, pervalues[0], pervalues[1], pervalues[2], pervalues[3], \
                      halfvalues[0], halfvalues[1], halfvalues[2], halfvalues[3], \
                      doublevalues[0], doublevalues[1], doublevalues[2], doublevalues[3], \
                      origstd, onestd, halfstd, doublestd, np.min([onestd/origstd, halfstd/origstd, doublestd/origstd])))
        
        
        
        

    fom.close()