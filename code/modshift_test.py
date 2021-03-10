# -*- coding: utf-8 -*-
"""
Perform a trapezoid model fit to the phased light curve
using the Pretty Good MCMC package
https://github.com/christopherburke/PGMCMC
Although the version on Github is out of date, so use the one included
in this package. Shame on me :(
Also, your not doing MCMC fitting with PGMCMC, just chisquare minimization.

Then run modshift test (via a system call)
The modshift test was written by Jeff Coughlin
https://github.com/JeffLCoughlin/Model-Shift
Modshift is a single c++ file.  The one included in TESS-ExoClass
is slightly modified to return negative values when appropriate.
compile modshift with
g++ -std=c++11 -Wno-unused-result -O3 -o modshift -O modshift.cpp

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
    depth = ioblk.physvals[2]
    bigT = ioblk.physvals[3]
    littleT = ioblk.physvals[4]

    ts = ioblk.normts
    phi = phaseData(ts, per, to)
    lc = np.ones_like(ioblk.normts)
    cadlen = ioblk.parm.cadlen
    samplen = ioblk.parm.samplen
    # Call trapezoid model for data points without any subsampling needed
    idx = np.where(np.logical_and(ioblk.fitdata, ioblk.sampleit == 1))[0]
    if idx.size > 0:
        ztmp = phi[idx] * per
        lctmp = trapezoid(ztmp, depth, bigT, littleT)
        lc[idx] = lctmp
    # Call trapezoid model for data points that need subsampling
    idx = np.where(np.logical_and(ioblk.fitdata, ioblk.sampleit > 1))[0]
    if idx.size > 0:
        ztmp = phi[idx] * per
        deltaXSmall = cadlen / np.float(samplen)
        smallBlock = np.linspace(-cadlen/2.0 + deltaXSmall/2.0,
                                  cadlen/2.0 - deltaXSmall/2.0, samplen)
        oN = ztmp.size
        ztmp_highres = np.tile(ztmp, samplen)
        ztmp_highres = np.reshape(ztmp_highres, (samplen, oN))
        smallBlock_highres = np.tile(smallBlock, oN)
        smallBlock_highres = np.reshape(smallBlock_highres, (oN, samplen))
        smallBlock_highres = np.transpose(smallBlock_highres)
        ztmp_highres = ztmp_highres + smallBlock_highres
        ztmp_highres = ztmp_highres.ravel(order='F')
        lctmp_highres = trapezoid(ztmp_highres, depth, bigT, littleT)
        nN = ztmp_highres.size
        lctmp = lctmp_highres.reshape([oN, nN/oN]).mean(1)
        lc[idx] = lctmp

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
        idx = np.where(ioblk.fitdata)[0]
        residuals = (ioblk.normlc[idx] - ioblk.modellc[idx])/(ioblk.errData[idx])
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
    depth = ioblk.physvals[2]
    bigts = ioblk.physvals[3]
    littlets = ioblk.physvals[4]
    
    ioblk.calcvals[0] = per
    ioblk.calcvals[1] = to
    ioblk.calcvals[2] = depth
    ioblk.calcvals[3] = bigts
    ioblk.calcvals[4] = littlets / bigts
    
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
    depth = ioblk.calcvals[2]
    bigts = ioblk.calcvals[3]
    lilbigratio = ioblk.calcvals[4]
    
    ioblk.physvals[0] = per
    ioblk.physvals[1] = to
    ioblk.physvals[2] = depth
    ioblk.physvals[3] = bigts
    ioblk.physvals[4] = lilbigratio * bigts
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
    #  Directory storing the ses mes time series
    sesMesDir = '/pdo/users/cjburke/spocvet/sector34'
    SECTOR = 34
    OVERWRITE = False
    doPNGs = True
#    pngFolder = '/pdo/users/cjburke/spocvet/sector2/pngs'
    # Run twice once with alt detrend and once with DV median detrend
    medianInputFlux = False
    fileOut = 'spoc_modshift_sector34_20210303.txt'
    #medianInputFlux = True
    #fileOut = 'spoc_modshift_med_sector34_20210303.txt'
    rerun = False    
    if os.path.exists(fileOut) and (not OVERWRITE):
        # Read in previous output to get last TICvalue
        dtypeseq=['i4','i4']
        dtypeseq.extend(['f8']*17)
        dtypeseq.extend(['i4']*6)
        dtypeseq.extend(['f8','f8'])
        dataBlock = np.genfromtxt(fileOut, dtype=dtypeseq)
        modTic = dataBlock['f0']
        lstTic = modTic[-1]
        fom = open(fileOut, 'a+')
        rerun = True
    else:
        fom = open(fileOut, 'w')
    vetFile = 'spoc_fluxtriage_sector34_20210303.txt'
    #vetFile = 'junk.txt'
    tceSeedInFile = 'sector34_20210303_tce.h5'
    
    badTic = np.array([], dtype=np.int64);

    # Load the tce data h5
    tceSeedInFile = 'sector34_20210303_tce.h5'
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=np.int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=np.int)
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
            
    if rerun:
        idx = np.where((alltic == lstTic))[0]
        idxUse = np.arange(idx[0]+1, len(alltic))
        alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
                allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
                allatrpdrstar, allatrpdrstare, allatadrstar = idx_filter(idxUse, \
                alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
                allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
                allatrpdrstar, allatrpdrstare, allatadrstar)
        
    # These lines can be used for debugging
    #idx = np.where((alltic == alltic[530]))[0]
    #alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
    #        allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
    #        allatrpdrstar, allatrpdrstare, allatadrstar = idx_filter(idx, \
    #        alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
    #        allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
    #        allatrpdrstar, allatrpdrstare, allatadrstar)
            
    # Calculate modshift over flux triage passing TCEs
    for i, curTic in enumerate(alltic):
        print('{:d} of {:d} - {:d}'.format(i, len(alltic), curTic))
        curPn = allpn[i]
    
        fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(curTic,curPn))
        f = h5py.File(fileInput,'r')
        altDetrend = np.array(f['altDetrend'])
        validData = np.array(f['validData'])
        time = np.array(f['time'])
        f.close()
        if medianInputFlux:
            fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_dvts_{0:016d}_{1:02d}.h5d'.format(curTic,curPn))
            f = h5py.File(fileInput,'r')
            altDetrend = np.array(f['lc_med_detrend'])+1.0
            print('hello')
            f.close()
            
            
        # Prepare for trapezoid minimization fitting
        # Instantiate pgmcmc_ioblk class and fill in values
        ioblk = pgmcmc_ioblk()
        ioblk.parm.samplen = 15
        ioblk.parm.cadlen = np.abs(np.median(np.diff(time)))
        print('Cadence Length: {:7.4f}'.format(ioblk.parm.cadlen))
        ioblk.parm.fitregion = 4.0
        ioblk.parm.debugLevel = 3
        idxGd = np.where(validData)[0]
        ioblk.normlc = altDetrend[idxGd]
        ioblk.normes = robust.mad(ioblk.normlc)
        ioblk.normts = time[idxGd]
        ioblk.modellc = np.copy(ioblk.normlc)
        ioblk.yData = np.copy(ioblk.normlc)
        ioblk.errData = np.full_like(ioblk.normlc, ioblk.normes)
        durday = alldur[i] / 24.0
        phidur = durday / allper[i]

        # Normalize the time series
        ts = ioblk.normts
        #medianEvent = np.median(np.round((ts - eph)/per))
        #ioblk.timezpt = eph + (medianEvent * per)
        #ioblk.normts = ioblk.normots - ioblk.timezpt
        # identify in transit data to over sample and fitting region
        phi = phaseData(ioblk.normts, allper[i], allatepoch[i])
        ioblk.sampleit = np.where(abs(phi) < (phidur * 1.5), ioblk.parm.samplen, 1)
        ioblk.fitdata = np.where(abs(phi) < (phidur * ioblk.parm.fitregion),\
                             True, False)
        # always fit less than a 0.25 of phase space for stability
        #  and efficiency reasons
        ioblk.fitdata = np.where(abs(phi) > 0.25, False, ioblk.fitdata)

        
        ioblk.physval_names = ['Per', 'To', 'Depth', 'BigT', 'littlet']
        ioblk.calcval_names = ['Per_c', 'To_c', 'Depth_c', 'BigT_c', 'lilBigRatio']
        # Give seed starting values for the minimization
        ioblk.origests = np.array([allper[i], allatepoch[i], allatdep[i]/1.0e6, alldur[i]/24.0, alldur[i]/24.0*0.3])
        # Give integer array for variables you want fixed during fit
        # 0 - not fixed (solved for) ; 1 - fixed (not solved for)
        ioblk.fixed = np.array([1, 0, 0, 0, 0])
        # Give upper and lower limits of variables
        #  Should comfortably constrain values
        #  and in this case will be used as hard constraints in the prior
        #  allowed range 
        delto = alldur[i]/24.0
        ioblk.physval_mins = np.array([allper[i]-0.1, allatepoch[i]-delto, 0.0, alldur[i]/24.0*0.1, 0.0])
        ioblk.physval_maxs = np.array([allper[i]+0.1, allatepoch[i]+delto, allatdep[i]*10.0/1.0e6, alldur[i]/24.0*3.0, alldur[i]/24.0*3.0])
        ioblk.calcval_mins = np.copy(ioblk.physval_mins)
        ioblk.calcval_maxs = np.copy(ioblk.physval_maxs)
        ioblk.calcval_maxs[4] = 1.0

        # Assign the users problem specific functions
        ioblk.func_physvals = pgmcmc_physvals
        ioblk.func_calcvals = pgmcmc_calcvals
        ioblk.func_prior = pgmcmc_prior
        ioblk.func_likehood = pgmcmc_likehood
        ioblk.func_showmodel = []
        ioblk.parm.debugLevel = 0
        ioblk.parm.likehoodmoddisplay = 100
        ioblk.likecount = 0



        # setup some more variables
        ioblk = pgmcmc_setup(ioblk)

        ioblk = pgmcmc_run_minimizer(ioblk, 6)
        
        #print(ioblk.bestphysvals)
        #print(ioblk.bestcalcvals)
        ioblk.physvals = np.copy(ioblk.bestphysvals)
        ioblk.calcvals = np.copy(ioblk.bestcalcvals)
        ioblk, err = pgmcmc_model(ioblk)
        bigT = ioblk.bestphysvals[3]
        trp_per = ioblk.bestphysvals[0]
        phi = phaseData(ioblk.normts, allper[i], ioblk.bestphysvals[1])
        #plt.plot(phi, ioblk.normlc, '.')
        #plt.plot(phi, ioblk.modellc, '.')
        #plt.show()
        
        # Write out the trapezoid fit data and time series in order to run the
        # modshift
        fileOutput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_trpzdfit_{0:016d}_{1:02d}.txt'.format(curTic,curPn))
        fo = open(fileOutput, 'w')
        # Write out the best fit parameters as a header
        #for j in range(len(ioblk.physvals)): 
        #    fo.write('# {:s} {:f}\n'.format(ioblk.physval_names[j], ioblk.bestphysvals[j]))

        for j in range(len(ioblk.normts)):
            strout = '{:f} {:f} {:f}\n'.format(ioblk.normts[j], ioblk.normlc[j]-1.0, ioblk.modellc[j]-1.0)
            fo.write(strout)
        fo.close()
        
        pngOutputPrefix = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_{0:016d}_{1:02d}'.format(curTic,curPn))
        if medianInputFlux:
            pngOutputPrefix = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_{0:016d}_{1:02d}_med'.format(curTic,curPn))
        # Build argument list
        syscall = '/pdo/users/cjburke/spocvet/sector4/modshift {:s} {:s} {:016d}_{:02d} {:f} {:f} 1'.format(\
                            fileOutput, pngOutputPrefix, curTic, curPn, ioblk.bestphysvals[0], ioblk.bestphysvals[1])
        p = Popen(syscall.split(), stdin=None, stdout=PIPE, stderr=PIPE)
        sysreturn, err = p.communicate()
        rc = p.returncode
        #print('alpha')
        #print(sysreturn)
        #print(err)
        #print('alpha')
        #print('data valid')
        #print(np.all(np.isfinite(ioblk.normlc)))
        #print('time valid')
        #print(np.all(np.isfinite(ioblk.normts)))
        #print('model valid')
        #print(np.all(np.isfinite(ioblk.modellc)))
        # Check if there was a bad return from modshift
        if not sysreturn:
            print('Bad modshift return values')
            # Did not get good values from modshift try one hail marry where the period is halved just
            #  to get some results
            #  Add tic to bad list
            badTic = np.append(badTic, curTic)
            syscall = '/pdo/users/cjburke/spocvet/sector4/modshift {:s} {:s} {:016d}_{:02d} {:f} {:f} 1'.format(\
                            fileOutput, pngOutputPrefix, curTic, curPn, ioblk.bestphysvals[0]/2.0, ioblk.bestphysvals[1])
            p = Popen(syscall.split(), stdin=None, stdout=PIPE, stderr=PIPE)
            sysreturn, err = p.communicate()
            rc = p.returncode
            
        if sysreturn:
            retlist = sysreturn.split()[1:]
            primsig = float(retlist[0])
            secsig = float(retlist[1])
            tersig = float(retlist[2])
            possig = float(retlist[3])
            oesig = float(retlist[4])
            depmnmedtest = float(retlist[5])
            shptest = float(retlist[6])
            asymtest = float(retlist[7])
            threshany = float(retlist[8])
            threshdiff = float(retlist[9])
            fred = float(retlist[10])
            primphi = float(retlist[11])
            secphi = float(retlist[12])
            terphi = float(retlist[13])
            posphi = float(retlist[14])
            secdep = float(retlist[15])
            secdeperr = float(retlist[16])
            
            # determine if primary is significant relative to the other detections
            primaryGd = True
            primaryBdRsn = 0
            # Check for fred == 0.0 issue in really bad light curves
            if fred == 0.0:
                fred = 1.0
                primaryGd = False
                primaryBdRsn = 8
            if primsig/fred < threshany:
                primaryGd = False
                primaryBdRsn = 1
            if primsig - tersig < threshdiff:
                primaryGd = False
                primaryBdRsn = primaryBdRsn + 2
            if primsig - possig < threshdiff:
                primaryGd = False
                primaryBdRsn = primaryBdRsn + 4
            # Determine if there is significant secondary
            secondaryGd = True
            secondaryBdRsn = 0
            if secsig/fred < threshany:
                secondaryGd = False
                secondaryBdRsn = 1
            if secsig - tersig < threshdiff:
                secondaryGd = False
                secondaryBdRsn = secondaryBdRsn + 2
            if secsig - possig < threshdiff:
                secondaryGd = False
                seconaryBdRsn = secondaryBdRsn + 4
            # If ther is a significant secondary look for overrided scenarios
            secOvrRid = False
            secOvrRidRsn = 0
            prcalbedo = np.array([99.9, 99.9])
            if secondaryGd:
                # period is incorrect actually shouldbe half
                if primsig - secsig < threshdiff:
                  phidiff = np.abs((primphi+0.5)-secphi)
                  if phidiff < bigT/trp_per*0.25:
                      secOvrRid = True
                      secOvrRidRsn = 1
                # Look for secondary that could be planet albedo detect
                if not secOvrRid:
                    Ntrials = 1000
                    ranrpDrstar = np.random.normal(allatrpdrstar[i], allatrpdrstare[i], size=(Ntrials,))
                    ransecdep = np.random.normal(secdep, secdeperr, size=(Ntrials,))
                    ranrpDa = ranrpDrstar / allatadrstar[i]
                    ranalbedo = ransecdep / ranrpDa / ranrpDa
                    prcalbedo = np.percentile(ranalbedo, [5.0, 50.0])
                    if prcalbedo[0] < 1.0:
                        secOvrRid = True
                        secOvrRidRsn = 2
                        
            
            fom.write('{:016d} {:02d} {:s} {:d} {:d} {:d} {:d} {:d} {:d} {:f} {:f}\n'.format( \
                          curTic, curPn, ' '.join(sysreturn.split()[1:]), int(primaryGd), primaryBdRsn, \
                          int(secondaryGd), secondaryBdRsn, int(secOvrRid), secOvrRidRsn, \
                          prcalbedo[0], prcalbedo[1]))
        else:
            print('Still Bad TIC: {0:d}'.format(curTic))
            badTic = np.append(badTic, curTic)
        
        
        

    fom.close()

    print('There were {0:d} Bad TICS'.format(len(badTic)))
    print(badTic)