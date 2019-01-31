""" Module for Pretty Good MCMC parameter estimation
    Heavily based upon Philip Gregory's awesome Bayesian book
    Author: Christopher J Burke 
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pickle
import copy
import scipy.optimize as opt


class pgmcmc_ioblk:
    """Define a class that contains all the data needed to perform
       an mcmc parameter estimate.  Numerous functions will use as input this
       class.  This is purely a storage class
       See pgmcmc_setup to illustrate how to iinitialize these
       storage classes
       CONTENTS:
       parm - [class] Storage class pgmcmc_parameters for algorithm parameters
       mcmc - [class] Storage class pgmcmc_mcmc for runtime variables
       pt - [class] Storage class pgmcmc_pt for parallel tempering runtime vars
    """
    def __init__(self):
        self.parm = pgmcmc_parameters()
        self.mcmc = pgmcmc_mcmc()
        self.pt = pgmcmc_pt()
        self.origests = np.array([0.0])
        self.physval_names = ['']
        self.calcval_names = ['']
        self.fixed = np.array([0])
        self.nparm = 0

        self.physval_mins = np.array([0.0])
        self.physval_maxs = np.array([0.0])
        self.physvals = np.array([0.0])
        self.physvalsavs = np.array([0.0])
        self.bestphysvals = np.array([0.0])
        self.calcvals = np.array([0.0])
        self.calcvalsavs = np.array([0.0])
        self.bestcalcvals = np.array([0.0])
        self.calcval_mins = np.array([0.0])
        self.calcval_maxs = np.array([0.0])
        self.scls = np.array([0.0])
        self.yModel = np.array([0.0])
        self.chi2min = 0.0
        self.yData = np.array([0.0])
        self.errData = np.array([0.0])
        self.xData = np.array([0.0])
        self.likenorm = 0.0
        self.expchi2 = 1.0
        self.func_physvals = []
        self.func_calcvals = []
        self.func_prior = []
        self.func_likehood = []
        self.func_showmodel = []
        self.fighandle = []
        self.axhandle = []
        
    #def __str__(self):
    #    for k in self.__dict__:
    #        print(k, self.__dict__[k])
    #    return ''

def pgmcmc_one_mcmc_step(ioblk):
    """ Take one iteration through all parameters one at a time
        this makes up one mcmc step
    """
    curp = np.copy(ioblk.calcvals[ioblk.mcmc.paridx])
    sigs = ioblk.scls[ioblk.mcmc.paridx]
    curchi2 = np.copy(ioblk.mcmc.like)
    curprior = np.copy(ioblk.mcmc.prior)
    nexp = np.copy(curp)
    tempfactor = 1.0
    if ioblk.parm.dopartemp:
        tempfactor = ioblk.pt.curtemp
        
    # iterate through each free to vary parameter
    for ip in np.arange(curp.size):
        ioblk.mcmc.attempts[ip] += 1
        ioblk.mcmc.curpar = ioblk.mcmc.paridx[ip]
        gdstp = 0
        # Jump in this parameter only
        nexp[ip] = np.random.randn() * sigs[ip] + curp[ip]
        ioblk.mcmc.pars = np.copy(nexp)
        # Check prior for new set of parameters
        ioblk, err = ioblk.func_prior(ioblk)
        if err == 0:
            # Passes prior now calculate likelihood of new parameters
            ioblk, err = ioblk.func_likehood(ioblk)
            ioblk.mcmc.like *= tempfactor
            nexchi2 = np.copy(ioblk.mcmc.like)
            nexprior = np.copy(ioblk.mcmc.prior)
            rat = nexprior - curprior + nexchi2 - curchi2
            if (rat >= 0.0): # Pass by large margin
                curp = nexp
                ioblk.mcmc.accepts[ip] += 1
                gdstp = 1
            elif (np.log(np.random.rand()) < rat): # check if randomly passes
                curp = nexp
                ioblk.mcmc.accepts[ip] += 1
                gdstp = 1
            else: # Did not pass repeat previous parameter values
                nexchi2 = curchi2
                nexprior = curprior
                nexp = curp
        else: # Did not pass prior or prior failed repeat previous parameters
            nexchi2 = curchi2
            nexprior = curprior
            nexp = curp
        curchi2 = np.copy(nexchi2)
        curprior = np.copy(nexprior)
        curp = np.copy(nexp)
    if gdstp == 0: # That last parameter did not pass. update ioblk
        ioblk.mcmc.like = np.copy(curchi2)
        ioblk.mcmc.prior = np.copy(curprior)
        ioblk.mcmc.pars = np.copy(curp)
        ioblk.calcvals[ioblk.mcmc.paridx] = np.copy(curp)
        ioblk, err = ioblk.func_physvals(ioblk)
    return ioblk
                
def pgmcmc_burnit(ioblk, nBurn):
    """ Burn though several chain steps just returning
        acceptance fraction and ignoring parameter outputs
    """
    ioblk.mcmc.accepts = ioblk.mcmc.accepts * 0
    ioblk.mcmc.attempts = ioblk.mcmc.attempts * 0
    for istp in range(nBurn):
        ioblk = pgmcmc_one_mcmc_step(ioblk)
    fracs = np.double(ioblk.mcmc.accepts) / \
            np.double(ioblk.mcmc.attempts)
    return ioblk, fracs
    
def pgmcmc_iterate_proposals(ioblk):
    """ Automatically determine the proposal step sizes
        This uses a gradual relaxation of values to ensure
        the fraction of accepted jumps is within a desired
        range.  It starts with a coarse/large adjustments
        of the step sizes then has a round of refined/small
        adjustments.  Parameters that run this routine are
        in pgmcmc_parameters class
    """
    if (ioblk.parm.debugLevel > 2): ioblk.func_showmodel(ioblk)

    # Do initial burn
    print("Start Initial Burn")
    ioblk, fracs = pgmcmc_burnit(ioblk, ioblk.parm.initNSteps)
    print("Initial Burn Finished")
    print(fracs)
    #print ioblk.scls[ioblk.mcmc.paridx]

    if (ioblk.parm.debugLevel > 2): ioblk.func_showmodel(ioblk)

    # Start Coarse Proposal Step size corrections
    converg = 0
    i = 0
    while converg == 0 and i < ioblk.parm.maxPropTries :
        idxlow = np.where(fracs < ioblk.parm.coarseLowFrac)[0]
        idxhgh = np.where(fracs > ioblk.parm.coarseHghFrac)[0]
        if idxlow.size + idxhgh.size == 0:
            converg = 1
        else:
            if (not idxlow.size == 0):
                j = ioblk.mcmc.paridx[idxlow]
                ioblk.scls[j] = ioblk.scls[j] / 2.0
            if (not idxhgh.size == 0):
                j = ioblk.mcmc.paridx[idxhgh]
                ioblk.scls[j] = ioblk.scls[j] * 2.0
            i += 1
            ioblk, fracs = pgmcmc_burnit(ioblk, ioblk.parm.coarseNSteps)
            print("Coarse iteration: ", i)
            print(fracs)
            print(ioblk.scls[ioblk.mcmc.paridx])
            if (ioblk.parm.debugLevel > 2): ioblk.func_showmodel(ioblk)
    print("Done Coarse proposal iteration")
    if (i == ioblk.parm.maxPropTries):
        print("Coarse Proposal Iteration did not converge in Steps: ", \
                ioblk.parm.maxPropTries)
        input("Please Ctrl-C to exit")


    # Start Refined Proposal Step Size corrections
    converg = 0
    i = 0
    while converg == 0 and i < ioblk.parm.maxPropTries :
        idxlow = np.where(fracs < ioblk.parm.refineLowFrac)[0]
        idxhgh = np.where(fracs > ioblk.parm.refineHghFrac)[0]
        if idxlow.size + idxhgh.size == 0:
            converg = 1
        else:
            if (not idxlow.size == 0):
                j = ioblk.mcmc.paridx[idxlow]
                ioblk.scls[j] = ioblk.scls[j] / (1.0 + \
                    np.abs(fracs[idxlow] - ioblk.parm.fracwant))
            if (not idxhgh.size == 0):
                j = ioblk.mcmc.paridx[idxhgh]
                ioblk.scls[j] = ioblk.scls[j] * (1.0 + \
                    np.abs(fracs[idxhgh] - ioblk.parm.fracwant))
            i += 1
            ioblk, fracs = pgmcmc_burnit(ioblk, ioblk.parm.refineNSteps)
            print("Refine iteration: ", i)
            print(fracs)
            print(ioblk.scls[ioblk.mcmc.paridx])
            if (ioblk.parm.debugLevel > 2): ioblk.func_showmodel(ioblk)
    print("Done Refined proposal iteration")
    if (i == ioblk.parm.maxPropTries):
        print("Refined Proposal Iteration did not converge in Steps: ", \
                ioblk.parm.maxPropTries)
        input("Please Ctrl-C to exit")
    
    pgmcmc_save_state(ioblk,'ioblk_props')
    
    return ioblk

def pgmcmc_check_pt_swap(ioblk):
    """Swap state between a randomly selected temperature values
    """
    
    # See if we should even attempt PT swap
    if np.random.rand() <= 1.0 / ioblk.parm.avgpt:
        # Try swap
        ioblk.mcmc.curpar = 0
        ioblk.pt.nswaps += 1
        ioblk.pt.swapat[ioblk.pt.nswaps] = ioblk.mcmc.pos
        ioblk.pt.swapgd[ioblk.pt.nswaps] = 0
        # Which Temp index to tray and swap
        wantt = int(np.floor(np.random.rand() * (ioblk.pt.ntemp-1)))
        ioblk.pt.swaptemp[ioblk.pt.nswaps] = wantt
        ioblk.pt.swapattempts[wantt] += 1
        # Get data where wantt is currently
        ioblk.pt.curtempidx = wantt
        ioblk = pgmcmc_loadtemp(ioblk)
        curp0 = np.copy(ioblk.calcvals[ioblk.mcmc.paridx])
        temp0 = np.copy(ioblk.pt.curtemp)
        ioblk.mcmc.pars = curp0
        ioblk, err = ioblk.func_prior(ioblk)
        ioblk, err = ioblk.func_likehood(ioblk)
        ioblk.mcmc.like *= temp0
        curchi200 = np.copy(ioblk.mcmc.like)
        curprior00 = np.copy(ioblk.mcmc.prior)
        # Get data where wantt+1 is currently
        ioblk.pt.curtempidx = wantt+1
        ioblk = pgmcmc_loadtemp(ioblk)
        curp1 = np.copy(ioblk.calcvals[ioblk.mcmc.paridx])
        temp1 = np.copy(ioblk.pt.curtemp)
        ioblk.mcmc.pars = curp1
        ioblk, err = ioblk.func_prior(ioblk)
        ioblk, err = ioblk.func_likehood(ioblk)
        ioblk.mcmc.like *= temp1
        curchi211 = np.copy(ioblk.mcmc.like)
        curprior11 = np.copy(ioblk.mcmc.prior)
        # Get data with wantt where wantt+1 is located
        ioblk.pt.curtempidx = wantt
        ioblk = pgmcmc_loadtemp(ioblk)
        ioblk.mcmc.pars = curp1
        ioblk, err = ioblk.func_prior(ioblk)
        ioblk, err = ioblk.func_likehood(ioblk)
        ioblk.mcmc.like *= temp0
        curchi201 = np.copy(ioblk.mcmc.like)
        curprior01 = np.copy(ioblk.mcmc.prior)
        # Get data with wantt+1 where wantt is located
        ioblk.pt.curtempidx = wantt+1
        ioblk = pgmcmc_loadtemp(ioblk)
        ioblk.mcmc.pars = curp0
        ioblk, err = ioblk.func_prior(ioblk)
        ioblk, err = ioblk.func_likehood(ioblk)
        ioblk.mcmc.like *= temp1
        curchi210 = np.copy(ioblk.mcmc.like)
        curprior10 = np.copy(ioblk.mcmc.prior)
        
        # Check whether to actually accept swap
        ratio = (curchi201+curprior01) + (curchi210+curprior10) - \
                (curchi200+curprior00) - (curchi211+curprior11)
        ratio = np.min([0.0, ratio])
        if np.log(np.random.rand()) < ratio: 
            # Successful Swap!
            ioblk.pt.swapaccepts[wantt] += 1
            ioblk.pt.swapgd[ioblk.pt.nswaps] = 1
            ioblk.pt.curtempidx = wantt
            ioblk = pgmcmc_loadtemp(ioblk)
            tempfactor = ioblk.pt.curtemp
            ioblk.mcmc.pars = curp1
            ioblk, err = ioblk.func_prior(ioblk)
            ioblk, err = ioblk.func_likehood(ioblk)
            ioblk.mcmc.like *= tempfactor
            ioblk = pgmcmc_savetemp(ioblk)
            ioblk.pt.curtempidx = wantt + 1
            ioblk = pgmcmc_loadtemp(ioblk)
            tempfactor = ioblk.pt.curtemp
            ioblk.mcmc.pars = curp0
            ioblk, err = ioblk.func_prior(ioblk)
            ioblk, err = ioblk.func_likehood(ioblk)
            ioblk.mcmc.like *= tempfactor
            ioblk = pgmcmc_savetemp(ioblk)
    return ioblk
    
def pgmcmc_run_mcmc(ioblk, postPropStart=False):    

    # reset things
    ioblk.mcmc.pos = 1 # a value of 1 turns off plotting
    ioblk.mcmc.curpar = 0
    ioblk.mcmc.attempts = ioblk.mcmc.attempts * 0
    ioblk.mcmc.accepts = ioblk.mcmc.accepts * 0

    if not postPropStart:
        # Start proposal step size iteration across all temperatures
        for t in range(ioblk.pt.ntemp):
            ioblk.pt.curtempidx=t
            ioblk = pgmcmc_loadtemp(ioblk)
            print("Start Proposal Iterations for Temperature: {0:d} {1:f}".format(t, ioblk.pt.curtemp))
            if t > 0 and ioblk.pt.temps[0] != 0.0:
                ioblk.scls = np.copy(ioblk.pt.allscls[:, t-1] * ioblk.pt.temps[t-1] \
                                / ioblk.pt.temps[t])
                ioblk.calcvals = np.copy(ioblk.pt.allcalcval[:, t-1])
                ioblk.mcmc.like = np.copy(ioblk.pt.alllike[t-1] * ioblk.pt.temps[t] \
                                / ioblk.pt.temps[t-1])
                ioblk.mcmc.prior = np.copy(ioblk.pt.allprior[t-1])
            ioblk = pgmcmc_iterate_proposals(ioblk)
            ioblk = pgmcmc_savetemp(ioblk)
        print("Done with Proposal Iterations")
    else:
        # Load a set of proposal steps from a previous run
        previoblk = pgmcmc_load_state('ioblk_props')
        if ioblk.parm.debugLevel > 2:
            fh = ioblk.fighandle
            ah = ioblk.axhandle
        func1 = ioblk.func_calcvals
        func2 = ioblk.func_physvals
        func3 = ioblk.func_likehood
        func4 = ioblk.func_prior
        func5 = ioblk.func_showmodel
        ioblk = copy.deepcopy(previoblk)
        
        if ioblk.parm.debugLevel > 2:
            ioblk.fighandle = fh
            ioblk.axhandle = ah
        ioblk.func_calcvals = func1
        ioblk.func_physvals = func2
        ioblk.func_likehood = func3
        ioblk.func_prior = func4
        ioblk.func_showmodel = func5
        print("Loaded a previous proposal iteration")
    print("Start MCMC Run")

    # reset things
    ioblk.mcmc.pos = 0
    ioblk.mcmc.curpar = 0
    ioblk.mcmc.attempts = ioblk.mcmc.attempts * 0
    ioblk.mcmc.accepts = ioblk.mcmc.accepts * 0
    # Make storage data arrays
    pvals = np.zeros((ioblk.parm.maxstps, ioblk.mcmc.paridx.size))
    cvals = np.copy(pvals)
    bvals = np.zeros((ioblk.parm.maxstps, 3))
    if ioblk.parm.saveAltTemp:
        altpvals = np.zeros((ioblk.parm.maxstps, ioblk.mcmc.paridx.size))
        altcvals = np.copy(altpvals)
        altbvals = np.zeros((ioblk.parm.maxstps, 3))

    for i in range(ioblk.parm.maxstps):
        for t in range(ioblk.pt.ntemp):
            ioblk.pt.curtempidx = t
            if (ioblk.parm.dopartemp):
                ioblk = pgmcmc_loadtemp(ioblk)
                ioblk = pgmcmc_one_mcmc_step(ioblk)
                ioblk = pgmcmc_savetemp(ioblk)
                # See if we want to save a particular temperatures result
                if ioblk.parm.saveAltTemp and t == ioblk.parm.savedTempIdx:
                           altpvals[i,:] = ioblk.physvals[ioblk.mcmc.paridx]
                           altcvals[i,:] = ioblk.calcvals[ioblk.mcmc.paridx]
                           altbvals[i,0] = ioblk.mcmc.like
                           altbvals[i,1] = ioblk.mcmc.prior
                           altbvals[i,2] = ioblk.mcmc.chi2 
            else:
                ioblk = pgmcmc_one_mcmc_step(ioblk)
        # On coolest temperate save result
        pvals[i,:] = ioblk.physvals[ioblk.mcmc.paridx]
        cvals[i,:] = ioblk.calcvals[ioblk.mcmc.paridx]
        bvals[i,0] = ioblk.mcmc.like
        bvals[i,1] = ioblk.mcmc.prior
        bvals[i,2] = ioblk.mcmc.chi2
        # Now attempt temperature swap
        if (ioblk.pt.ntemp > 1 and ioblk.parm.dopartemp):
            ioblk = pgmcmc_check_pt_swap(ioblk)
            
        if np.mod(i, ioblk.parm.saveNSteps) == 0:
            pgmcmc_save_state(ioblk, 'run', pvals, cvals, bvals)
            if ioblk.parm.saveAltTemp:
                pgmcmc_save_state(ioblk, 'runalt', altpvals, altcvals, altbvals)                
        if np.mod(i, 100) == 0:
            print("Step: ", i)
        ioblk.mcmc.pos += 1
        
    pgmcmc_save_state(ioblk, 'run', pvals, cvals, bvals)
    if ioblk.parm.saveAltTemp:
        pgmcmc_save_state(ioblk, 'runalt', altpvals, altcvals, altbvals)
    return
        
def pgmcmc_save_state(ioblk, prefix, pvals=[], cvals=[], bvals=[]):
    # Handles cannot be pickled so save them, undefine them, then put back
    fh = ioblk.fighandle
    ah = ioblk.axhandle
    func1 = ioblk.func_calcvals
    func2 = ioblk.func_physvals
    func3 = ioblk.func_likehood
    func4 = ioblk.func_prior
    func5 = ioblk.func_showmodel
    ioblk.fighandle = []
    ioblk.axhandle = []
    ioblk.func_calcvals = []
    ioblk.func_physvals = []
    ioblk.func_likehood = []
    ioblk.func_prior = []
    ioblk.func_showmodel =[]
    
    # These are specific to my occurrence rate problem
    if hasattr(ioblk, 'cudaper2d'):
        d1 = ioblk.cudaper2d
        d2 = ioblk.cudarp2d
        d3 = ioblk.cudasfuncret
        d4 = ioblk.hostsfuncret
        ioblk.cudaper2d = []
        ioblk.cudarp2d = []
        ioblk.cudasfuncret = []
        ioblk.hostsfuncret = []
    if hasattr(ioblk, 'cudaplansamps'):
        e1 = ioblk.cudaplansamps 
        e2 = ioblk.cudaplansamppers 
        e3 = ioblk.cudaplansfuncret 
        e4 = ioblk.cudaplanhostsfuncret 
        ioblk.cudaplansamps = []
        ioblk.cudaplansamppers = []
        ioblk.cudaplansfuncret  = []
        ioblk.cudaplanhostsfuncret = []
        
       
    pickle.dump(ioblk, open(prefix+'.pkl', 'wb'), protocol=0)#, fix_imports=True)
    ioblk.fighandle = fh
    ioblk.axhandle = ah
    ioblk.func_calcvals = func1
    ioblk.func_physvals = func2
    ioblk.func_likehood = func3
    ioblk.func_prior = func4
    ioblk.func_showmodel = func5
   
    if hasattr(ioblk, 'cudaper2d'):
        ioblk.cudaper2d = d1
        ioblk.cudarp2d = d2
        ioblk.cudasfuncret = d3
        ioblk.hostsfuncret = d4
    if hasattr(ioblk, 'cudaplansamps'):
        ioblk.cudaplansamps = e1
        ioblk.cudaplansamppers = e2
        ioblk.cudaplansfuncret  = e3
        ioblk.cudaplanhostsfuncret = e4
     
    f = h5py.File(prefix+'.hd5', 'w')
    if (not len(pvals) == 0):
        pvaldset = f.create_dataset('Pvals', data = pvals, compression='gzip')
    if (not len(cvals) == 0):
        cvaldset = f.create_dataset('Cvals', data = cvals, compression='gzip')
    if (not len(bvals) == 0):
        bvaldset = f.create_dataset('Bvals', data = bvals, compression='gzip')
    f.close()
    
    return
    
def pgmcmc_load_state(prefix):
    
    ioblk = pickle.load(open(prefix+'.pkl', 'rb'))
    return ioblk
    
def pgmcmc_setup(ioblk):
    """Setup various data products before minimizing
       INPUT:
       ioblk - [class] pgmcmc_ioblk class structure
       OUTPUT: 
       ioblk - [class] modified ioblk
    """
    ioblk.physvals = np.copy(ioblk.origests)
    ioblk.calcvals = np.copy(ioblk.origests)
    ioblk.nparm = np.size(ioblk.fixed)

    ioblk.model = np.full_like(ioblk.yData, 1.0)
    ioblk, err = ioblk.func_calcvals(ioblk)

    # physvalsavs and calcvalsavs are used to store parameters
    #  that are fixed during the calculation
    #  ***They must be populated with fixed values before moving forward
    ioblk.physvalsavs = np.copy(ioblk.physvals)
    ioblk.calcvalsavs = np.copy(ioblk.calcvals)

    ioblk.bestphysvals = np.copy(ioblk.physvals)
    ioblk.bestcalcvals = np.copy(ioblk.calcvals)
    
    ioblk.mcmc.pos = 0
    ioblk.mcmc.paridx = np.where(ioblk.fixed == 0)[0]
    ioblk.mcmc.pars = ioblk.calcvals[ioblk.mcmc.paridx]
    ioblk.mcmc.like = 0.0
    ioblk.mcmc.prior = 0.0
    ioblk.mcmc.chi2 = 0.0
    ioblk.mcmc.attempts = np.zeros_like(ioblk.mcmc.pars, dtype=np.uint64)
    ioblk.mcmc.accepts = np.zeros_like(ioblk.mcmc.pars, dtype=np.uint64)
    ioblk.mcmc.curpar = 0
    ioblk.mcmc.nfreepar = ioblk.mcmc.pars.size
    
    ioblk.chi2min = ioblk.yData.size * 2000.0
    ioblk.expchi2 = ioblk.yData.size - ioblk.mcmc.paridx.size
    ndat = ioblk.yData.size
    ioblk.likenorm = -np.sum(np.log(ioblk.errData)) - (ndat * np.log(np.sqrt(2.0*np.pi)))

    # if ioblk.scls not correct length (most likely because just doing minizer)
    #  set it to ones of correct length
    if not len(ioblk.scls) == ioblk.nparm:
        ioblk.scls = np.ones_like(ioblk.origests)

    # Setup parallel tempering things
    # Redefine temparature vector of parallel tempering turned off
    #  This way if user sets dopartemp == False but defines temperatures
    #  We won't do parallel tempering
    if not ioblk.parm.dopartemp:
        ioblk.pt.temps = np.array([1.0])
    ioblk.pt.curtempidx = 0
    ioblk.pt.curtemp = ioblk.pt.temps[0]
    ioblk.pt.ntemp = ioblk.pt.temps.size
    ioblk.pt.swapat = np.zeros((np.floor(ioblk.parm.maxstps/2).astype(int),))
    ioblk.pt.swapgd = np.copy(ioblk.pt.swapat)
    ioblk.pt.swaptemp = np.copy(ioblk.pt.swapat)
    ioblk.pt.swapattempts = np.zeros_like(ioblk.pt.temps)
    ioblk.pt.swapaccepts = np.zeros_like(ioblk.pt.temps)
    ioblk.pt.allphysval = np.tile(ioblk.physvals.reshape(ioblk.nparm,1), \
                                    (1, ioblk.pt.ntemp))
    ioblk.pt.allcalcval = np.tile(ioblk.calcvals.reshape(ioblk.nparm,1), \
                                    (1, ioblk.pt.ntemp))
    ioblk.pt.allscls = np.tile(ioblk.scls.reshape(ioblk.nparm,1), \
                                    (1, ioblk.pt.ntemp))
    ioblk.pt.alllike = np.ones_like(ioblk.pt.temps) * ioblk.mcmc.like
    ioblk.pt.allprior = np.ones_like(ioblk.pt.temps) * ioblk.mcmc.prior
    ioblk.pt.allattempts = np.zeros((ioblk.mcmc.nfreepar, ioblk.pt.ntemp))
    ioblk.pt.allaccepts = np.zeros((ioblk.mcmc.nfreepar, ioblk.pt.ntemp))

    if ioblk.parm.debugLevel > 2:
        # Setup  figures for first time
        ioblk.fighandle = plt.figure(figsize=(3,2),dpi=300,facecolor='white')
        ioblk.axhandle = plt.gca()
        ioblk.axhandle.set_position([0.125, 0.125, 0.825, 0.825])
        #ioblk.axhandle.set_axis_bgcolor('white')                            
    
    # Test prior
    ioblk, err = ioblk.func_prior(ioblk)
    if ioblk.parm.debugLevel > 0:
        print("Prior Test: ",err)
    
    # Test likelihood
    ioblk, err = ioblk.func_likehood(ioblk)
    ioblk.mcmc.like = ioblk.mcmc.like * ioblk.pt.curtemp
    if ioblk.parm.debugLevel > 0:
        print("Likelihood Test: ",err)
    ioblk.mcmc.pos = 1 # Suppress plotting temporarily
    # Test mcmc step
    ioblk = pgmcmc_one_mcmc_step(ioblk)
    
    # Refill data for all temps
    ioblk.pt.allphysval = np.tile(ioblk.physvals.reshape(ioblk.nparm,1), \
                                    (1, ioblk.pt.ntemp))
    ioblk.pt.allcalcval = np.tile(ioblk.calcvals.reshape(ioblk.nparm,1), \
                                    (1, ioblk.pt.ntemp))
    ioblk.pt.alllike = np.ones_like(ioblk.pt.temps) * ioblk.mcmc.like
    ioblk.pt.allprior = np.ones_like(ioblk.pt.temps) * ioblk.mcmc.prior
                                    
    pgmcmc_save_state(ioblk, 'ioblk_runstart')

    return ioblk

class pgmcmc_parameters:
    """Storage class for the parameters of pcmcmc routines
       CONTENTS: 
         likehoodmoddisplay - [Int] If debugLevel > =3 display likelihood call
                              model and residual every iteration mod of
                              this parameter
         debugLevel - [Int] 
         maxstps - [Int] Maximum number of mcmc steps until stopping
         avgpt - [Int] Average time between Parallel Temperature swaps
    """
    def __init__(self):
        self.dopartemp = False
        self.likehoodmoddisplay = 200
        self.debugLevel = 3
        self.maxstps = 50000
        self.fracwant = 0.25
        self.initNSteps = 200
        self.coarseNSteps = 100
        self.coarseLowFrac = 0.1
        self.coarseHghFrac = 0.8
        self.refineNSteps = 400
        self.refineLowFrac = 0.15
        self.refineHghFrac = 0.42
        self.maxPropTries = 50
        self.saveNSteps = 100
        self.avgpt = 6
        self.saveAltTemp = False
        self.savedTempIdx = 0
        

    def __str__(self):
        for k in self.__dict__:
            print(k, self.__dict__[k])
        return ''

class pgmcmc_mcmc:
    """Storage class for the runtime mcmc variables
       CONTENTS:
          pos - [Int] current step position along chain
          pars - [Array] current variable values
          paridx - [index array] index array from physvals into pars
          like - [Float] current log likelihood value
          prior - [Float] current log prior value
          chi2- [Float] current chi2 goodness of fit metric
          attempts - [Int array] - track each parameters acceptance rate with
          accepts - [Int Array]  - track "
          curpar - [Int] - Current parameter stepping in
          nfreepar - [Int] - Number of parameters being fit
    """
    def __init__(self):
        self.pos = 0
        self.pars = np.array([0.0])
        self.paridx = np.array([])
        self.like = 0.0
        self.prior = 0.0
        self.chi2 = 0.0
        self.attempts = np.array([0])
        self.accepts = np.array([0])
        self.curpar = 0
        self.nfreepar = 1
        
    def __str__(self):
        for k in self.__dict__:
            print(k, self.__dict__[k])
        return ''
        
class pgmcmc_pt:
    """Storage class for the runtime mcmc variables involving 
        parallel tempering
       CONTENTS:
          temps - [Array] List of temperatures in use
          curtempidx - [Int] Current temperature index
          curtemp - [Float] Current temperature value
          nswap - [Int] Running total of temperature swaps
          swapt - [Int] Chain step that a temp swap was attempted
          swaptemp - [Int] temperature selected for swap
          swapattempts - [Int Array] Number of times each temperature had
                          a swap attempt
          swapaccepts - [Int Array] Number of succesful temperature swaps
          all... - [Array] Storage for physval, calcval, scls,
                           likelihood, prior, attempts and accepts
                           state at each temperature
    """
    def __init__(self):
        self.temps = np.array([1.0])
        self.ntemp = 1
        self.curtempidx = 0
        self.curtemp = 1.0
        self.nswaps = 0
        self.swapat = np.array([])
        self.swapgd = np.array([])
        self.swaptemp = np.array([])
        self.swapattempts = np.array([0])
        self.swapaccepts = np.array([0])
        self.allphysval = np.array([0])
        self.allcalcval = np.array([0])
        self.allscls = np.array([0])
        self.alllike = np.array([0])
        self.allprior = np.array([0])
        self.allattempts = np.array([0])
        self.allaccepts = np.array([0])
                
    def __str__(self):
        for k in self.__dict__:
            print(k, self.__dict__[k])
        return ''
        
def pgmcmc_loadtemp(ioblk):
    """Load everything saved at current temperature to working values
    """
    ii = ioblk.pt.curtempidx
    ioblk.pt.curtemp = ioblk.pt.temps[ii]

    ioblk.physvals = np.copy(ioblk.pt.allphysval[:, ii])
    ioblk.calcvals = np.copy(ioblk.pt.allcalcval[:, ii])
    ioblk.scls = np.copy(ioblk.pt.allscls[:, ii])
    ioblk.mcmc.attempts = np.copy(ioblk.pt.allattempts[:, ii])
    ioblk.mcmc.accepts = np.copy(ioblk.pt.allaccepts[:, ii])
    ioblk.mcmc.like = np.copy(ioblk.pt.alllike[ii])
    ioblk.mcmc.prior = np.copy(ioblk.pt.allprior[ii])
    return ioblk

def pgmcmc_savetemp(ioblk):
    """ Save working values into storage for current temperature
    """
    ii = ioblk.pt.curtempidx
    ioblk.pt.allphysval[:,ii] = np.copy(ioblk.physvals)
    ioblk.pt.allcalcval[:,ii] = np.copy(ioblk.calcvals)
    ioblk.pt.allscls[:,ii] = np.copy(ioblk.scls)
    ioblk.pt.allattempts[:,ii] = np.copy(ioblk.mcmc.attempts)
    ioblk.pt.allaccepts[:,ii] = np.copy(ioblk.mcmc.accepts)
    ioblk.pt.alllike[ii] = np.copy(ioblk.mcmc.like)
    ioblk.pt.allprior[ii] = np.copy(ioblk.mcmc.prior)
    return ioblk            

def on_key_event(event):
    '''Keyboard interaction
    http://central.scipy.org/item/84/1/simple-interactive-matplotlib-plots
    courtesy thomas haslwanter'''

    #print('you pressed %s'%event.key)        
    key = event.key

    # In Python 2.x, the key gets indicated as "alt+[key]"
    # Bypass this bug:
    if key.find('alt') == 0:
        key = key.split('+')[1]

    curAxis = plt.gca()
    if key in 'aeiou':
        curAxis.set_title('Well done!')
        plt.pause(0.01)
        plt.close()
    else:
        curAxis.set_title(key + ' is not a vowel: try again to find a vowel ....')
        plt.draw()

def press_key_to_close_figure():
    '''This works in concert with on_key_event to open a figure
    that can be closed and calculation proceed with selection
    of a vowel key'''
    # Define colors, font sizes, line widths, and marker sizes
    myblack = tuple(np.array([0.0, 0.0, 0.0]) / 255.0)
    mynearblack = tuple(np.array([75.0, 75.0, 75.0]) / 255.0)
    myblue = tuple(np.array([0.0, 109.0, 219.0]) / 255.0)
    myred = tuple(np.array([146.0, 0.0, 0.0]) / 255.0)
    myorange = tuple(np.array([219.0, 209.0, 0.0]) / 255.0)
    myskyblue = tuple(np.array([182.0, 219.0, 255.0]) / 255.0)
    myyellow = tuple(np.array([255.0, 255.0, 109.0]) / 255.0)
    mypink = tuple(np.array([255.0, 182.0, 119.0]) / 255.0)
    labelfontsize = 3.0
    tickfontsize = 3.0
    datalinewidth = 1.0
    plotboxlinewidth = 1.0
    markersize = 0.5
    bkgcolor = 'white'
    axiscolor = myblack
    labelcolor = myblack
    fig = plt.figure(figsize=(8,8), facecolor=bkgcolor)
    ax = plt.gca()
#    fig, ax = plt.subplots()
    fig.canvas.mpl_connect('key_press_event', on_key_event)
    
    # Disable default Matplotlib shortcut keys:
    keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
    for key in keymaps:
        plt.rcParams[key] = ''
    
    ax.set_title('Interactive Window, select a vowel on the keyboard to close window and proceed:')

    figstydict={'labelfontsize':labelfontsize, 'tickfontsize':tickfontsize, \
                'datalinewidth':datalinewidth, 'plotboxlinewidth':plotboxlinewidth, \
                'markersize':markersize, 'bkgcolor':bkgcolor, \
                'axiscolor':axiscolor, 'labelcolor':labelcolor, \
                'myblack':myblack, 'mynearblack':mynearblack, \
                'myblue':myblue, 'myred':myred, 'myorange':myorange, \
                'myskyblue':myskyblue, 'myyellow':myyellow, 'mypink':mypink}
    return fig, ax, figstydict

#### The following functions are for running the minimzer rather than MCMC    
def boundedvals(ioblk):
    """Convert parameters to bounded versions that the minimzer will use
       INPUT:
         ioblk - [class] trp_ioblk class
       OUTPUT:
         ioblk - [class]
         err - [0 ok ; 1 not ok]
    """
    err = 0 # Error flag
    maxmindelta = ioblk.calcval_maxs - ioblk.calcval_mins
    datamindelta = np.copy(ioblk.calcvals) - ioblk.calcval_mins
    ioblk.boundedvals = -np.log( maxmindelta / datamindelta - 1.0)
    if ~np.isfinite(ioblk.boundedvals).all():
        print("Bounded Vals Bad")
        print(ioblk.boundedvals)
        print(ioblk.calcvals)
        print(ioblk.physvals)
        err = 1
    return ioblk, err

def unboundedvals(ioblk):
    """Convert bounded parameter values that the minimizer uses to physvals
       INPUT:
         ioblk - [class] trp_ioblk class
       OUTPUT:
         ioblk - [class]
         err - [0 ok ; 1 not ok]
    """
    err = 0 # Error flag
    maxmindelta = ioblk.calcval_maxs - ioblk.calcval_mins
    ioblk.calcvals = ioblk.calcval_mins + \
                     (maxmindelta / (1.0 + np.exp( -np.copy(ioblk.boundedvals) )))
    #if np.sum( np.isfinite(ioblk.physvals) ) != np.size(ioblk.boundedvals) :
    if ~np.isfinite(ioblk.calcvals).all():
        print("UnBounded Vals Bad")
        print(ioblk.boundedvals)
        print(ioblk.calcvals)
        print(ioblk.physvals)
        err = 1
    return ioblk, err

def pgmcmc_run_minimizer(ioblk, nIter, startLocs=None):
    """Peform multiple iterations starting from random initial conditions
       return the best solution in a chi2 sense among the nIter iterations
    """
    bestChi2s = np.zeros(nIter)
    bestParameters = np.zeros((ioblk.physvals.size, nIter))
    gdFits = np.zeros(nIter, dtype=np.bool)
    ioblk, err = boundedvals(ioblk)

    # physvalsavs and boundedvalsavs are used to store parameters
    #  that are fixed during the calculation
    #  ***They must be populated with fixed values before moving forward

    ioblk.boundedvalsavs = np.copy(ioblk.boundedvals)
    ioblk.bestboundedvals = np.copy(ioblk.boundedvals)
    ioblk.minimized = False
    ioblk.mcmc.curpar = 0

    for i in range(nIter):
        if startLocs is None:
            ioblk.calcvals = ioblk.calcval_mins + \
                         np.random.rand(ioblk.calcvals.size) * \
                         (ioblk.calcval_maxs - ioblk.calcval_mins)
        else:
            ioblk.calcvals = startLocs[i,:]
        # Replace random starts with parameters values that are fixed
        ioblk.calcvals = np.where(ioblk.fixed == 1, ioblk.calcvalsavs, \
                                    ioblk.calcvals)
        ioblk, err = boundedvals(ioblk)
        startParameters = ioblk.boundedvals[ioblk.mcmc.paridx]
        usemethod = 'Nelder-Mead'
        #usemethod = 'Powell'
        useoptions = {'xtol': 1e-5, 'ftol': 1e-5, 'maxiter': 2000, 'maxfev': 2000}
        #usemethod = 'CG'
        #useoptions = {'gtol': 1e-5, 'maxiter': 2000}
        allOutput = opt.minimize(minimizer_likehood, startParameters, args=(ioblk,), \
                                 method=usemethod, options=useoptions)
        ioblk.boundedvals[ioblk.mcmc.paridx] = allOutput['x']
        ioblk.boundedvals = np.where(ioblk.fixed == 1, ioblk.boundedvalsavs, \
                                     ioblk.boundedvals)
        ioblk, err = unboundedvals(ioblk)
        ioblk, err = ioblk.func_physvals(ioblk)
        chi2min = allOutput['fun']
        if ioblk.parm.debugLevel > 0:
            strout = "%s %d %s %f" % ("It: ",i," Chi2: ",chi2min)
            print(strout)
            print(ioblk.physvals)
        if np.isfinite(ioblk.physvals).all():
            gdFits[i] = True
            bestChi2s[i] = chi2min
            bestParameters[:,i] = ioblk.physvals

    # Done with iterations find the best one by chi2min
    bestMaskedIdx = np.argmin(bestChi2s[gdFits])
    ioblk.chi2min = bestChi2s[gdFits][bestMaskedIdx]
    ioblk.bestphysvals = bestParameters[:,gdFits][:,bestMaskedIdx]
    ioblk.physvals = np.copy(ioblk.bestphysvals)
    ioblk, err = ioblk.func_calcvals(ioblk)
    ioblk, err = boundedvals(ioblk)
    ioblk.bestboundedvals = ioblk.boundedvals
    if ioblk.parm.debugLevel > 0:
        strout = "%s %f" % ("Overall Best Chi2 Min: ",ioblk.chi2min)
        print(strout)
        print(ioblk.physvals)
    ioblk.minimized = True
    return ioblk

def minimizer_likehood(sParms, ioblk):
    ioblk.mcmc.pos += 1
    ioblk.boundedvals[ioblk.mcmc.paridx] = np.copy(sParms)
    ioblk, err = unboundedvals(ioblk)
    ioblk, err = ioblk.func_physvals(ioblk)
    ioblk, err = ioblk.func_likehood(ioblk)
   
    return np.log10(ioblk.mcmc.chi2) 