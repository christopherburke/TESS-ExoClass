"""  Module for Pretty Good MCMC parameter estimation
    Heavily based upon Philip Gregory's awesome Bayesian book
    Author: Christopher J Burke 
    These are routines to make diagnostic plots 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import h5py
import pickle
from pgmcmc import pgmcmc_ioblk, pgmcmc_parameters, pgmcmc_mcmc
from pgmcmc import press_key_to_close_figure, on_key_event

def acf(x, lag=20):
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1] \
        for i in range(1,lag)])
            
def pgmcmc_load_state(prefix):
    f = h5py.File(prefix+'.hd5','r')
    pvals = np.array(f['Pvals'])
    cvals = np.array(f['Cvals'])
    bvals = np.array(f['Bvals'])
    
    ioblk = pickle.load(open(prefix+'.pkl', 'rb'))
    return ioblk, pvals, cvals, bvals

def parameter_diagnostics(nameString, data, likes, priors, chi2s, bestVal):
    percents=np.array([0.5-0.9973/2.0, 0.5-0.9545/2.0, 0.5-0.6827/2.0, \
                       0.5, 0.5+0.6827/2.0, 0.5+0.9545/2.0, 0.5+0.9973/2.0]) \
                     * 100.0
    limits=np.array([1.0-0.9973, 1.0-0.9545, 1.0-0.6827, \
                     0.6827, 0.9545, 0.9973]) * 100.0

    prcresult = np.percentile(data, percents)
    limresult = np.percentile(data, limits)

    # Plot for raw data series    
    fig, ax, fsd = press_key_to_close_figure()
    # Turn off default axis
    ax.set_visible(False)
    gs = gridspec.GridSpec(5,5)
    ax1 = fig.add_subplot(gs[3,:])
    ax1.plot(data,'.')
    #ax1.set_xlabel('Chain Step', size=fsd['labelfontsize'])
    ax1.set_ylabel(nameString, size=fsd['labelfontsize'])
    ax1.tick_params('both', labelsize=fsd['tickfontsize'],
                           width=fsd['plotboxlinewidth']/3.0, 
                           color=fsd['axiscolor'], length=fsd['plotboxlinewidth'])


    # Histogram plot
    nBins = 35
    xmin = limresult[0]
    xmax = limresult[-1]
    xedges = np.linspace(xmin, xmax, nBins+1)
    #midx = xedges[:-1] + np.diff(xedges)/2.0
    ax2 = fig.add_subplot(gs[0:3,:])    
    n, xedgesused, patches = ax2.hist(data, xedges,
                                          histtype='bar',
                                          edgecolor='k', color=fsd['myskyblue'],
                                          normed=True)
    ax2.axvline(bestVal)
    ax2.set_title(nameString, size=fsd['labelfontsize'])
    ax2.set_ylabel('Prob dbin', size=fsd['labelfontsize'])
    ax2.tick_params('both', labelsize=fsd['tickfontsize'],
                           width=fsd['plotboxlinewidth']/3.0, 
                           color=fsd['axiscolor'], length=fsd['plotboxlinewidth'])

    # Autocorrelation
    acfdata = acf(data, 30)
    ax3 = fig.add_subplot(gs[4,:])
    ax3.plot(acfdata,'-k')
    ax3.set_xlabel('Lag', size=fsd['labelfontsize'])
    ax3.set_ylabel('ACF', size=fsd['labelfontsize'])
    ax3.tick_params('both', labelsize=fsd['tickfontsize'],
                           width=fsd['plotboxlinewidth']/3.0, 
                           color=fsd['axiscolor'], length=fsd['plotboxlinewidth'])

    plt.show()

    return

def parameter_correlations(data, names):
    sz = data.shape
    
    for i in np.arange(sz[1]):
        for j in np.arange(i+1,sz[1]):
            xData = data[:,i]
            xString = names[i]
            yData = data[:,j]
            yString = names[j]
            fig, ax, fsd = press_key_to_close_figure()
            ax.hexbin(xData, yData, gridsize=50, bins='log', \
                      mincnt = 3, cmap=plt.cm.YlOrRd_r)
            ax.set_xlabel(xString, size=fsd['labelfontsize'])
            ax.set_ylabel(yString, size=fsd['labelfontsize'])
            ax.tick_params('both', labelsize=fsd['tickfontsize'],
                           width=fsd['plotboxlinewidth']/3.0, 
                           color=fsd['axiscolor'], length=fsd['plotboxlinewidth'])
            plt.show()
    
    return
    
def mcmc_stats(ioblk):
    # Show proposal step acceptances
    nfreepar = ioblk.mcmc.paridx.size
    print('Acceptance Fractions')
    for i in range(nfreepar):
        j = ioblk.mcmc.paridx[i]
        curname = ioblk.physval_names[j]
        curaccept = ioblk.mcmc.accepts[i]
        curattempt = ioblk.mcmc.attempts[i]
        print("Parameter: {0:s} {1:f}".format( \
                 curname, curaccept/curattempt))
                 
    if ioblk.parm.dopartemp:
        for t in range(ioblk.pt.ntemp):
            print("Accept Fraction By Temperature: {0:d}".format(t))
            for i in range(nfreepar):
                j = ioblk.mcmc.paridx[i]
                curname = ioblk.physval_names[j]
                curaccept = ioblk.pt.allaccepts[i,t]
                curattempt = ioblk.pt.allattempts[i,t]
                print ("Parameter: {0:s} {1:f}".format( \
                     curname, curaccept/curattempt))
        for t in range(ioblk.pt.ntemp-1):
            curattempt = ioblk.pt.swapattempts[t]
            curaccept = ioblk.pt.swapaccepts[t]
            print( "Temp: {0:d} Attmpts: {1:f} Accpt Fraction: {2:f}".format( \
                    t, curattempt, curaccept/curattempt))
                
    
if __name__ == "__main__":
    
    ioblk, pvals, cvals, bvals  = pgmcmc_load_state('run')
    mcmc_stats(ioblk)
    nfreepar = ioblk.mcmc.paridx.size
    Clipout = 200
    currentPos = ioblk.mcmc.pos
    bestPhysVals = ioblk.bestphysvals
    bestCalcVals = ioblk.bestcalcvals
    for i in range(nfreepar):
        j = ioblk.mcmc.paridx[i]
        currentName = ioblk.physval_names[j]
        currentBest = bestPhysVals[j]
        currentData = pvals[Clipout:currentPos-1,i]
        currentLikes = bvals[Clipout:currentPos-1,0]
        currentPriors = bvals[Clipout:currentPos-1,1]
        currentChi2s = bvals[Clipout:currentPos-1,2]
        parameter_diagnostics(currentName, currentData, \
                 currentLikes, currentPriors, currentChi2s, currentBest)
                 
    parameter_correlations(pvals[Clipout:currentPos-1,:], ioblk.physval_names[ioblk.mcmc.paridx])
    
    for i in range(nfreepar):
        j = ioblk.mcmc.paridx[i]
        currentName = ioblk.calcval_names[j]
        currentBest = bestCalcVals[j]
        currentData = cvals[Clipout:currentPos-1,i]
        currentLikes = bvals[Clipout:currentPos-1,0]
        currentPriors = bvals[Clipout:currentPos-1,1]
        currentChi2s = bvals[Clipout:currentPos-1,2]
        parameter_diagnostics(currentName, currentData, \
                 currentLikes, currentPriors, currentChi2s, currentBest)
                 
    parameter_correlations(cvals[Clipout:currentPos-1,:], ioblk.calcval_names[ioblk.mcmc.paridx])
                 
    