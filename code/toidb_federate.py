# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:20:21 2017
This is the main routines to perform detailed federation or ephemeris matching
This was implemented as the official algorithm in the Kepler science office
for matching TCEs and KOIs from different pipeline runs.  It matches using
both period and phase.  It works quite well for finding matches in data
with ephemerides determined from the same source.  However, it is too 
restrictive and not recommended for matching
between ephemerides from different data. For instance matching previously
known planets and TESS data does not work at this level of detail and
one should revert to just period and host ID matching.  In other words, 
this routine is appropriate for planet ephemerides found in TESS data
between different pipelines, but both using TESS data.

The algorithm is briefly described in Sec. 4.1 of
http://adsabs.harvard.edu/abs/2015ApJS..217...31M
@author: Christopher J. Burke
"""


import argparse
import numpy as np
import sys

import cjb_utils as cjb


class timeseries:
    tres = 0
    deltaT = 0.0
    ts = np.array([])
    nt = 0
    tsmjd = np.array([])
    tskjd = np.array([])
    tsjd = np.array([])
    
    def __init__(self, tstrt, tend):
        self.tres = 10 # time series resolution in MIN
        self.deltaT = self.tres / 1440.0 # resolution in Days
        self.ts = np.arange(tstrt, tend, self.deltaT)
        self.nt = len(self.ts)
        self.tsmjd = np.copy(self.ts)
        self.tskjd = np.copy(self.ts - 54832.5)
        self.tsjd = np.copy(self.ts + 2400000.5)
        self.ts = np.copy(self.ts) # whatever is left in self.ts is system 
                                        # used in calculation
        self.ephemCentral = np.zeros((self.nt,), dtype=np.int8)
        self.ephemFull = np.zeros_like(self.ephemCentral)
        
    def makeEphemVector(self, period, epoch, duration):
        """ Make the ephemeris arrays """
        self.ephemCentral = np.zeros((self.nt,), dtype=np.int8)
        self.ephemFull = np.zeros_like(self.ephemCentral)
        phase = np.mod((self.ts - epoch), period) / period
        phaseDuration = (duration / 24.0) / period / 2.0
        idx = np.where(phase > 0.5)[0]
        phase[idx] = phase[idx] - 1.0
        idx = np.where(np.abs(phase) <= phaseDuration)[0]
        self.ephemFull[idx] = 1
        # Detect edges of sequences of ones
        offsetdiff = np.diff(self.ephemFull)
        offsetdiff = np.insert(offsetdiff, 0, self.ephemFull[0])
        offsetdiff = np.insert(offsetdiff, len(offsetdiff), self.ephemFull[-1])
        # Get indices that edges occur on 
        idx = np.where(np.abs(offsetdiff) == 1)[0]
        # Get mean of the edge indices in order to find central index
        b = np.reshape(idx, (int((len(idx)/2)), 2))
        c = np.zeros_like(b)
        c[:,1] = 1
        idx = np.int64((np.mean((b-c), 1)).flatten())
        self.ephemCentral[idx] = 1
        return self
        

def federateFunction(period, epoch, ts, tcePeriods, tceEpochs, tceDurations):
    # Here are some algorithm parameters
    # Only consider periods within minmult and maxmult factors
    minMult = 1.0/3.0
    maxMult = 3.0
    # minimum period to attempt a match [day]
    minPeriod = 0.11
    # 1:1 period matching tolerance
    periodMatchTol = 0.006
    # duration minimum in hours
    durFloor = 1.0
    durFloorShort = 0.5 # If period < 1 day use smaller durFloor
    
    minPer = minMult * period
    maxPer = maxMult * period
    curntce = len(tcePeriods)
    # Make ephemeris vector for TOI
    ts = ts.makeEphemVector(period, epoch, 1.0)
    toiephem = ts.ephemCentral
    toiephem2d = np.tile(toiephem, (curntce, 1))
    nindatabase = np.sum(toiephem2d, 1)
    idx = np.where(nindatabase == 0)[0]
    nindatabase[idx] = 1
    # Build the 2D ephem matrix for the TCEs
    curngd = np.zeros((curntce,))
    curngd2 = np.zeros_like(curngd)
    allpers = np.zeros_like(curngd)
    matchres = np.zeros_like(curngd)
    tceephem2d = np.zeros_like(toiephem2d)
    for kk in range(curntce):
        tcePeriod = tcePeriods[kk]
        allpers[kk] = tcePeriod
        tceEpoch = tceEpochs[kk]
        tceDuration = tceDurations[kk]
        tceDuration = tceDuration * 2.0
        if (tceDuration/24.0 / tcePeriod > 0.25):
            tceDuration = tcePeriod * 0.25 * 24.0
        if tcePeriod < 1.0:
            durFloor = durFloorShort
        tceDuration = np.max([tceDuration, durFloor])
        ts = ts.makeEphemVector(tcePeriod, tceEpoch, tceDuration)
        tceephem2d[kk,:] = ts.ephemFull
        curngd2[kk] = np.max([1.0,np.sum(ts.ephemCentral)])
        curngd[kk] = np.sum(ts.ephemFull)
    inverttceephem2d = np.logical_not(tceephem2d)
    corr1 = tceephem2d * toiephem2d
    corr2 = inverttceephem2d * toiephem2d
    stats = np.sum(corr1,1) / curngd2 - np.sum(corr2,1) / nindatabase
    idx = np.where((allpers < minPer) | (allpers > maxPer))[0]
    stats[idx] = -1.0
    if (period < minPeriod):
        stats = -np.ones_like(stats)
    idx = np.where(stats >= 0.8)[0]
    matchres[idx] = 1
    idx = np.where((stats >=0.4) & (stats < 0.8))[0]
    matchres[idx] = 2
    idx = np.where((stats >=0.15) & (stats < 0.4))[0]
    matchres[idx] = 3
    idx = np.where((stats >=-0.5) & (stats < 0.15) & (np.abs(1.0 - \
                            period / allpers) < periodMatchTol))[0]
    matchres[idx] = 4
    idxsort = np.argsort(stats)[::-1]
    sortedstats = stats[idxsort]
    bstIdx = idxsort[0]
    bstMatch = matchres[bstIdx]
    bstStat = sortedstats[0]
    bstPeriodRatio = period / allpers[bstIdx]
    bstPeriodRatioFlag = 0
    bstFederate = 0
    if (abs(1.0 - bstPeriodRatio) < periodMatchTol):
        bstPeriodRatioFlag = 1
    # overall federation condition
    if (bstPeriodRatioFlag == 1):
        if (bstMatch >= 1) and (bstMatch <= 4):
            bstFederate = 1
    else:
        if (bstMatch >=1) and (bstMatch <=2) and (bstPeriodRatio>0.45):
            bstFederate = 1

    # Also return the number of inputs that have a federation
    allPeriodRatio = period / allpers
    allPeriodRatioFlag = np.zeros_like(allPeriodRatio, dtype=np.int)
    idx = np.where(np.abs(1.0 - allPeriodRatio) < periodMatchTol)[0]
    allPeriodRatioFlag[idx]
    allFederate = np.zeros_like(allPeriodRatio, dtype=np.int)
    idx = np.where((allPeriodRatioFlag == 1) & (matchres>=1) & (matchres<=4))[0]
    allFederate[idx] = 1
    idx = np.where((allPeriodRatioFlag == 0) & (matchres>=1) & (matchres<=2) & (allPeriodRatio>0.45))[0]
    allFederate[idx] = 1
    nFed = len(np.where(allFederate == 1)[0])
    return bstIdx, bstMatch, bstStat, bstPeriodRatio, bstPeriodRatioFlag, bstFederate, nFed
    

def federate2TOI(tce_data, db, sector, quiet=False, debugFigures=False):
    """ Copy Kepler federation code to perform similar federation with TESS
    """
    
    # get the tce data into more convenient np arrays
    tceTics = np.array([obj.ticNumber for obj in tce_data])
    tcePeriods = np.array([obj.period for obj in tce_data])
    tceEpochs = np.array([obj.epoch for obj in tce_data])
    tceDurations = np.array([obj.duration for obj in tce_data])
    tcePlanNs = np.array([obj.planetNumberCreate for obj in tce_data])
    
    # Get the latest toi information from DB
    tois = np.array([dict_['toiNumber'] for dict_ in db])
    tics = np.array([dict_['ticNumber'] for dict_ in db])
    epochs = np.array([dict_['epoch'] for dict_ in db])
    periods = np.array([dict_['period'] for dict_ in db])
    tes = np.array([dict_['timeEnd'] for dict_ in db])
    # Only get entries with latest data
    idx = np.where(tes == 0)[0]
    tois, tics, epochs, periods = cjb.idx_filter(idx, tois, tics, epochs, periods)
    # Sort the entries in the off chance that the order of data returned from db is 
    #  random
    idx = np.argsort(tois)
    tois, tics, epochs, periods = cjb.idx_filter(idx, tois, tics, epochs, periods)
    nToi = len(tois)
    
    # Setup some output variables in the TOI space
    toi2Tce = np.zeros((nToi,))
    toiToi = np.zeros_like(toi2Tce)
    toiTic = np.zeros((nToi,), dtype=np.int64)
    toiPlanN = np.zeros((nToi,), dtype=np.int)
    toiMatchQual = np.zeros((nToi,), dtype=np.int)
    toiPeriod = np.zeros_like(toi2Tce)
    toiTcePeriodRatio = np.zeros_like(toi2Tce)
    toiTcePeriodRatioFlag = np.zeros_like(toiMatchQual)
    toiMatchStat = np.zeros_like(toiMatchQual)
    toiFederate = np.zeros_like(toiMatchQual)
    
    # Setup  ephemeris time series
    # The sector is used to get the start and end times from toidb_addtoi.py 
    #  unitofworkdata class
    uow = toiadd.unitOfWorkData()
    idx = np.where(uow.sectorNumbers == sector)[0]
    tStart = uow.sectorStarts[idx]
    tEnd = uow.sectorEnds[idx]
    ts = timeseries(tStart, tEnd)
    
    # Now go through TOIs and try to match them
    for i in range(len(tois)):
        curToi = tois[i]
        toiToi[i] = tois[i]
        curTic = tics[i]
        toiTic[i] = curTic
        curEpoch = epochs[i]
        curPeriod = periods[i]
        toiPeriod[i] = curPeriod
        # Try to find this target in TCEs        
        idxintce = np.where(tceTics == curTic)[0]
        curntce = len(idxintce)
        if (curntce > 0):
            if not quiet:
                print("Start match {0:f} TIC: {1:d} NTCEs: {2:d}".format(curToi, curTic, curntce))
            federateResult = federateFunction(curPeriod, curEpoch, ts, \
                    tcePeriods[idxintce], tceEpochs[idxintce], tceDurations[idxintce])
            toiPlanN[i] = tcePlanNs[idxintce][federateResult[0]]
            toiMatchQual[i] = federateResult[1]
            toiMatchStat[i] = federateResult[2]
            toiTcePeriodRatioFlag[i] = federateResult[4]
            toiTcePeriodRatio[i] = federateResult[3]
            toiFederate[i] = federateResult[5]
            print("{0:d} {1:f} {2:d} Matches {3:d} Prat: {4:f}".format( i, curToi, \
                        curTic, idxintce[federateResult[0]], toiTcePeriodRatio[i]))
    return toiToi, toiTic, toiPlanN, toiMatchQual, toiMatchStat, toiTcePeriodRatioFlag, \
                toiTcePeriodRatio, toiFederate
                
if  __name__ == "__main__":
    
    # Parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("sector", type=int,
                        help="TESS Sector Number (integer)", choices=[1])
    parser.add_argument("inputFile",
                        help="File with entries to federate to TOI list")
    parser.add_argument("outputFile",
                        help="Output file describing federation result for lines in inputFile")
    parser.add_argument("-databaseFile",
                        help="Choose an alternative database file (default - tst_db.json)", \
                        default='tst.db')
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress messages")
    args = parser.parse_args()


    #connect to database file
    if not args.quiet:
        print("Using {:s} for database".format(args.databaseFile))
    db = TinyDB(args.databaseFile)
    
    # Validate  input file
    if not args.quiet:
        print("Validating {:s} input file".format(args.inputFile))
    in_data, valid = toiadd.read_validate_toi_file(args.inputFile)
    if not valid:
        print("Input file failed to validate")
        sys.exit(1)
    if not args.quiet:
        print("Validating complete")
    
    # Determine federation 
    if not args.quiet:
        print("Start federating input against TOI database")
    federateResults = federate2TOI(in_data, db, args.sector)
    if not args.quiet:
        print("End federation against database")

        
    if not args.quiet:
        print("Writing output file of federation")
    toiToi = federateResults[0]
    toiTic = federateResults[1]
    toiPlanN = federateResults[2]
    toiMatchQual = federateResults[3]
    toiMatchStat = federateResults[4]
    toiTcePeriodRatioFlag = federateResults[5]
    toiTcePeriodRatio = federateResults[6]
    toiFederate = federateResults[7]

    fid = open(args.outputFile, 'w')
    for i in range(len(toiToi)):
        fid.write("{0:7.2f} {1:9d} {2:2d} {3:2d} {4:6.3f} {5:2d} {6:10.5f} {7:2d}\n".format( \
                    toiToi[i], toiTic[i], toiPlanN[i], toiMatchQual[i], toiMatchStat[i], \
                    toiTcePeriodRatioFlag[i], toiTcePeriodRatio[i], toiFederate[i]))
    
        
    