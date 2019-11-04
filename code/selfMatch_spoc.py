#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:51:04 2018

@author: cjburke
"""

import numpy as np
import csv
import cjb_utils as cjb
import scipy.special as spec
import toidb_federate as fed
from gather_tce_fromdvxml import tce_seed
import time
import json
import cjb_utils as cjb
try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve
try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib  

from astropy import units as u
from astropy.coordinates import SkyCoord
import pickle
import os

def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version 
    version = '3.6.3' #".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content

# Do cone search around TIC
def getTICSep(tic1, tic2):
    
    # find the position of this tic
    startTime = time.time()
    ticStringList = ['{0:d}'.format(tic1),'{0:d}'.format(tic2)]

    request = {'service':'Mast.Catalogs.Filtered.Tic', \
               'params':{'columns':'*', 'filters':[{ \
                        'paramName':'ID', 'values':ticStringList}]}, \
                'format':'json', 'removenullcolumns':True}
    while True:    
        headers, outString = mastQuery(request)
        try:
            outObject = json.loads(outString)
        except:
            outObject = []
            break
        if outObject['status'] != 'EXECUTING':
                break
        if time.time() - startTime > 30:
                print('Working...')
                startTime = time.time()
        time.sleep(5)
    # Protect against missing TICs
    try:
        curRa1 = outObject['data'][0]['ra']
        curDec1 = outObject['data'][0]['dec']
        try:
            curRa2 = outObject['data'][1]['ra']
            curDec2 = outObject['data'][1]['dec']
            c1 = SkyCoord(curRa1, curDec1, frame='icrs', unit='deg')
            c2 = SkyCoord(curRa2, curDec2, frame='icrs', unit='deg')
            sep = c1.separation(c2)
            ticSep = np.log10(sep.arcsecond)
        except:
            ticSep = -99.0
    except:
        ticSep = -99.0
    return ticSep




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


def genericFed(per, epc, tryper, tryepc, trydur, trypn, trytic, tStart, tEnd):
    ts = fed.timeseries(tStart, tEnd)
    federateResult = fed.federateFunction(per, epc, ts, \
                    tryper, tryepc, trydur)
    bstpn = trypn[federateResult[0]]
    bsttic = trytic[federateResult[0]]
    bstMatch = int(federateResult[1])
    bstStat = federateResult[2]
    bstPeriodRatio = federateResult[3]
    bstPeriodRatioFlag = federateResult[4]
    bstFederateFlag = federateResult[5]
    nFed = federateResult[6]
    
    return bstpn, bsttic, bstMatch, bstStat, bstPeriodRatio, bstPeriodRatioFlag, bstFederateFlag, nFed


if __name__ == '__main__':
    fout = open('selfMatch_sector16_20191029.txt', 'w')
    dataSpan = 27.0
    # Load the tce data pickle    
    tceSeedInFile = 'sector16_20191029_tce.pkl'
    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()
    # Check to see if cadence to time mappting is available
    hasCadTimeMap = False
    if os.path.exists('cadnoVtimemap.txt'):
        dataBlock = np.genfromtxt('cadnoVtimemap.txt', dtype=['i4','f8','i4','i4'])
        cadmap = dataBlock['f0']
        timemap = dataBlock['f1']
        print('Cadence time map available')
        hasCadTimeMap = True


    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=np.int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=np.int)
    allrp = np.array([x.at_rp for x in all_tces])
    allper = np.array([x.at_period for x in all_tces])
    alldur = np.array([x.at_dur for x in all_tces])
    allepc = np.array([x.at_epochbtjd for x in all_tces])
    alltrpvalid = np.array([x.trp_valid for x in all_tces], dtype=np.int)
    alltrpdur = np.array([x.trp_dur for x in all_tces])
    alltrpepc = np.array([x.trp_epochbtjd for x in all_tces])
    alltcedur = np.array([x.pulsedur for x in all_tces])
    alltceepc = np.array([x.tce_epoch for x in all_tces])
    alltceper = np.array([x.tce_period for x in all_tces])
    alltceCadStrt = np.array([x.data_start for x in all_tces], dtype=np.int64)
    alltceCadEnd = np.array([x.data_end for x in all_tces], dtype=np.int64)
    # Go through each tce and use valid fits from dv, trpzd, tce in that order for matching
    useper = np.zeros_like(allper)
    useepc = np.zeros_like(allper)
    usedur = np.zeros_like(allper)
    idx = np.where((allatvalid == 0) & (alltrpvalid == 0) )[0]
    useper[idx] = alltceper[idx]
    useepc[idx] = alltceepc[idx]
    usedur[idx] = alltcedur[idx]
    idx = np.where((alltrpvalid == 1))[0]
    useper[idx] = alltceper[idx]
    useepc[idx] = alltrpepc[idx]
    usedur[idx] = alltrpdur[idx]
    idx = np.where((allatvalid == 1))[0]
    useper[idx] = allper[idx]
    useepc[idx] = allepc[idx]
    usedur[idx] = alldur[idx]
    # Left over variable names from QLP version
    gtTIC = alltic
    gtPer = useper
    gtEpc = useepc
    gtDur = usedur
    gtPN = allpn
    gtTOI = allpn
    gtCadStrt = alltceCadStrt
    gtCadEnd = alltceCadEnd


    uowStart = np.min(useepc)-1.0
    uowEnd = np.max(useepc) + dataSpan + 1.0

    # put an upper limit on duration
    idx = np.where(gtDur>5.0)[0]
    gtDur[idx] = 5.0
    # Left over variable names from QLP version where only PCs were evaluated
    # Evaluate them all!!
    pcTIC = gtTIC
    pcTOI = gtPN
    pcPer = gtPer
    pcEpc = gtEpc
    pcDur = gtDur
    pcCadStrt = gtCadStrt
    pcCadEnd = gtCadEnd
    
    # Use this for debugging
#    idx = np.where(pcTIC == 220432563)[0]
#    pcTIC, pcTOI, pcPer, pcEpc, pcDur = cjb.idx_filter(idx, pcTIC, pcTOI,\
#                                                       pcPer, pcEpc, pcDur)
    for i in range(len(pcTIC)):
        curTic = pcTIC[i]
        curper = pcPer[i]
        curepc = pcEpc[i]
        curToi = pcTOI[i]
        curCadStrt = pcCadStrt[i]
        curCadEnd = pcCadEnd[i]
        # Tailor unit of work to this TCE
        if hasCadTimeMap:
            minCadStrt = np.min(curCadStrt)
            maxCadStrt = np.max(curCadEnd)
            idxfnd = np.where(minCadStrt == cadmap)[0]
            uowStartUse = timemap[idxfnd]
            idxfnd = np.where(maxCadStrt == cadmap)[0]
            uowEndUse = timemap[idxfnd]
            print('Time: {0:f} {1:f}'.format(uowStartUse[0], uowEndUse[0]))
        else:
            uowStartUse = uowStart
            uowEndUse = uowEnd
        
        # Remove all signals on the current tic
        idx = np.where(np.logical_not((curTic == gtTIC)))[0]
        usePer = gtPer[idx]
        useEpc = gtEpc[idx]
        useTic = gtTIC[idx]
        useToi = gtTOI[idx]
        useDur = gtDur[idx]
        sigMatch = np.zeros((len(usePer),), dtype=np.float)
        for j in range(len(usePer)):
            sigMatch[j] = coughlin_sigmap(curper, usePer[j])
        idxSig = np.where(sigMatch > 3.0)[0]
        if len(idxSig) > 0:
            # Potential match
            tryper = usePer[idxSig]
            tryepc = useEpc[idxSig]
            trypn = np.ones_like(idxSig)
            trydur = useDur[idxSig]
            trytic = useTic[idxSig]
            bstpn, bsttic, bstMatch, bstStat, bstPeriodRatio, bstPeriodRatioFlag, bstFederateFlag, nFed = \
                    genericFed(curper, curepc, tryper, tryepc, trydur, trypn, trytic, uowStartUse, uowEndUse)
            # If federation found calculate spatial separation between TICs
            bstSep = -99.0
            if bstFederateFlag == 1:
                bstSep = getTICSep(curTic, bsttic)

            str = '{:12d} {:2d} {:12d} {:2d} {:2d} {:6.3f} {:2d} {:10.5f} {:2d} {:6.2f} {:d}\n'.format(curTic, \
                   curToi, bsttic, bstpn, bstMatch, bstStat, bstPeriodRatioFlag, \
                   bstPeriodRatio, bstFederateFlag, bstSep, nFed)
            fout.write(str)
            print(str)
        
#        for j in range(len(idxSig)):
#            k = idxSig[j]
#            print('TIC {:d} Match w/ TIC: {:d} {:7.2f}'.format(curTic, useTic[k], sigMatch[k] ))
        
#        str = '{:12d} {:8.2f} {:2s} {:12d} {:2d} {:2d} {:6.3f} {:2d} {:10.5f} {:2d}\n'.format(curTic, \
#                   curToi, gtDisp[i], bsttic, bstpn, bstMatch, bstStat, bstPeriodRatioFlag, \
#                   bstPeriodRatio, bstFederateFlag)
#        fout.write(str)
#        print(str)
    fout.close()
    print('hello world')