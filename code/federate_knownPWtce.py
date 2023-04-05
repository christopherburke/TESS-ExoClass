#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:51:04 2018
Federate/ephemeris match known confirmed planets hosted
at NExScI with the TCEs.  Must be online for the search to
work.  Also queries MAST for TIC's that may be near.

@author: Christopher J. Burke
"""

import numpy as np
import toidb_federate as fed
from gather_tce_fromdvxml import tce_seed
import scipy.special as spec
import csv
import sys
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
import urllib
# Removign python 2 support
#import urllib2
import os
import ssl

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

# Do cone search around TIC
def query_ticcone(curRa, curDec, searchRad):
    

    # Do cone search around this position
    startTime = time.time()
    request = {'service':'Mast.Catalogs.Filtered.Tic.Position', \
               'params':{'columns':'c.*', \
                         'filters':[ \
                                    {'paramName':'Tmag',\
                                     'values':[{'min':0, 'max':20.0}]}], \
                         'ra':'{:10.5f}'.format(curRa),\
                         'dec':'{:10.5f}'.format(curDec),\
                         'radius':'{:10.7f}'.format(searchRad/3600.0) \
                         }, \
                'format':'json', 'removenullcolumns':False}
    while True:    
        headers, outString = mastQuery(request)
        try:
            outObject = json.loads(outString)
        except:
            print("Error in json load")
            print("Request")
            print(request)
            print("outString")
            print(outString)
        if outObject['status'] != 'EXECUTING':
                break
        if time.time() - startTime > 30:
                print('Working...')
                startTime = time.time()
        time.sleep(5)
     
    ticList = [x['ID'] for x in outObject['data']]


    return ticList

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


if __name__ == '__main__':
    fout = open('federate_knownP_sector-62_20230404.txt', 'w')
    dataSpan = 27.0
    wideSearch = True
    searchRad = 180.0 # Arcsecond search radius for other TICs
    # Check to see if cadence to time mappting is available
    hasCadTimeMap = False
    if os.path.exists('cadnoVtimemap.txt'):
        dataBlock = np.genfromtxt('cadnoVtimemap.txt', dtype=['i4','f8','i4','i4','i4'])
        cadmap = dataBlock['f0']
        timemap = dataBlock['f1']
        print('Cadence time map available')
        hasCadTimeMap = True


    # Load known transiting planet table from NEXSCI   
    #whereString = 'pl_tranflag = 1 and st_elat>0.0'
    selectString = 'select pl_name,pl_orbper,pl_tranmid,pl_trandur,ra,dec,elat'
    fromString = 'from pscomppars'
    whereString = 'where tran_flag = 1'
    url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?'
    data = {'query':'{0} {1} {2}'.format(selectString, fromString, whereString), \
                    'format':'csv'}

    url_values = urllib.parse.urlencode(data)
    print(url_values)
    
    queryData = urllib.request.urlopen(url + url_values)
 

#    url_values = urllib.parse.urlencode(data)
    #print url_values
#    queryData = urllib.request.urlopen(url + url_values)
    returnPage = queryData.read()
    dtypeseq = ['U40']
    dtypeseq.extend(['f8'] * 6)
    dataBlock = np.genfromtxt(returnPage.splitlines(), delimiter=',', skip_header=1, \
                        dtype=dtypeseq)
    gtName = dataBlock['f0']
    gtPer = dataBlock['f1']
    gtEpc = dataBlock['f2']-2457000.0
    gtDur = dataBlock['f3'] # Used to be in days now its in hours*24.0 # Is in days convert to hours
    gtRa = dataBlock['f4']
    gtDec = dataBlock['f5']
    gtEclipLat = dataBlock['f6']
    gtTIC = np.arange(len(gtName))
    gtTOI = np.arange(len(gtName))
    
    # Filter for planets in the correct ecliptic area to speed this up
    # Ecliptic pointing
    #idx = np.where((gtEclipLat > -20.0) & (gtEclipLat < 20.0))[0]
    # North Ecliptic pointing
    #idx = np.where((gtEclipLat > 3.0))[0]
    # South Ecliptic pointing
    idx = np.where((gtEclipLat < -3.0))[0]
    gtName = gtName[idx]
    gtPer, gtEpc, gtDur, gtRa, gtDec, gtTIC, gtTOI = cjb.idx_filter(idx, \
        gtPer, gtEpc, gtDur, gtRa, gtDec, gtTIC, gtTOI                                                            )

    
    # Check for missing ephemeris values
    idxBd = np.where((np.logical_not(np.isfinite(gtPer))) | \
                     (np.logical_not(np.isfinite(gtEpc))))[0]
#    idxBd = np.where((np.logical_not(np.isfinite(gtPer))))[0]
    if len(idxBd) > 0:
        for curIdxBd in idxBd:
            print('Bad Ephemeris for Known Planet')
            #print(gtName[curIdxBd])
            # Load known planet data from NEXSCI
            # Using TAP interface
            selectString = 'select pl_name,pl_orbper,pl_tranmid,pl_trandur'
            fromString = 'from ps'
            #whereString = 'where pl_name like \'%{0}%\''.format(gtName[curIdxBd])
            curName = gtName[curIdxBd].replace('"','')
            print(curName)
            whereString = 'where pl_name = \'{0}\''.format(curName)
            url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?'
            data = {'query':'{0} {1} {2}'.format(selectString, fromString, whereString), \
                    'format':'csv'}
            url_values = urllib.parse.urlencode(data)
            #print(url_values)
            queryData = urllib.request.urlopen(url + url_values)
            returnPage = queryData.read()
            dtypeseq = ['U40']
            dtypeseq.extend(['f8'] * 3)
            dataBlock = np.genfromtxt(returnPage.splitlines(), delimiter=',', skip_header=1, \
                        dtype=dtypeseq)
            curName = dataBlock['f0']
            curPer = dataBlock['f1']
            curEpc = dataBlock['f2']
            curDur = dataBlock['f3']
            if curPer.size > 1:
                gotValues = False
                idxGd = np.where(np.isfinite(curPer))[0]
                for ip in range(len(curPer)):
                    if np.isfinite(curPer[ip]) and np.isfinite(curEpc[ip]) and not gotValues:
                        gotValues = True
                        gtPer[curIdxBd] = curPer[ip]
                        gtEpc[curIdxBd] = curEpc[ip]
                        print('Found {0} {1} {2}'.format(gtPer[curIdxBd], gtEpc[curIdxBd], np.mean(curPer[idxGd])))
            else:
                print('Not Matching to {0}'.format(gtName[curIdxBd]))
                print(curPer, curEpc)
    # This removes planets with no valid period found
    idx = np.where((np.isfinite(gtPer)))[0]
    gtName = gtName[idx]
    gtPer, gtEpc, gtDur, gtRa, gtDec, gtTIC, gtTOI = cjb.idx_filter(idx, \
        gtPer, gtEpc, gtDur, gtRa, gtDec, gtTIC, gtTOI                                                            )
    print("Tot # PCs: {0:d}".format(len(gtName)))


    # Load the tce data h5
    tceSeedInFile = 'sector-62_20230404_tce.h5'
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=int)
    allrp = np.array([x.at_rp for x in all_tces])
    allper = np.array([x.at_period for x in all_tces])
    alldur = np.array([x.at_dur for x in all_tces])
    allepc = np.array([x.at_epochbtjd for x in all_tces])
    alltrpvalid = np.array([x.trp_valid for x in all_tces], dtype=int)
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
#    ia = np.argsort(alltic)
#    alltic, allpn, useper, useepc, usedur = cjb.idx_filter(ia, alltic, \
#                                allpn, useper, useepc, usedur)

    # write to header the inputs
    fout.write('# Match {:s}\n'.format('Exoplanet Archive Confirmed Planets'))
    fout.write('# To {:s}\n'.format(tceSeedInFile))

#    if np.min(useepc)-1.0 < uowStart:
    uowStart = np.min(useepc)-1.0
    uowEnd = np.max(useepc) + dataSpan + 1.0
    # Go  the ground truth data (ground truth acts like KOIs in kepler federation)
    nGT = len(gtTIC)
    for i in range(len(gtTIC)):
        curTic = gtTIC[i]
        curper = gtPer[i]
        curepc = gtEpc[i]
        curToi = gtTOI[i]
        curName = gtName[i]
        curRa = gtRa[i]
        curDec = gtDec[i]
        if np.mod(i,20) == 0:
            print('Done {0:d} of {1:d}'.format(i, nGT))
        # If wideSearch True then query MAST for 
        #  all TICs within searchRad arcsec of this target
        if wideSearch:
            otherTICs = np.array(query_ticcone(curRa, curDec, searchRad), dtype=np.int64)
            #otherTICs = np.sort(otherTICs)
            idx = np.array([], dtype=np.int64)
            for ot in otherTICs:
                idxtmp = np.where(ot == alltic)[0]
                idx = np.append(idx, idxtmp)                
        else:
            # find this tic in the tces
            idx = np.where(curTic == alltic)[0]
        if len(idx) > 0:
            # Potential match
            tryper = useper[idx]
            tryepc = useepc[idx]
            trypn = allpn[idx]
            trydur = usedur[idx]
            trytic = alltic[idx]
            tryCadStrt = alltceCadStrt[idx]
            tryCadEnd = alltceCadEnd[idx]
            # Tailor unit of work to longest duration among potential TCEs
            if hasCadTimeMap:
                minCadStrt = np.min(tryCadStrt)
                maxCadStrt = np.max(tryCadEnd)
                idxfnd = np.argmin(np.abs(minCadStrt-cadmap))
                #idxfnd = np.where(minCadStrt == cadmap)[0]
                uowStartUse = timemap[idxfnd]
                idxfnd = np.argmin(np.abs(maxCadStrt - cadmap))
                #idxfnd = np.where(maxCadStrt == cadmap)[0]
                uowEndUse = timemap[idxfnd]
                print('Time: {0:f} {1:f}'.format(uowStartUse, uowEndUse))
            else:
                uowStartUse = uowStart
                uowEndUse = uowEnd
                
            bstpn, bsttic, bstMatch, bstStat, bstPeriodRatio, bstPeriodRatioFlag, bstFederateFlag, nFed = \
                    genericFed(curper, curepc, tryper, tryepc, trydur, trypn, trytic, uowStartUse, uowEndUse)

#            sigMatch = np.zeros((len(tryper),), dtype=float)
#            for j in range(len(tryper)):
#                sigMatch[j] = coughlin_sigmap(curper, tryper[j])
#            idxSig = np.where(sigMatch > 3.0)[0]
#            if len(idxSig)>0:
#                tryper = tryper[idxSig]
#                tryepc = tryepc[idxSig]
#                trypn = trypn[idxSig]
#                trydur = trydur[idxSig]
#                trytic = trytic[idxSig]
#                idxSigia = np.argsort(trypn)[0]
#                bstpn = trypn[idxSigia]
#                bsttic = trytic[idxSigia]
#                bstStat = sigMatch[idxSigia]
#                bstPeriodRatio = curper / tryper[idxSigia]
#                bstPeriodRatioFlag = 1
#                bstFederateFlag = 1
#                bstMatch = 1
#            else:
#                bstpn = 0
#                bstMatch = -1
#                bstStat = -1.0
#                bstPeriodRatio = -1.0
#                bstPeriodRatioFlag = 0
#                bstFederateFlag = 0
#                bsttic = 0
                
            
        else:
            #print("No match in Ground truth")
            bstpn = 0
            bstMatch = -1
            bstStat = -1.0
            bstPeriodRatio = -1.0
            bstPeriodRatioFlag = 0
            bstFederateFlag = 0
            bsttic = 0
        str = '{:12d} {:8.2f} {:s} {:12d} {:2d} {:2d} {:6.3f} {:2d} {:10.5f} {:2d}\n'.format(curTic, \
                   curToi, curName.astype('str').replace(" ",""), bsttic, bstpn, bstMatch, bstStat, bstPeriodRatioFlag, \
                   bstPeriodRatio, bstFederateFlag)
        if bstpn > 0:
            fout.write(str)
            print(str)
    fout.close()
    print('hello world')
