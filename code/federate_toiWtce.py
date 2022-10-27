#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:51:04 2018

@author: cjburke
"""

import numpy as np
import toidb_federate as fed
from gather_tce_fromdvxml import tce_seed
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
import os
import astropy.units as u
from astropy.coordinates import SkyCoord


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
def query_othertics(ticWant, searchRad):
    
    # find the position of this tic
    startTime = time.time()
    request = {'service':'Mast.Catalogs.Filtered.Tic', \
               'params':{'columns':'*', 'filters':[{ \
                        'paramName':'ID', 'values':['{:d}'.format(ticWant)]}]}, \
                'format':'json', 'removenullcolumns':True}
    while True:    
        headers, outString = mastQuery(request)
        outObject = json.loads(outString)
        if outObject['status'] != 'EXECUTING':
                break
        if time.time() - startTime > 30:
                print('Working...')
                startTime = time.time()
        time.sleep(5)
    # Protect against missing TICs
    try:
        curRa = outObject['data'][0]['ra']
        curDec = outObject['data'][0]['dec']

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
            outObject = json.loads(outString)
            if outObject['status'] != 'EXECUTING':
                    break
            if time.time() - startTime > 30:
                    print('Working...')
                    startTime = time.time()
            time.sleep(5)
         
        ticList = [x['ID'] for x in outObject['data']]
    except:
        # No tic or somehow fail
        print('TIC Search {:d} Fail!'.format(ticWant))
        ticList = [ticWant]

    return ticList



if __name__ == '__main__':
    fout = open('federate_toiWtce_sector-56_20221023.txt', 'w')
    dataSpan = 27.0

    wideSearch = True # Do MASTTIC query if true to search
                        # for nearby matches in case the
                        # signal is found on a different TIC
                        # than in the TOIcatalog
    searchRad = 180.0 # Arcsecond search radius for other TICs
    # Check to see if cadence to time mappting is available
    hasCadTimeMap = False
    if os.path.exists('cadnoVtimemap.txt'):
        dataBlock = np.genfromtxt('cadnoVtimemap.txt', dtype=['i4','f8','i4','i4','i4'])
        cadmap = dataBlock['f0']
        timemap = dataBlock['f1']
        print('Cadence time map available')
        hasCadTimeMap = True
    
    # load the TOI data
    #MAST Format
#    qlpfile = 'hlsp_tess-data-alerts_tess_phot_alert-summary-s01+s02+s03+s04_tess_v9_spoc.csv'
#    dtypeseq = ['i4','f8','U2']
#    dtypeseq.extend(['f8']*10)
#    dataBlock = np.genfromtxt(qlpfile, \
#                              dtype=dtypeseq, delimiter=',',skip_header=1)
#    gtTIC = dataBlock['f0']
#    gtTOI = dataBlock['f1']
#    gtDisp = dataBlock['f2']
#    gtPer = dataBlock['f9']
#    gtEpc = dataBlock['f7']
#    gtDur = dataBlock['f11']
#    # Need to fix single transit TOIs with nan for period
#    idx = np.where(np.logical_not(np.isfinite(gtPer)))[0]
#    gtPer[idx] = 1000.0
    
    # TEV version
    # TEV csv has commas in strings
    # Use this sed 's/,"\(.*\),\(.*\)",/,"\1;\2",/' toi-plus-2019-03-15.csv
    # To fix string before reading in
    # As of Oct. 2019 I needed to use this to fix commas in strings
    # sed -e 's/""//g' -e 's/,"[^"]*/,"NOCOMMENT/g' csv-file-2019-10-29.csv > toi-plus-2019-10-29-fixed.csv
    qlpfile = 'csv-file-toi-catalog-FIXED-20221023.csv'
    dtypeseq = ['U20','U20','i4','f8','U2']
    dtypeseq.extend(['f8']*14)
    dtypeseq.extend(['U20','U80'])
    dtypeseq.extend(['f8']*12)
    dtypeseq.extend(['U20'])
    dtypeseq.extend(['i4']*7)
    dtypeseq.extend(['U40','U40'])
    dataBlock = np.genfromtxt(qlpfile, \
                              dtype=dtypeseq, delimiter=',',skip_header=5)
    gtTIC = dataBlock['f2']
    gtTOI = dataBlock['f3']
    gtDisp = dataBlock['f4']
    gtRA = dataBlock['f5']
    gtDec = dataBlock['f6']
    gtPer = dataBlock['f13']
    gtEpc = dataBlock['f11']
    print('Last TOI {0}'.format(np.max(gtTOI)))
    # Need to fix single transit TOIs with nan for period
    idx = np.where(np.logical_not(np.isfinite(gtPer)))[0]
    gtPer[idx] = 1000.0

    # Use the following to debug a particular target
#    idx = np.where(gtTIC == 279741379)[0]
#    gtTIC, gtTOI, gtDisp, gtPer, gtEpc, gtDur = cjb.idx_filter(idx, gtTIC, \
#                                gtTOI, gtDisp, gtPer, gtEpc, gtDur)

    # Load the tce data h5
    tceSeedInFile = 'sector-56_20221023_tce.h5'
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
    
    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=int)
    allatvalid = np.array([x.at_valid for x in all_tces], dtype=int)
    allra = np.array([x.raDeg for x in all_tces])
    alldec = np.array([x.decDeg for x in all_tces])
    allrow = np.array([x.row for x in all_tces])
    allcol = np.array([x.col for x in all_tces])
    allcam = np.array([x.camera for x in all_tces])
    allccd = np.array([x.ccd for x in all_tces])
    allpixvalid = np.array([x.pixposvalid for x in all_tces])
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
    fout.write('# Match {:s}\n'.format(qlpfile))
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
        if np.mod(i,20) == 0:
            print('Done {0:d} of {1:d}'.format(i, nGT))
        # If wideSearch True then query MAST for 
        #  all TICs within searchRad arcsec of this target
        print('Current Tic Search {0:d}'.format(curTic))
        if wideSearch:
            otherTICs = np.array(query_othertics(curTic, searchRad), dtype=np.int64)
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
#                idxfnd = np.where(minCadStrt == cadmap)[0]
                uowStartUse = timemap[idxfnd]
                idxfnd = np.argmin(np.abs(maxCadStrt - cadmap))
#                idxfnd = np.where(maxCadStrt == cadmap)[0]
                uowEndUse = timemap[idxfnd]
                print('Time: {0:f} {1:f}'.format(uowStartUse, uowEndUse))
            else:
                uowStartUse = uowStart
                uowEndUse = uowEnd
            bstpn, bsttic, bstMatch, bstStat, bstPeriodRatio, bstPeriodRatioFlag, bstFederateFlag, nFed = \
                    genericFed(curper, curepc, tryper, tryepc, trydur, trypn, trytic, uowStart, uowEnd)
            
        else:
            #print("No match in Ground truth")
            bstpn = 0
            bstMatch = -1
            bstStat = -1.0
            bstPeriodRatio = -1.0
            bstPeriodRatioFlag = 0
            bstFederateFlag = 0
            bsttic = 0
        str = '{:12d} {:8.2f} {:2s} {:12d} {:2d} {:2d} {:6.3f} {:2d} {:10.5f} {:2d}\n'.format(curTic, \
                   curToi, gtDisp[i], bsttic, bstpn, bstMatch, bstStat, bstPeriodRatioFlag, \
                   bstPeriodRatio, bstFederateFlag)
        if bstpn > 0:
            fout.write(str)
            print(str)
    fout.close()
    print('hello world')
