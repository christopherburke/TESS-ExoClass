# -*- coding: utf-8 -*-
"""
Pulls the difference image information out of the full DV report
and makes them available as individual pages for the TEC report

AUTHOR: Christopher J. Burke
"""

import numpy as np
from gather_tce_fromdvxml import tce_seed
import os
from subprocess import Popen, PIPE
import math
import glob

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



if __name__ == '__main__':
    # These are for parallel procoessing
    wID = 0
    nWrk = 1
    
    summaryFolder = '/pdo/spoc-data/sector-046/dv-reports'
    summaryPrefix = 'tess2021337012456-'
    summaryPostfix = '-00547_dvr.pdf'
    SECTOR1 = 46
    SECTOR2 = 46
    multiRun = False
    if SECTOR2 - SECTOR1 > 0:
        multiRun = True
    tceSeedInFile = 'sector46_20220118_tce.h5'
    sesMesDir = '/pdo/users/cjburke/spocvet/sector46'
    SECTOR = 46
    overwrite = False
    
    # Load the tce data h5
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)    

    alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
    allpn = np.array([x.planetNum for x in all_tces], dtype=int)
    #idx = np.where(alltic == 167600516)[0]
    #alltic = alltic[idx[0]:]
    #allpn = allpn[idx[0]:]
    for i in range(len(alltic)):
        if np.mod(i, nWrk) == wID:
            curTic = alltic[i]
            print(curTic, i, len(alltic))
            curPN = allpn[i]
            srchstr = '{0}s{1:04d}-s{2:04d}-{3:016d}{4}'.format(summaryPrefix,SECTOR1,SECTOR2,curTic,'*dvr.pdf')
            dvReportFileList = glob.glob(os.path.join(summaryFolder,srchstr))
            if not len(dvReportFileList)==1:
                if len(dvReportFileList) == 0:
                    print('EXITING! worker {0:d} of {1:d} Could not find {2}'.format(wID, nWrk, srchstr))
                    exit()
                else:
                    print('EXITING! worker {0:d} of {1:d} Found multiple files from {2}'.format(wID, nWrk, srchstr))
                    exit()
            dvReportFile = dvReportFileList[0]
    #        comstring = 'pdftotext -layout {0} - | grep -A 12 \"Difference image for target {1:d}, planet candidate {2:d}\" | tail -n 1'.format(dvReportFile, curTic, curPN)
    
            # Need to also determine number of contents pages before page 1
            pdftotext_com = 'pdftotext -layout {0} - '.format(dvReportFile)
            grep_com = ['grep', '-B', '1', 'SUMMARY']
            p1 = Popen(pdftotext_com.split(), stdout=PIPE)
            p2 = Popen(grep_com, stdin=p1.stdout, stdout=PIPE)
            p1.stdout.close()
            sysreturn, err = p2.communicate()
            rc = p2.returncode
            retlist = sysreturn.split(b'\n')
            pgistr = retlist[0].strip(b' ').decode('ascii')
            prePages = 0
            if pgistr == 'ii':
                prePages = 2
            if pgistr == 'iii':
                prePages = 3
            if pgistr == 'iv':
                prePages = 4
            if pgistr == 'v':
                prePages = 5
            if pgistr == 'vi':
                prePages = 6
            if pgistr == 'vii':
                prePages = 7
            if pgistr == 'viii':
                prePages = 8
            if prePages == 0:
                print('EXITING! worker {0:d} of {1:d} found unexpected number of prepages {2}'.format(wID, nWrk, srchstr))
                exit()

            #prePages = len(retlist[0].split('i'))-1
            if multiRun: # There is a summary centroid plot get its page and save it out
                grep_com = ['grep','-A','5','planet-{0:02d}/difference-image/{1:016d}-{0:02d}-difference-image-centroid-offsets.fig'.format(curPN,curTic)]
                p1 = Popen(pdftotext_com.split(), stdout=PIPE)
                p2 = Popen(grep_com, stdin=p1.stdout, stdout=PIPE)
                p1.stdout.close()      
                sysreturn, err = p2.communicate()
                rc = p2.returncode
                retlist = sysreturn.split(b'\n')
                pageWant = -1
                if len(retlist) > 1: # If has difference image
        #            print(retlist)
                    pageWant = int(retlist[-2])
                    pageWant = prePages + pageWant
                    
                    dvDiffFile = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_diffImg_{0:016d}_{1:02d}_centsum.pdf'.format(curTic,curPN))
                    if (not os.path.isfile(dvDiffFile)) and (not overwrite):
                        gs_com = 'gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage={0:d} -dLastPage={0:d} -sOutputFile={2} {1}'.format(pageWant, dvReportFile, dvDiffFile)
                        p1 = Popen(gs_com.split(), stdout=PIPE)
                        sysreturn, err = p1.communicate()
                        rc = p1.returncode
            for curSector in np.arange(SECTOR1, SECTOR2+1):
                grep_com = ['grep','-A','5','planet-{0:02d}/difference-image/{1:016d}-{0:02d}-difference-image-{2:02d}'.format(curPN,curTic,curSector)]
                p1 = Popen(pdftotext_com.split(), stdout=PIPE)
                p2 = Popen(grep_com, stdin=p1.stdout, stdout=PIPE)
                p1.stdout.close()
                sysreturn, err = p2.communicate()
                rc = p2.returncode
                retlist = sysreturn.split(b'\n')
                pageWant = -1
                if len(retlist) > 1: # If has difference image
                    pageWant = int(retlist[-2])
                    pageWant = prePages + pageWant
                    
                    dvDiffFile = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_diffImg_{0:016d}_{1:02d}_{2:02d}.pdf'.format(curTic,curPN,curSector))
                    if (not os.path.isfile(dvDiffFile)) and (not overwrite):
                        gs_com = 'gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage={0:d} -dLastPage={0:d} -sOutputFile={2} {1}'.format(pageWant, dvReportFile, dvDiffFile)
                        p1 = Popen(gs_com.split(), stdout=PIPE)
                        sysreturn, err = p1.communicate()
                        rc = p1.returncode
    
