# -*- coding: utf-8 -*-
"""
reports status and file counts for the TEC steps.
Used to make sure all files exist as expected

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
    SECTOR1 = 1
    SECTOR2 = 46
    multiRun = False
    if SECTOR2 - SECTOR1 > 0:
        multiRun = True
    tceSeedInFile = 'sector1-46_20220328_tce.h5'
    sesMesDir = '/pdo/users/cjburke/spocvet/sector1-46'
    SECTOR = -1
    skyline_out = 'skyline_data_sector1-46_20220328.txt'
    fed_knownP_out = 'federate_knownP_sector1-46_20220328.txt'
    fed_toi_out = 'federate_toiWtce_sector1-46_20220328.txt'
    fed_self_out = 'selfMatch_sector1-46_20220328.txt'
    fluxVetOut = 'spoc_fluxtriage_sector1-46_20220328.txt'
    SWEETMAXPER = 5.0
    sweet_out = 'spoc_sweet_sector1-46_20220328.txt'
    modump_out = 'spoc_modump_sector1-46_20220328.txt'    
    fileOut1 = 'spoc_ranking_Tier1_sector1-46_20220328.txt'
    fileOut2 = 'spoc_ranking_Tier2_sector1-46_20220328.txt'
    fileOut3 = 'spoc_ranking_Tier3_sector1-46_20220328.txt'
    if multiRun:
        useSector = 1000+SECTOR2
    else:
        useSector = SECTOR2

    
    # Check that TCE Seed File exists
    prereq = 0
    if (os.path.isfile(tceSeedInFile)):
        # Load the tce data h5
        tcedata = tce_seed()
        all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
        alltic = np.array([x.epicId for x in all_tces], dtype=np.int64)
        all_pns = np.array([x.planetNum for x in all_tces], dtype=int)
        allsolarflux = np.array([x.at_effflux for x in all_tces])
        alltrpvalid = np.array([x.trp_valid for x in all_tces])
        allatvalid = np.array([x.at_valid for x in all_tces], dtype=int)
        allper = np.array([x.at_period for x in all_tces])

        allUnqTic = np.unique(alltic)
        prereq = 1
        print('TCE Seed File Exists Num TCE: {0:d}'.format(len(alltic)))

    else:
        print('TCE Seed Does Not Exist!')
        
    # Check That Light Curve files exist
    NLC = 0
    missing = []
    if prereq > 0: # Need TCE Seed file
        NLCExp = len(alltic)
        for i, curTic in enumerate(alltic):
            curPn = all_pns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_dvts_{0:016d}_{1:02d}.h5d'.format(curTic, curPn))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} LC files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_lc_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()
    
    # Check that skyline file exists
    if prereq > 0:
        if os.path.isfile(skyline_out):
            print('Skyline File Exists')
        else:
            print('Skyline File Does NOT Exist!')

    # Check that federation with Known planets exists
    if prereq > 0:
        if os.path.isfile(fed_knownP_out):
            print('Federation with Known Planets File Exists')
        else:
            print('Federation with Known Planets Does NOT Exist!')
            
    # Check that federation with TOIs exists
    if prereq > 0:
        if os.path.isfile(fed_toi_out):
            print('Federation with TOIs File Exists')
        else:
            print('Federation with TOIs Does NOT Exist!')

    # Check that federation with self TCEs exists
    if prereq > 0:
        if os.path.isfile(fed_self_out):
            print('Federation with TCEs File Exists')
        else:
            print('Federation with TCEs Does NOT Exist!')

    # Check That SES time series files exist
    NLC = 0
    missing = []
    if prereq > 0: # Need TCE Seed file
        NLCExp = len(alltic)
        for i, curTic in enumerate(alltic):
            curPn = all_pns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(curTic,curPn))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} SES time series files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_ses_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()

    # Check That DV difference image pdf exists
    NLC = 0
    missing = []
    if prereq > 0: # Need TCE Seed file
        NLCExp = len(alltic)
        for i, curTic in enumerate(alltic):
            curPn = all_pns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_diffImg_{0:016d}_{1:02d}_*.pdf'.format(curTic,curPn))
            fileList = glob.glob(fileInput)
            if len(fileList) > 0:
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} DV Diff Image files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_dvdiff_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()

    # Check That target pixel files exists
    NLC = 0
    missing = []
    if prereq > 0: # Need TCE Seed file
        NLCExp = len(allUnqTic)
        for i, curTic in enumerate(allUnqTic):
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_tpf_{0:016d}_*.h5d'.format(curTic))
            fileList = glob.glob(fileInput)
            if len(fileList) > 0:
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} TPF files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_tpf_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()

    # Check that triage file exists
    if prereq > 0:
        if os.path.isfile(fluxVetOut):
            print('Triage File Exists')
            prereq = 2
            # Load the  flux vetting results and merge into tce seed vector
            dataBlock = np.genfromtxt(fluxVetOut, dtype=[int,int,int,'S1'])
            fvtic = dataBlock['f0']
            fvpn = dataBlock['f1']
            fvvet = dataBlock['f2']
            
            allvet = np.zeros_like(all_pns)
            for i in range(len(allvet)):
                idx = np.where((alltic[i] == fvtic) & (all_pns[i] == fvpn))[0]
                if len(idx) > 0:
                    allvet[i] = fvvet[idx]
        else:
            print('Triage File Does NOT Exist!')

    # Check that modshift alt detrend exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have modshift run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        for i, curTic in enumerate(curTics):
            curPn = curPns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_{0:016d}_{1:02d}-modshift.pdf'.format(curTic,curPn))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} Modshift Alt Detrend files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_modalt_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()

    # Check that modshift median detrend exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have modshift run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        for i, curTic in enumerate(curTics):
            curPn = curPns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_{0:016d}_{1:02d}_med-modshift.pdf'.format(curTic,curPn))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} Modshift Median Detrend files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_modmed_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()
        
    # Check that sweet test exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have sweet test run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1) & (allper<SWEETMAXPER))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        # Read in Sweet test results
        if os.path.isfile(sweet_out):
            dtypeseq=['i4','i4']
            dtypeseq.extend(['f8']*17)
            dataBlock = np.genfromtxt(sweet_out, dtype=dtypeseq)
            swTic = dataBlock['f0']
            NLC = len(swTic)
            print('{0:d} of {1:d} sweet test results found'.format(NLC, NLCExp))
        else:
            print('No Sweet test results found')
    # Check that the flux pdc statistics were gathered exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have modshift run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        for i, curTic in enumerate(curTics):
            curPn = curPns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_flxwcent_{0:016d}_{1:02d}_*.h5d'.format(curTic,curPn))
            fileList = glob.glob(fileInput)
            if len(fileList) > 0:
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} Flux PDC stats files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_flxwcent_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()
    
    # Check that twexo exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have sweet test run
        idx = np.where(allvet == 1)[0]
        curTics = alltic[idx]
        curUnqTics = np.unique(curTics)
        NLCExp = len(curUnqTics)
        for curTic in curUnqTics:
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'twexo_{0:016d}.pdf'.format(curTic))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} Twexo files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_twexo_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()

    # Check that momentum dump data exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have sweet test run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1) )[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        # load the momentum dump transit fraction data
        if os.path.isfile(modump_out):
            dtypeseq=['i4','i4','f8']
            dataBlock = np.genfromtxt(modump_out, dtype=dtypeseq)
            mdTic = dataBlock['f0']
            NLC = len(mdTic)
            print('{0:d} of {1:d} momentum dump results found'.format(NLC, NLCExp))
        else:
            print('No momentum dump file available')
    # Check that basic centroidingexists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have basic centroid
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        for i, curTic in enumerate(curTics):
            curPn = curPns[i]
            localDir = make_data_dirs(sesMesDir,SECTOR,curTic)
            fileInput = os.path.join(localDir, 'tess_bsc_diffImg_{0:016d}_{1:02d}_*.pdf'.format(curTic,curPn))
            fileList = glob.glob(fileInput)
            if len(fileList) > 0:
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} Basic Centroid files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_bsccen_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()

    # Check that Tier file exists
    if prereq > 1:
        if os.path.isfile(fileOut1):
            print('Tier1 File Exists')
        else:
            print('Tier1 File Does NOT Exist!')
    if prereq > 1:
        if os.path.isfile(fileOut2):
            print('Tier2 File Exists')
        else:
            print('Tier2 File Does NOT Exist!')
    if prereq > 1:
        if os.path.isfile(fileOut3):
            print('Tier3 File Exists')
        else:
            print('Tier3 File Does NOT Exist!')
            
    # Check that TEC report exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have modshift run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        NLCExp = len(curTics)
        for i, curTic in enumerate(curTics):
            curPn = curPns[i]
            fileInput = os.path.join(sesMesDir,'pdfs', 'tec-s{3:04d}-{1:016d}-{2:02d}.pdf'.format(i, curTic, curPn, SECTOR2))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} TEC reports files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_tecrep_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()
    
    # Check that TEV merged report exists
    NLC = 0
    missing = []
    if prereq > 1:
        # Only subset of targets have modshift run
        idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1))[0]
        curTics = alltic[idx]
        curPns = all_pns[idx]
        curUnqTics = np.unique(curTics)
        NLCExp = len(curUnqTics)
        for i, curTic in enumerate(curUnqTics):
            fileInput = os.path.join(sesMesDir,'tevpdfs', 'tec-s{0:04d}-{1:016d}-00001_dvm.pdf'.format(useSector, curTic))
            if os.path.isfile(fileInput):
                NLC += 1
            else:
                missing.append(fileInput)
        print('{0:d} of {1:d} TEV Merged reports files found'.format(NLC, NLCExp))
        if not NLC == NLCExp: # If missing files write out their names
            fm = open('tec_tevmrg_missing.txt', 'w')
            for missFile in missing:
                fm.write('{0}\n'.format(missFile))
            fm.close()
