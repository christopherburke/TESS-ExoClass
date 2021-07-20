# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
from astropy.visualization import simple_norm
import fluxts_conditioning as flux_cond
import glob
import cjb_utils as cjb

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list

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

def assignEvents(t, epc, phi, per, phiDur):
    """Phase the data at period per and centered at to
       INPUT:
         phi - data phase running from -0.5<phi<=0.5
       OUTPUT:
         events - transit event number first appearing in data is event 1
     """
    events = np.zeros_like(phi)
    idx = np.where(np.abs(phi) <= phiDur)[0]
    if len(idx)>0:
        runcad = np.arange(len(t))
        minidx = np.min(idx)
        minto = t[minidx]
        idx2 = np.where(np.abs(t-minto) < 2.0*phiDur*per)[0]
        tmpcad = runcad[idx2]
        tmpphi = phi[idx2]
        ia = np.argmin(np.abs(tmpphi))
        cntrto = t[tmpcad[ia]]
        oldeve = np.round((cntrto-epc)/per)
        newcntrto = epc + oldeve*per
        events = np.round((t-newcntrto)/per) + 1.0
    return events


if __name__ == '__main__':
    # These are for parallel procoessing
    wID = 5
    nWrk = 6
    OVERWRITE = True
    #  Directory storing the ses mes time series
    sesMesDir = '/pdo/users/cjburke/spocvet/sector38'
    SECTOR = 6
    SECTOR1 = 384
    SECTOR2 = 38
#    sesMesDir = '/pdo/users/cjburke/spocvet/sector1-2'
#    SECTOR=-1

    #vetFile = 'spoc_sector1_early_fluxvet_20180904.txt'
    vetFile = 'spoc_fluxtriage_sector38_20210719.txt'
    tceSeedInFile = 'sector38_20210719_tce.pkl'
#    vetFile = 'spoc_sector1_2_fluxtriage_20181019.txt'
#    tceSeedInFile = 'sector1_2_20181019_tce.pkl'

    # Max number cadences closest to midtransit to go into  median depth estiamte
    MEDDEPN = 15
    # Search and filter parameters
    cadPerHr = 6
    firstFilterScaleFac = 10 # low frequency median filter will be
                            # firstFilterScaleFac*searchDurationHours medfilt window


    fin = open(tceSeedInFile, 'rb')
    all_tces = pickle.load(fin)
    fin.close()
    
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
    #idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
    #               (allvet == 1) )[0]
# DEBUG SINGLE TIC
#    idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
#                   (allvet == 1) & (alltic == 12862119 ) & (allpn==1))[0]
    idx = np.where((allatvalid == 1) & (alltrpvalid == 1) & (allsolarflux > 0.0) & \
                   (allvet == 1) )[0]
    alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
            allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
            allatrpdrstar, allatrpdrstare, allatadrstar = idx_filter(idx, \
            alltic, allpn, allatvalid, allrp, allrstar, alllogg, allper, alltmags, \
            allmes, allsnr, alldur, allsolarflux, allatdep, allatepoch, \
            allatrpdrstar, allatrpdrstare, allatadrstar)
            
    # Calculate modshift over flux triage passing TCEs
    for i, curTic in enumerate(alltic):
        print('{:d} of {:d} tic: {:d}'.format(i, len(alltic), curTic))
        if np.mod(i, nWrk) == wID:
            curPn = allpn[i]

            outFile = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_mods_diffImg_{0:016d}_{1:02d}_*.pdf'.format(curTic,curPn))
            outFileList = glob.glob(outFile)
            if OVERWRITE or not (len(outFileList)>0):
                fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(curTic,curPn))
                f = h5py.File(fileInput,'r')
                altDetrend = np.array(f['altDetrend'])
                validData = np.array(f['validData'])
                time = np.array(f['time'])
                cadNo = np.array(f['cadNo'])
                f.close()
                
                # Get the modshift trapezoid model fit
                fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_trpzdfit_{0:016d}_{1:02d}.txt'.format(curTic,curPn))
                dataBlock = np.genfromtxt(fileInput, dtype=['f8','f8','f8'])
                trpzdModel = dataBlock['f2']
                
                # Get trapezoid model on same times at time, cadNo
                tmpTrpzdModel = np.zeros_like(time)
                tmpTrpzdModel[validData] = trpzdModel

                # For each sector of tpf files that exist do the centroid                
                # Load the target pixel file
                for k in range(SECTOR1, SECTOR2+1):
                    fileInput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_tpf_{0:016d}_{1:02d}.h5d'.format(curTic,k))
                    if os.path.isfile(fileInput):
                        f = h5py.File(fileInput, 'r')
                        tpf_cad = np.array(f['cadenceNo'])
                        tpf_time = np.array(f['timetbjd'])
                        tpf_array = np.array(f['flux_array'])
                        tpf_bkg_array = np.array(f['flux_bkg_array'])
                        tpf_dq = np.array(f['dq_flag'])
                        tpf_vd = np.array(f['valid_data_flag'])
                        tpf_sat = np.array(f['saturate_pixel'])
                        tpf_medimg = np.array(f['median_image'])
                        centRow0 = f['1CRV4P'][0]
                        centCol0 = f['2CRV4P'][0]
                        ia, ib = cjb.intersect(cadNo, tpf_cad)
                        useTime = time[ia]
                        useVD = validData[ia]
                        useTrpzdModel = tmpTrpzdModel[ia]
                        
                        hasCent = False
                        fileCentroid = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_flxwcent_{0:016d}_{1:02d}_{2:02d}.h5d'.format(curTic,curPn, k))
                        if os.path.isfile(fileCentroid):
                            f = h5py.File(fileCentroid,'r')
                            centTime = np.array(f['time'])
                            centFlxw1 = np.array(f['flxw_centr1'])
                            centFlxw2 = np.array(f['flxw_centr2'])
                            centDqFlg = np.array(f['dqflgs'])
                            idx = np.where(centDqFlg == 0)[0]
                            centTime, centFlxw1, centFlxw2 = cjb.idx_filter(idx, \
                                    centTime, centFlxw1, centFlxw2)
                            hasCent = True

                        
                        #norm = simple_norm(tpf_medimg, 'log', percent=99.)
                        #plt.imshow(tpf_medimg, norm=norm, origin='lower', cmap='viridis')
                        #plt.colorbar()    
                        #plt.show()
                        usePhase = phaseData(useTime, allper[i], allatepoch[i])
                        phaseDur = alldur[i] / 24.0 / allper[i]
                        useEvents = assignEvents(useTime, allatepoch[i], usePhase, allper[i], phaseDur)
                        eventTime = usePhase + useEvents
                        # Mark data in transit as deweighted during detrending
                        ootvd = np.full_like(tpf_vd, True)
                        idx = np.where(np.abs(usePhase) < alldur[i]/2.0/24.0/allper[i])[0]
                        ootvd[idx] = False
                        ootvdGd = np.copy(ootvd)
                        idx = np.where(np.logical_not(useVD))[0]
                        ootvdGd[idx] = False
                        intvd = np.logical_not(ootvd)
                        intvdGd = np.copy(intvd)
                        intvdGd[idx] = False
                        # Check to see if there are any in transit points
                        #  If none skip this sector diff image generation
                        if (np.sum(intvdGd) > 0):
                            searchDurationHours = alldur[i]
                    
                            # Cap duration at 15 hours
                            searchDurationHours = np.min([searchDurationHours, 15.0])
                    
                    
                            imgShp = tpf_medimg.shape
                            nr = imgShp[0]
                            nc = imgShp[1]
                            depthImg = np.zeros((nr,nc), dtype=np.float)
                            snrImg = np.zeros((nr,nc), dtype=np.float)
                            priposImg = np.zeros((nr,nc), dtype=np.float)
                            doDebug = False
                            for ii in range(nr):
                                for jj in range(nc):
                                    curFlx = tpf_array[:,ii,jj]
                                    # Need to add a bias level 
                                    mnFlx = np.min(curFlx[useVD])
                                    print('MinFlux: {:f}'.format(mnFlx))
                                    # Catch nan in TPF
                                    if np.isfinite(mnFlx) and (len(np.where(tpf_vd)[0])>0):
                                        if mnFlx < 0.0:
                                            curFlx = curFlx - mnFlx + 10.0
                                            #plt.plot(curFlx, '.')
                                            #plt.show()
                                        tmpDebug = doDebug
                        #                if ii == 8 and jj == 5:
                        #                    plt.plot(eventTime, curFlx[validData], '.')
                        #                    plt.show()
                        #                    tmpDebug = True
                                        final_smooth_flux, bad_edge_flag = flux_cond.detrend_with_smoothn_edgefix(\
                                            curFlx, tpf_vd, ootvd, int(np.ceil(cadPerHr*searchDurationHours)), fixEdge=True, \
                                             medfiltScaleFac=10, gapThreshold=5, edgeExamWindow=8, \
                                             edgeSig=6.0, edgeMinCad=50, debug=tmpDebug)
                        #                if ii == 8 and jj == 5:
                        #                    print(phaseDur)
                        #                    plt.plot(usePhase[validData], final_smooth_flux[validData]-1.0, '.')
                        #                    plt.show()
                                        # Put detrended flux and trapezoid model into a temporary file                
                                        fileOutput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_cent_modtemp_{0:016d}_{1:02d}.txt'.format(curTic, curPn))
                                        np.savetxt(fileOutput, np.transpose([useTime[useVD], \
                                                                final_smooth_flux[useVD]-1.0, \
                                                                useTrpzdModel[useVD]]))
                                        pngOutputPrefix = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_cnt_modtemp_{0:02d}_{1:02d}_{2:016d}_{3:02d}_{4:02d}'.format(ii,jj, curTic, curPn, k))
                                        # Build argument list
                                        syscall = '/pdo/users/cjburke/spocvet/modshift/modshift {:s} {:s} {:016d}_{:02d} {:f} {:f} 1'.format(\
                                                            fileOutput, pngOutputPrefix, curTic, curPn, allper[i], allatepoch[i])
                                        p = Popen(syscall.split(), stdin=None, stdout=PIPE, stderr=PIPE)
                                        sysreturn, err = p.communicate()
                                        rc = p.returncode
                                        #print('alpha')
                                        retlist = sysreturn.split()[1:]
                                        try:
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
                                        except:
                                            primsig = 0.0
                                            fred = 1.0
                                            possig = 0.0
                        
                        
                                        depthImg[ii,jj] = primsig
                                        if fred > 0.0:
                                            snrImg[ii,jj] = primsig/fred
                                        else:
                                            snrImg[ii,jj] = primsig
                                        priposImg[ii,jj] = np.max([primsig - possig, 0.0])
                        
                        #                plt.cla()
                        #                plt.plot(usePhase[validData], final_smooth_flux[validData], '.')
                        #                plt.plot(usePhase[intvdGd], intflx, '.')
                        #                tmpy = np.ones_like(usePhase)
                        #                tmpx = np.copy(usePhase)
                        #                tmpy[intvd] = 1.0-intdep/1.0e6
                        #                ia = np.argsort(tmpx)
                        #                plt.plot(tmpx[ia], tmpy[ia], '-')
                        #                plt.ylim([np.min([1.0-1.2*allatdep[i]/1.0e6,1.0-errflx/1.0e6*3.5 ]), 1.0+errflx*3.5/1.0e6])
                        #                plt.show()
                                        
                                        
                                        #print('Row: {:d} of {:d} Col: {:d} of {:d}'.format(ii, imgShp[0], jj, imgShp[1]) )
                                        #plt.plot(usePhase, final_smooth_flux[validData], '.')
                                        #plt.show()
                            if (np.sum(np.isfinite(tpf_medimg.ravel())) > 10):                          
                                norm = simple_norm(tpf_medimg, 'log', percent=99.)
                                plt.subplot(2,2,1)
                                plt.imshow(tpf_medimg, norm=norm, origin='lower', cmap='viridis')
                                plt.colorbar() 
                                if hasCent:
                                    plt.plot(centFlxw1-centRow0, centFlxw2-centCol0,   '-c')
                                plt.title('Median Image S{0:02d}'.format(k))
                            
                            norm = simple_norm(depthImg, 'linear', percent=99.)
                            plt.subplot(2,2,3)
                            plt.imshow(depthImg, norm=norm, origin='lower', cmap='viridis')
                            plt.colorbar()   
                            if hasCent:
                                plt.plot(centFlxw1-centRow0, centFlxw2-centCol0,   '-c')
                            plt.title('Primary Sig S{0:02d}'.format(k))
                        
                            norm = simple_norm(snrImg, 'linear', percent=99.)
                            plt.subplot(2,2,4)
                            plt.imshow(snrImg, norm=norm, origin='lower', cmap='viridis')
                            plt.colorbar()   
                            plt.title('Primary Sig / Fred')
                            
                            norm = simple_norm(priposImg, 'linear', percent=99.)
                            plt.subplot(2,2,2)
                            plt.imshow(snrImg, norm=norm, origin='lower', cmap='viridis')
                            plt.colorbar()   
                            plt.title('Primary - Positive Sig')
                    
                            outFile = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'tess_mods_diffImg_{0:016d}_{1:02d}_{2:02d}.pdf'.format(curTic,curPn,k))
                            plt.savefig(outFile, format='pdf')
                            #plt.show()
                            plt.close()
                            print('alpha')
                            f.close()
            else:
                print('Skipping {0}'.format(outFile))
                    