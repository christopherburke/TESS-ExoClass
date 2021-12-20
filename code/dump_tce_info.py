#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 09:25:25 2018
This routine performs detrending of the flux time series (see smoothn.py)
The calculates the single event statistics (SES) using a python port
of the Kepler pipeline TPS wavelet search methodology.  The SES timeseries
allows recalculating the transit SNR (or MES in TPS parlance) under alternative
detrending and allows metrics to be calculated to determine whether the
signal is consistent a transit and not a systematic.
In particular from the Kepler robovetter the SES/MES test is implemented
here and the CHASES test.  


@author: Christopher J. Burke
"""

import numpy as np
from gather_tce_fromdvxml import tce_seed


if __name__ == "__main__":
    # Load the h5 file that contains TCE seed information
    # The h5 file is created by gather_tce_fromdvxml.py
    tceSeedInFile = 'sector45_20211220_tce.h5'
    outFile = 'sector45_20211220_tce.txt'
    delim = ' | '
    tcedata = tce_seed()
    all_tces = tcedata.fill_objlist_from_hd5f(tceSeedInFile)
    
    # Loop over tces and write their info out
    cnt = 0
    fo = open(outFile, 'w')
    hdrout2 = """#Python Code to read in output
# import numpy as np
# dtypeseq = ['i4','i4','i4']
# dtypeseq.extend(['f8']*8)
# dtypeseq.extend(['i4']*3)
# dtypeseq.extend(['f8']*5)
# dtypeseq.extend(['i4'])
# dtypeseq.extend(['f8']*22)
# dtypeseq.extend(['i4'])
# dtypeseq.extend(['f8']*7)
# dtypeseq.extend(['i4'])
# dtypeseq.extend(['f8']*9)
# dataBlock = np.genfromtxt('junk.txt', dtype=dtypeseq, delimiter='|', skip_header=75)
# tic = dataBlock['f0']
# pn = dataBlock['f1']
# totnp = dataBlock['f2']
# dec = dataBlock['f3']
# ra = dataBlock['f4']
# teff = dataBlock['f5']
# teff_e = dataBlock['f6']
# feh = dataBlock['f7']
# logg = dataBlock['f8']
# rstar = dataBlock['f9']
# rstar_e = dataBlock['f10']
# sector = dataBlock['f11']
# ccd = dataBlock['f12']
# camera = dataBlock['f13']
# row = dataBlock['f14']
# column = dataBlock['f15']
# pmra = dataBlock['f16']
# pmdec = dataBlock['f17']
# tmag = dataBlock['f18']
# trn_valid = dataBlock['f19']
# trn_snr = dataBlock['f20']
# trn_epc = dataBlock['f21']
# trn_epc_e = dataBlock['f22']
# trn_rp = dataBlock['f23']
# trn_rp_e = dataBlock['f24']
# trn_imp = dataBlock['f25']
# trn_dur = dataBlock['f26']
# trn_dep = dataBlock['f27']
# trn_dep_e = dataBlock['f28']
# trn_per = dataBlock['f29']
# trn_per_e = dataBlock['f30']
# trn_rpdrstar = dataBlock['f31']
# trn_rpdrstar_e = dataBlock['f32']
# trn_adrstar = dataBlock['f33']
# trn_eqtemp = dataBlock['f34']
# trn_effflx = dataBlock['f35']
# ghst_coreval = dataBlock['f36']
# ghst_coresig = dataBlock['f37']
# ghst_haloval = dataBlock['f38']
# ghst_halosig = dataBlock['f39']
# trn_chi = dataBlock['f40']
# tce_gof = dataBlock['f41']
# tce_ntran = dataBlock['f42']
# tce_chi = dataBlock['f43']
# tce_epc = dataBlock['f44']
# tce_mes = dataBlock['f45']
# tce_maxsesmes = dataBlock['f46']
# tce_per = dataBlock['f47']
# tce_robstat = dataBlock['f48']
# tce_pulsedur = dataBlock['f49']
# trp_valid = dataBlock['f50']
# trp_snr = dataBlock['f51']
# trp_epc = dataBlock['f52']
# trp_dur = dataBlock['f53']
# trp_dep = dataBlock['f54']
# cent_oot_off = dataBlock['f55']
# cent_oot_off_e = dataBlock['f56']
# cent_tic_off = dataBlock['f57']
# cent_tic_off_e = dataBlock['f58']
# oddevensig = dataBlock['f59']
# --------------------------------------------------
"""
    fo.write(hdrout2)
    
    
    
    hdrout = '#TIC | PN | TotNumPlan | Dec  | RA | Teff[K]  | Teff_e[K] | '
    hdrout = hdrout + '[Fe/H] | logg | Rstar[Rsun] | Rstar_e | Sector | '
    hdrout = hdrout + 'CCD | Camera | Row[pix] | Column[pix] | pmRA | '
    hdrout = hdrout + 'pmDec | Tmag | TrnFitValid | Trn_SNR | Trn_Epc[TBJD] | '
    hdrout = hdrout + 'Trn_Epc_e | Trn_Rp[Rearth] | Trn_Rp_e | Trn_Imp | '
    hdrout = hdrout + 'Trn_Dur[Hr] | Trn_Dep[ppm] | Trn_Dep_e | Trn_Period | '
    hdrout = hdrout + 'Trn_Period_e | Trn_RpDRstar | Trn_RpDRstar_e | '
    hdrout = hdrout + 'Trn_aDRstar | Trn_EqTemp[K] | Trn_EffFlx[InsolEarth] | '
    hdrout = hdrout + 'Ghst_CoreVal | Ghst_CoreSig | Ghst_HaloSig | Ghst_HaloVal | '
    hdrout = hdrout + 'Trn_Chi | TCE_GOF | TCE_Ntran | TCE_Chi | TCE_Epoch | '
    hdrout = hdrout + 'MES | maxMESinSES | TCE_Period | RobStat | PulseDur | '
    hdrout = hdrout + 'Trp_Valid | Trp_SNR | Trp_Epoch[TBJD] | Trp_Dur[Hr] | '
    hdrout = hdrout + 'Trp_Dep[ppm] | Cent_OOT_Off | Cent_OOT_Off_e | '
    hdrout = hdrout + 'Cent_TIC_Off | Cent_TIC_Off_e | OddEvenSig'
    fo.write('{0}\n'.format(hdrout))



    for td in all_tces:        
        sout = '{0:16d}'.format(td.epicId)  #1
        print(cnt)
        cnt = cnt+1
        sout = sout + delim + '{0:2d}'.format(td.planetNum)
        
        sout = sout + delim + '{0:2d}'.format(td.totPlanetNum )
        # Target information
        sout = sout + delim + '{0:10.6f}'.format(td.decDeg )
        sout = sout + delim + '{0:10.6f}'.format(td.raDeg) #5 
        sout = sout + delim + '{0:6.1f}'.format(td.teff)
        sout = sout + delim + '{0:5.1f}'.format(td.teff_e) 
        sout = sout + delim + '{0:6.3f}'.format(td.feh )
        sout = sout + delim + '{0:6.3f}'.format(td.logg)
        sout = sout + delim + '{0:8.3f}'.format(td.rstar)  #10
        sout = sout + delim + '{0:7.3f}'.format(td.rstar_e) 
        sout = sout + delim + '{0:3d}'.format(td.sector )
        sout = sout + delim + '{0:1d}'.format(td.ccd )
        sout = sout + delim + '{0:1d}'.format(td.camera) 
        sout = sout + delim + '{0:7.2f}'.format(td.row ) #15
        sout = sout + delim + '{0:7.2f}'.format(td.col )
        sout = sout + delim + '{0:9.3f}'.format(td.pmra )
        sout = sout + delim + '{0:9.3f}'.format(td.pmdec )
        sout = sout + delim + '{0:6.3f}'.format(td.tmag )
        # All transit fit
        sout = sout + delim + '{0:1d}'.format(td.at_valid)  #20
        sout = sout + delim + '{0:8.3f}'.format(td.at_snr )
        sout = sout + delim + '{0:11.5f}'.format(td.at_epochbtjd) 
        sout = sout + delim + '{0:9.5f}'.format(td.at_epochbtjd_e )
        sout = sout + delim + '{0:8.2f}'.format(td.at_rp )
        sout = sout + delim + '{0:7.2f}'.format(td.at_rp_e) #25
        sout = sout + delim + '{0:6.3f}'.format(td.at_imp )
        sout = sout + delim + '{0:7.3f}'.format(td.at_dur )
        sout = sout + delim + '{0:8.1f}'.format(td.at_depth )
        sout = sout + delim + '{0:7.1f}'.format(td.at_depth_e) 
        sout = sout + delim + '{0:11.6f}'.format(td.at_period ) #30
        sout = sout + delim + '{0:9.6f}'.format(td.at_period_e) 
        sout = sout + delim + '{0:9.6f}'.format(td.at_rpDrstar )
        sout = sout + delim + '{0:9.6f}'.format(td.at_rpDrstar_e) 
        sout = sout + delim + '{0:8.3f}'.format(td.at_aDrstar )
        sout = sout + delim + '{0:7.1f}'.format(td.at_eqtemp ) #35
        sout = sout + delim + '{0:7.2f}'.format(td.at_effflux )
        # Ghost diag
        sout = sout + delim + '{0:8.3f}'.format(td.ghostcoreval) 
        sout = sout + delim + '{0:8.3f}'.format(td.ghostcoresig )
        sout = sout + delim + '{0:8.3f}'.format(td.ghosthalosig )
        sout = sout + delim + '{0:8.3f}'.format(td.ghosthaloval ) #40
        # TCE Infor
        tmp = td.modchi2
        if not np.isfinite(tmp):
            sout = sout + delim + '{0:8.2f}'.format(99999.99)
        else:
            sout = sout + delim + '{0:8.2f}'.format(td.modchi2 ) #41
        sout = sout + delim + '{0:8.3f}'.format(td.gof )
        sout = sout + delim + '{0:4d}'.format(td.ntran)
        sout = sout + delim + '{0:8.3f}'.format(td.chi2 )
        sout = sout + delim + '{0:11.5f}'.format(td.tce_epoch) #45
        sout = sout + delim + '{0:9.3f}'.format(td.mes )
        sout = sout + delim + '{0:8.3f}'.format(td.maxsesinmes) 
        sout = sout + delim + '{0:11.6f}'.format(td.tce_period )
        sout = sout + delim + '{0:9.3f}'.format(td.robstat )
        sout = sout + delim + '{0:7.3f}'.format(td.pulsedur ) #50
        # Trpzed fit
        sout = sout + delim + '{0:1d}'.format(td.trp_valid )
        sout = sout + delim + '{0:9.3f}'.format(td.trp_snr )
        sout = sout + delim + '{0:11.5f}'.format(td.trp_epochbtjd) 
        sout = sout + delim + '{0:7.3f}'.format(td.trp_dur )
        sout = sout + delim + '{0:8.1f}'.format(td.trp_depth)  #55
        # Centroid Fits
        sout = sout + delim + '{0:7.2f}'.format(td.cent_oot_offset) 
        sout = sout + delim + '{0:6.2f}'.format(td.cent_oot_offset_e) 
        sout = sout + delim + '{0:7.2f}'.format(td.cent_tic_offset )
        sout = sout + delim + '{0:6.2f}'.format(td.cent_tic_offset_e )
        # Odd/Even
        sout = sout + delim + '{0:8.3f}'.format(td.oe_signif) #60
        
        fo.write('{0}\n'.format(sout))
