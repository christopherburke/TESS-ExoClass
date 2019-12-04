#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 22:24:16 2019

@author: cjburke
"""

import numpy as np

if __name__ == '__main__':
    
    dtypeseq = ['i4','i4','i4']
    dtypeseq.extend(['f8']*8)
    dtypeseq.extend(['i4']*3)
    dtypeseq.extend(['f8']*5)
    dtypeseq.extend(['i4'])
    dtypeseq.extend(['f8']*22)
    dtypeseq.extend(['i4'])
    dtypeseq.extend(['f8']*7)
    dtypeseq.extend(['i4'])
    dtypeseq.extend(['f8']*9)
    dataBlock = np.genfromtxt('junk.txt', dtype=dtypeseq, delimiter='|', \
                              skip_header=1)
    tic = dataBlock['f0']
    pn = dataBlock['f1']
    totnp = dataBlock['f2']
    dec = dataBlock['f3']
    ra = dataBlock['f4']
    teff = dataBlock['f5']
    teff_e = dataBlock['f6']
    feh = dataBlock['f7']
    logg = dataBlock['f8']
    rstar = dataBlock['f9']
    rstar_e = dataBlock['f10']
    sector = dataBlock['f11']
    ccd = dataBlock['f12']
    camera = dataBlock['f13']
    row = dataBlock['f14']
    column = dataBlock['f15']
    pmra = dataBlock['f16']
    pmdec = dataBlock['f17']
    tmag = dataBlock['f18']
    trn_valid = dataBlock['f19']
    trn_snr = dataBlock['f20']
    trn_epc = dataBlock['f21']
    trn_epc_e = dataBlock['f22']
    trn_rp = dataBlock['f23']
    trn_rp_e = dataBlock['f24']
    trn_imp = dataBlock['f25']
    trn_dur = dataBlock['f26']
    trn_dep = dataBlock['f27']
    trn_dep_e = dataBlock['f28']
    trn_per = dataBlock['f29']
    trn_per_e = dataBlock['f30']
    trn_rpdrstar = dataBlock['f31']
    trn_rpdrstar_e = dataBlock['f32']
    trn_adrstar = dataBlock['f33']
    trn_eqtemp = dataBlock['f34']
    trn_effflx = dataBlock['f35']
    ghst_coreval = dataBlock['f36']
    ghst_coresig = dataBlock['f37']
    ghst_haloval = dataBlock['f38']
    ghst_halosig = dataBlock['f39']
    trn_chi = dataBlock['f40']
    tce_gof = dataBlock['f41']
    tce_ntran = dataBlock['f42']
    tce_chi = dataBlock['f43']
    tce_epc = dataBlock['f44']
    tce_mes = dataBlock['f45']
    tce_maxsesmes = dataBlock['f46']
    tce_per = dataBlock['f47']
    tce_robstat = dataBlock['f48']
    tce_pulsedur = dataBlock['f49']
    trp_valid = dataBlock['f50']
    trp_snr = dataBlock['f51']
    trp_epc = dataBlock['f52']
    trp_dur = dataBlock['f53']
    trp_dep = dataBlock['f54']
    cent_oot_off = dataBlock['f55']
    cent_oot_off_e = dataBlock['f56']
    cent_tic_off = dataBlock['f57']
    cent_tic_off_e = dataBlock['f58']
    oddevensig = dataBlock['f59']
    
    print('help')
    