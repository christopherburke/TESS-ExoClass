#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 21:13:43 2019

@author: cjburke
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import math
import os

if __name__ == '__main__':
    tic = 55092869
    pn = 1
    sector = 9
    sesDataDir = '/pdo/users/cjburke/spocvet/sector42/S09'

    epcDir = '{0:04d}'.format(int(math.floor(tic/1000.0)))
    localDir = os.path.join(sesDataDir,epcDir)
    fileInput = os.path.join(localDir, 'tess_sesmes_{0:016d}_{1:02d}.h5d'.format(tic,pn))
    f = h5py.File(fileInput,'r')

    phi = np.array(f['phi'])
    eve = np.array(f['events'])
    tphi = phi+eve
    vd = np.array(f['validData'])
    flx = np.array(f['initFlux'])
    altflx = np.array(f['altDetrend'])
    normTS = np.array(f['normTS'])
    corrTS = np.array(f['corrTS'])
    allSes = np.array(f['allSes'])
    allChases = np.array(f['allChases'])
    sesTS = corrTS/normTS
    cdppTS = 1.0e6/normTS
    print('SES')
    print(allSes)
    print('Chases')
    print(allChases)
    
    vd_ext = np.array(f['valid_data_flag_ext'])
    plt.plot(tphi, flx[vd], '.')
    plt.show()
    plt.plot(tphi, altflx[vd], '.')
    plt.show()
    plt.plot(tphi, sesTS[vd_ext], '.' )
    plt.show()
    plt.plot(tphi, cdppTS[vd_ext], '.')
    plt.show()
    
    