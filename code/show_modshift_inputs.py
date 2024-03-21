import numpy as np
import matplotlib.pyplot as plt
import argparse
from tec_used_params import tec_use_params
import os
import math

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
    parser = argparse.ArgumentParser()
    parser.add_argument("t", type=int, \
                        help="TIC ID to examine")
    parser.add_argument("n", type=int,\
                        default = 1, \
                        help="Planet Number to examine")
    args = parser.parse_args() 
    tID = args.t
    PN = args.n
    tp = tec_use_params()

    #  Directory storing the resampled dv time series data
    dvDataDir = '/pdo/users/cjburke/spocvet/{0}'.format(tp.tecdir)
    # Directory of output
    outputDir = dvDataDir
    SECTOR = tp.sector
    fileInput = os.path.join(make_data_dirs(outputDir, SECTOR, tID), 'tess_trpzdfit_{0:016d}_{1:02d}_1.txt'.format(tID,PN))
    print(fileInput)

    dataBlock = np.genfromtxt(fileInput, dtype=['f8']*3)
    ts = dataBlock['f0']
    flx = dataBlock['f1']
    fmod = dataBlock['f2']

    print('ts non-Finite: ',np.any(~np.isfinite(ts)))
    print('flx non-Finite: ',np.any(~np.isfinite(flx)))
    print('fmod non-Finite: ',np.any(~np.isfinite(fmod)))


    plt.plot(ts, flx, '.')
    plt.plot(ts, fmod, '.')
    plt.show()




