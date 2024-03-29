# -*- coding: utf-8 -*-
"""
Generate twexo pdfs for all targets
"""

import math
import os
from subprocess import Popen, PIPE
import numpy as np
import argparse
from tec_used_params import tec_use_params

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

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list


if __name__ == '__main__':
    # Parse the command line arguments for multiprocessing
    # With Gnu parallel with 3 cores
    # seq 0 2 | parallel --results gen_twexo_results python gen_twexo.py -w {} -n 3
    #  In practice, run on two separate pdo machine

    parser = argparse.ArgumentParser()
    parser.add_argument("-w", type=int,\
                        default = 0, \
                        help="Worker ID Number 0 through nWrk-1")
    parser.add_argument("-n", type=int,\
                        default = 1, \
                        help="Number of Workers")
    

    args = parser.parse_args() 
    # These are for parallel procoessing
    wID = int(args.w)
    nWrk = int(args.n)

    tp = tec_use_params()

    #  Directory storing the ses mes time series
    sesMesDir = '/pdo/users/cjburke/spocvet/{0}'.format(tp.tecdir)
    SECTOR = tp.sector

    doPDFs = True
    vetFile = 'spoc_fluxtriage_{0}.txt'.format(tp.tecfile)
    overwrite = False

    # Load the  flux vetting
    dataBlock = np.genfromtxt(vetFile, dtype=[int,int,int,'S1'])
    fvtic = dataBlock['f0']
    fvpn = dataBlock['f1']
    fvvet = dataBlock['f2']

    # Generate html for triage passing TICs
    idx = np.where(fvvet == 1)[0]
    useTic = np.unique(fvtic[idx])
    nTic = len(useTic)
    for i,curTic in enumerate(useTic):
        if np.mod(i, nWrk) == wID:
            if np.mod(i, 20) == 0:
                print('Done {0:d} of {1:d}'.format(i, nTic))
            htmlOutput = os.path.join(make_data_dirs(sesMesDir, SECTOR, curTic), 'twexo_{0:016d}'.format(curTic))
            if not os.path.isfile(htmlOutput + '.pdf') or overwrite:
                # Build argument list
                syscall = 'python twexo.py -t {0:d} -of {1} -nw'.format(\
                                curTic, htmlOutput)
                print(syscall)
                p = Popen(syscall.split(), stdin=None, stdout=PIPE, stderr=PIPE)
                sysreturn, err = p.communicate()
                if os.path.isfile(htmlOutput + '.html'):
                    syscall = 'wkhtmltopdf {0}.html {0}.pdf'.format(htmlOutput)
                    print(syscall)
                    p = Popen(syscall.split(), stdin=None, stdout=PIPE, stderr=PIPE)
                    sysreturn, err = p.communicate()
                else:
                    print(syscall + ' failed to produce output')
                    print(sysreturn)
                    print(err)
