import numpy as np
import argparse
import os
from datetime import datetime
import glob

class tec_param:
    def __init__(self):
        self.sector = 0
        self.sector1 = 0
        self.sector2 = 0
        self.tecdir = 'sector'
        self.spocdir = 'sector-000'
        self.tecfile = 'sector-00_230101'
        self.toifile = 'FIXED-230101'
        self.lcpre = 'tess2022357055054-s0060-'
        self.lcnum = '0249'
        self.dvpre = 'tess2019199201929-'
        self.dvnum = '00699'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sector1", type=int, nargs=1,\
                        help="Sector Number or First Sector Number if Multi-sector [int]")
    parser.add_argument("sector2", type=int, nargs='?',\
                        default=None,\
                        help="Not used for Single Sector. For Multi-sector,  last Sector Number [int]")
    parser.add_argument("-spocdir", default=None,\
                        help="Not used for single sector. For Multi-sector, specify specific spoc directory of data (i.e., sector-014-026+040-060)")

    args = parser.parse_args()
    spocPathPrefix = '/pdo/spoc-data'
    tecPathPrefix = '/pdo/users/cjburke/spocvet'

    # Get date string
    current_datetime = datetime.now()
    usedate = current_datetime.strftime("%Y%m%d")

    tp = tec_param()
    # Determine whether it is a single or multi-sector run
    if args.sector2 is None:
        # Single Sector
        print('Single Sector')
        tp.sector = args.sector1[0] # nargs with integer produces list
        tp.sector1 = args.sector1[0]
        tp.sector2 = args.sector1[0]
        tp.tecdir = 'sector{0:d}'.format(args.sector1[0])
        tp.spocdir= 'sector-{0:03d}'.format(args.sector1[0])
        tp.tecfile = 'sector-{0:d}_{1}'.format(args.sector1[0], usedate)

        # Find LC file prefix and suffix number     
        lclist = glob.glob(os.path.join(spocPathPrefix, tp.spocdir,\
                        'light-curve','*lc.fits.gz'))
        exlc = os.path.basename(lclist[0])
        exlclist = exlc.split('-')
        tp.lcpre = '{0}-{1}-'.format(exlclist[0],exlclist[1])
        tp.lcnum = exlclist[3]

    else:
        # Multi sector
        print('Multi-Sector')
        tp.sector = -1
        tp.sector1 = args.sector1[0]
        tp.sector2 = args.sector2 # nargs with ? produces single item NOT list
        tp.tecdir = 'sector{0:d}-{1:d}'.format(args.sector1[0], args.sector2)
        if args.spocdir is None:
            print("For multi-sector you must specify SPOC data directory! EXITING!")
            exit()
        tp.spocdir = args.spocdir
        tp.tecfile = 'sector-{0:d}-{1:d}_{2}'.format(args.sector1[0],\
                            args.sector2, usedate)

        # Multi-sector does not have LC files
        # Find LC file prefix and suffix number from last sectyor  
        lclist = glob.glob(os.path.join(spocPathPrefix, 'sector-{0:03d}'.format(tp.sector2),\
                                                        'light-curve','*lc.fits.gz'))
        exlc = os.path.basename(lclist[0])
        exlclist = exlc.split('-')
        tp.lcpre = '{0}-{1}-'.format(exlclist[0],exlclist[1])
        tp.lcnum = exlclist[3]


    tp.toifile = 'FIXED-{0}'.format(usedate)

    # Find DV report prefix and suffix number                                 
    dvrlist = glob.glob(os.path.join(spocPathPrefix, tp.spocdir,\
                        'dv-reports','*_dvr.pdf'))
    exdvr = os.path.basename(dvrlist[0])
    exdvrlist = exdvr.split('-')
    tp.dvpre = '{0}-'.format(exdvrlist[0])
    tp.dvnum = exdvrlist[4][0:5]

    print('Sector = {0:d}'.format(tp.sector))
    print('Sector1 = {0:d}'.format(tp.sector1))
    print('Sector2 = {0:d}'.format(tp.sector2))
    print('TEC Dir {0}'.format(tp.tecdir))
    print('SPOC Dir {0}'.format(tp.spocdir))
    print('TEC Filename prefix {0}'.format(tp.tecfile))
    print('TOI File {0}'.format(tp.toifile))
    print('LC File Prefix {0}'.format(tp.lcpre))
    print('LC File Number {0}'.format(tp.lcnum))
    print('DV File Prefix {0}'.format(tp.dvpre))
    print('DV File Number {0}'.format(tp.dvnum))
