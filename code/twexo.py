"""
twexo.py - TESS Web Explode
  TESS Input Catalog web browser data exploder.
  WARNING: This routines opens tabs on your browser!
     This routine saves a local file for the html page!
     Always sanity check the target is correctly resolved!
     High proper motion targets may not be resolved correctly.
     
  Creates a local html page tailored for the target of interest
    and loads the page in your web browser with the links to the major
    catalogs to get information about the target.  Targets can be specified
    by TIC id, TOI number, common name, and ra/dec coordinates.
    
    
USAGE: 
        -to display command line arguments:
            python twexo.py -h
            
        -query by Name
            python twexo.py -n 'Pi Mensae'
            
        -For the impatient and the EXPLODE experience add -E
            python twexo.py -n 'Pi Mensae' -E
            
        -query by TOI number
            python twexo.py -toi 144.01
        -query by TIC number
            python twexo.py -t 261136679
        -query by Coordinates in decimal degrees
            python twexo.py -c 84.291198 -80.469143
       
      
AUTHORS: Christopher J. Burke (MIT)
 Testing and Advice from Susan Mullally (STScI) and Jennifer Burt (MIT)

VERSION: 0.5

NOTES: This routine opens tabs on your browser!
    This routine saves a local file for the html!
    
DEPENDENCIES:
    python 3+
    astropy
    astroquery
    numpy
    
SPECIAL THANKS TO:
    Includes code from the python MAST query examples 
    https://mast.stsci.edu/api/v0/pyex.html
    Brett Morris (UW) for help with the name resolving from another the tess-point project

"""

import webbrowser
import numpy as np
import os
import argparse
from astropy.coordinates import SkyCoord, ICRS
from astropy.time import Time, TimezoneInfo
from astroquery.gaia import Gaia
import astropy.units as u
import sys
import datetime
import json
try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
    from urllib.parse import urlencode as dict_urlencode
    from urllib.request import urlopen
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve
    from urllib import urlencode as dict_urlencode
    from urllib import urlopen
try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib
    
import csv
import time
# Look for optional tess-point for predicting which TESS sectors data is avail
foo_have_tp = True
try:
    from tess_stars2px import tess_stars2px_function_entry
except:
    foo_have_tp = False


def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list

## [Mast Query]
def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

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
## [Mast Query]

def safe_char(input):
    """ https://stackoverflow.com/questions/5843518/remove-all-special-characters-punctuation-and-spaces-from-string
    """
    if input:

        import re
        input = re.sub('\W+', '', input)

    return input


if __name__ == '__main__':
    # Parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--ticId", type=int, \
                        help="TIC Id [int] ")
    parser.add_argument("-c", "--coord", type=float, nargs=2, \
                        help="RA and Dec of target [deg]")
    parser.add_argument("-n", "--name", nargs=1, type=str, \
                        help="Search for a target by resolving its name with SESAME")
    parser.add_argument("-toi", "--toi", type=float, \
                        help="TOI number of target")
    parser.add_argument("-E", "--explode", action='store_true', \
                        help="Pre-load all URLs into tabs of browser rather than just the link page")
    parser.add_argument("-nw", "--noweb", action='store_true', \
                        help="Do not open in web.  Only generate html that plays well with wkhtmltopdf conversion tool")
    parser.add_argument("-of", "--outfile", nargs=1, type=str,\
                        help="Override default html output filename [default (without this argument) is based on query].  DO NOT include .html in this argument")
                        
    args = parser.parse_args()

# DEBUG BLOCK for hard coding input parameters and testing
 #   class test_arg:
 #       def __init__(self):
#            self.ticId = 366443426
#           self.ticId = None
#            self.coord = None
#            self.name = None
#            self.name = 'Pi Mensae'
#            self.toi = None
#            self.explode = False
#            self.toi = 383.01
#    args = test_arg()

    #Search radius for reporting nearby TIC targets
    SEARCH_RAD = 60.0 # arcsec
    # At least one Mode -t -c -n -toi must have been specified
    if (args.ticId is None) and (args.coord is None) and (args.name is None) and (args.toi is None):
        print('You must specify one and only one mode -t, -c, -n, -toi')
        print('`python twex.py -h\' for help')
        sys.exit(1)

    useTIC = 0
    # If TIC specified assign it to useTIC
    # and query MAST for RA and Dec and other catalog identifiers
    if args.ticId is not None:
        useTIC = int(args.ticId)
        starTics = np.array([useTIC], dtype=np.int32)
        ticStringList = ['{0:d}'.format(x) for x in starTics]    
        # Setup mast query
        request = {'service':'Mast.Catalogs.Filtered.Tic', \
           'params':{'columns':'*', 'filters':[{ \
                    'paramName':'ID', 'values':ticStringList}]}, \
            'format':'json', 'removenullcolumns':True}
        headers, outString = mastQuery(request)
        outObject = json.loads(outString)
        if len(outObject['data']) > 0:
            oo = outObject['data'][0]
            starRa = outObject['data'][0]['ra']
            starDec = outObject['data'][0]['dec']
            star2mass = outObject['data'][0]['TWOMASS']
            stargaia = outObject['data'][0]['GAIA']
            if 'pmRA' in oo:
                starPmRa = outObject['data'][0]['pmRA']
            else:
                starPmRa = 0.0
            if 'pmDEC' in oo:
                starPmDec = outObject['data'][0]['pmDEC']
            else:
                starPmDec = 0.0
            starPmTot = np.sqrt(starPmRa*starPmRa + starPmDec*starPmDec)
            if 'e_pmRA' in oo:
                starPmRaE = outObject['data'][0]['e_pmRA']
                starPmDecE = outObject['data'][0]['e_pmDEC']
                starPmTotE = np.sqrt(starPmRaE*starPmRaE + starPmDecE*starPmDecE)
            else:
                starPmRaE = 100.0
                starPmDecE = 100.0
                starPmTotE = 1000.0
            if 'Teff' in oo:
                starTeff = outObject['data'][0]['Teff']
            else:
                starTeff = 0.0
            if 'e_Teff' in oo:
                starTeffE = outObject['data'][0]['e_Teff']
            else:
                starTeffE = 0.0
            if 'logg' in oo:
                starLogg = outObject['data'][0]['logg']
            else:
                starLogg = 0.0
            if 'e_logg' in oo:
                starLoggE = outObject['data'][0]['e_logg']
            else:
                starLoggE = 0.0
            if 'rad' in oo:
                starRad = outObject['data'][0]['rad']
            else:
                starRad = 0.0
            if 'e_rad' in oo:
                starRadE = outObject['data'][0]['e_rad']
            else:
                starRadE = 0.0
            if 'mass' in oo:
                starMass = outObject['data'][0]['mass']
            else:
                starMass = 0.0
            if 'e_mass' in oo:
                starMassE = outObject['data'][0]['e_mass']
            else:
                starMassE = 0.0
            starTmag = outObject['data'][0]['Tmag']
            
            #print(outObject['data'][0])
        else:
            print('MAST search of TIC: {0} did not return an entry'.format(args.ticId))
            sys.exit(1)
        
    # TOI specified get toi -> tic mapping from table at exofop
    if args.toi is not None:
        toiURL = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv'
        tmp = []
        with urlopen(toiURL) as response:
            for line in response:
                line = line.decode('utf-8')
                for x in csv.reader([line]):
                    tmp.append("|".join(x))

        dtypeseq = ['i4','f8']
        dtypeseq.extend(['i4']*7)
        dtypeseq.extend(['U3','U3'])
        dtypeseq.extend(['f8','f8'])
        dtypeseq.extend(['U40','i4','U40','U40','U40'])
        dtypeseq.extend(['f8']*27)
        dtypeseq.extend(['U80']*4)
        dataBlock = np.genfromtxt(tmp, \
                                  dtype=dtypeseq, delimiter='|',skip_header=1)
        exoTic = dataBlock['f0']
        exoToi = dataBlock['f1']
        idx = np.where(exoToi == args.toi)[0]
        if len(idx) == 0:
            print('Could not find TOI {0} at ExoFop'.format(args.toi))
            sys.exit(1)
        useTIC = int(exoTic[idx])
        # Now that we have the TIC do a MAST query to get TIC values
        starTics = np.array([useTIC], dtype=np.int32)
        ticStringList = ['{0:d}'.format(x) for x in starTics]    
        # Setup mast query
        request = {'service':'Mast.Catalogs.Filtered.Tic', \
           'params':{'columns':'*', 'filters':[{ \
                    'paramName':'ID', 'values':ticStringList}]}, \
            'format':'json', 'removenullcolumns':True}
        headers, outString = mastQuery(request)
        outObject = json.loads(outString)
        oo = outObject['data'][0]
        starRa = outObject['data'][0]['ra']
        starDec = outObject['data'][0]['dec']
        star2mass = outObject['data'][0]['TWOMASS']
        stargaia = outObject['data'][0]['GAIA']
        if 'pmRA' in oo:
            starPmRa = outObject['data'][0]['pmRA']
        else:
            starPmRa = 0.0
        if 'pmDEC' in oo:
            starPmDec = outObject['data'][0]['pmDEC']
        else:
            starPmDec = 0.0
        starPmTot = np.sqrt(starPmRa*starPmRa + starPmDec*starPmDec)
        if 'e_pmRA' in oo:
            starPmRaE = outObject['data'][0]['e_pmRA']
            starPmDecE = outObject['data'][0]['e_pmDEC']
            starPmTotE = np.sqrt(starPmRaE*starPmRaE + starPmDecE*starPmDecE)
        else:
            starPmRaE = 100.0
            starPmDecE = 100.0
            starPmTotE = 1000.0

        if 'Teff' in oo:
            starTeff = outObject['data'][0]['Teff']
        else:
            starTeff = 0.0
        if 'e_Teff' in oo:
            starTeffE = outObject['data'][0]['e_Teff']
        else:
            starTeffE = 0.0
        if 'logg' in oo:
            starLogg = outObject['data'][0]['logg']
        else:
            starLogg = 0.0
        if 'e_logg' in oo:
            starLoggE = outObject['data'][0]['e_logg']
        else:
            starLoggE = 0.0
        if 'rad' in oo:
            starRad = outObject['data'][0]['rad']
        else:
            starRad = 0.0
        if 'e_rad' in oo:
            starRadE = outObject['data'][0]['e_rad']
        else:
            starRadE = 0.0
        if 'mass' in oo:
            starMass = outObject['data'][0]['mass']
        else:
            starMass = 0.0
        if 'e_mass' in oo:
            starMassE = outObject['data'][0]['e_mass']
        else:
            starMassE = 0.0
        starTmag = outObject['data'][0]['Tmag']

        
        #print('TOI: {0} is associated with TIC: {1} at ExoFop'.format(args.toi, useTIC))

    # Do name resolving to get ra and dec
    if args.name is not None:
        # Name resolve  in try except  for detecting problem
        try:
            coordinate = SkyCoord.from_name(args.name[0])
            print("Coordinates for {0}: ({1}, {2})"
              .format(args.name[0], coordinate.ra.degree,
                      coordinate.dec.degree))
            starRa = coordinate.ra.degree
            starDec = coordinate.dec.degree
        except:
            print("Could not resolve: {0}".format(args.name[0]))
            sys.exit(1)


    if args.coord is not None:
        starRa = args.coord[0]
        starDec = args.coord[1]
    # Do a MAST cone search around the ra and dec to get the nearby stars
    #  for all runs and to find the closest  TIC for the coordinate input 
    #   and name query
    # Protect against missing TICs
    try:
        # Do cone search around this position
        startTime = time.time()
        request = {'service':'Mast.Catalogs.Filtered.Tic.Position', \
                   'params':{'columns':'c.*', \
                             'filters':[ \
                                        {'paramName':'Tmag',\
                                         'values':[{'min':0, 'max':20.0}]}], \
                             'ra':'{:10.5f}'.format(starRa),\
                             'dec':'{:10.5f}'.format(starDec),\
                             'radius':'{:10.7f}'.format(SEARCH_RAD/3600.0) \
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
         
        ticList = np.array([x['ID'] for x in outObject['data']], dtype=np.int32)
        ticRas = np.array([x['ra'] for x in outObject['data']], dtype=np.float)
        ticDecs = np.array([x['dec'] for x in outObject['data']], dtype=np.float)
        ticTeffs = np.array([x['Teff'] for x in outObject['data']], dtype=np.float)
        ticTeffEs = np.array([x['e_Teff'] for x in outObject['data']], dtype=np.float)
        ticLoggs = np.array([x['logg'] for x in outObject['data']], dtype=np.float)
        ticLoggEs = np.array([x['e_logg'] for x in outObject['data']], dtype=np.float)
        ticRads = np.array([x['rad'] for x in outObject['data']], dtype=np.float)
        ticRadEs = np.array([x['e_rad'] for x in outObject['data']], dtype=np.float)
        ticMasss = np.array([x['mass'] for x in outObject['data']], dtype=np.float)
        ticMassEs = np.array([x['e_mass'] for x in outObject['data']], dtype=np.float)
        ticTmags = np.array([x['Tmag'] for x in outObject['data']], dtype=np.float)

    except:
        # Cone search failed
        print('TIC Search Fail at position {0} {1}'.format(starRa, starDec))
        sys.exit(1)

    # Check if no targets were returned in cone search
    if len(ticList) == 0:
        print('MAST cone search around TIC returned no objects')
        sys.exit(1)
    
    targCoord = SkyCoord(starRa, starDec, unit='deg')
    ticCoords = SkyCoord(ticRas, ticDecs, unit='deg')
    seps = targCoord.separation(ticCoords)
    ia = np.argsort(seps)
    ticList, ticRas, ticDes, seps, ticTeffs, ticTeffEs, ticLoggs, ticLoggEs, \
        ticRads, ticRadEs, ticMasss, ticMassEs, ticTmags = idx_filter(ia , \
         ticList, ticRas, ticDecs, seps, ticTeffs, ticTeffEs, ticLoggs, ticLoggEs, \
        ticRads, ticRadEs, ticMasss, ticMassEs, ticTmags)

    # If we hadn't defined which target is useTIC yet do it now
    if useTIC == 0:
        useTIC = ticList[0]
        starTics = np.array([useTIC], dtype=np.int32)
        ticStringList = ['{0:d}'.format(x) for x in starTics]    
        # Setup mast query
        request = {'service':'Mast.Catalogs.Filtered.Tic', \
           'params':{'columns':'*', 'filters':[{ \
                    'paramName':'ID', 'values':ticStringList}]}, \
            'format':'json', 'removenullcolumns':True}
        headers, outString = mastQuery(request)
        outObject = json.loads(outString)
        #print(outObject['data'][0])
        oo = outObject['data'][0]
        starRa = outObject['data'][0]['ra']
        starDec = outObject['data'][0]['dec']
        star2mass = outObject['data'][0]['TWOMASS']
        stargaia = outObject['data'][0]['GAIA']
        if 'pmRA' in oo:
            starPmRa = outObject['data'][0]['pmRA']
        else:
            starPmRa = 0.0
        if 'pmDEC' in oo:
            starPmDec = outObject['data'][0]['pmDEC']
        else:
            starPmDec = 0.0
        starPmTot = np.sqrt(starPmRa*starPmRa + starPmDec*starPmDec)
        if 'e_pmRA' in oo:
            starPmRaE = outObject['data'][0]['e_pmRA']
            starPmDecE = outObject['data'][0]['e_pmDEC']
            starPmTotE = np.sqrt(starPmRaE*starPmRaE + starPmDecE*starPmDecE)
        else:
            starPmRaE = 100.0
            starPmDecE = 100.0
            starPmTotE = 1000.0
        if 'Teff' in oo:
            starTeff = outObject['data'][0]['Teff']
        else:
            starTeff = 0.0
        if 'e_Teff' in oo:
            starTeffE = outObject['data'][0]['e_Teff']
        else:
            starTeffE = 0.0
        if 'logg' in oo:
            starLogg = outObject['data'][0]['logg']
        else:
            starLogg = 0.0
        if 'e_logg' in oo:
            starLoggE = outObject['data'][0]['e_logg']
        else:
            starLoggE = 0.0
        if 'rad' in oo:
            starRad = outObject['data'][0]['rad']
        else:
            starRad = 0.0
        if 'e_rad' in oo:
            starRadE = outObject['data'][0]['e_rad']
        else:
            starRadE = 0.0
        if 'mass' in oo:
            starMass = outObject['data'][0]['mass']
        else:
            starMass = 0.0
        if 'e_mass' in oo:
            starMassE = outObject['data'][0]['e_mass']
        else:
            starMassE = 0.0
        starTmag = outObject['data'][0]['Tmag']

    # If we weren't given TOI number check exofop to see if it has TOIs
    if args.toi is None:
        toiURL = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv'
        tmp = []
        response = urlopen(toiURL)
        for line in response:
            line = line.decode('utf-8')
            for x in csv.reader([line]):
                tmp.append("|".join(x))

        dtypeseq = ['i4','f8']
        dtypeseq.extend(['i4']*7)
        dtypeseq.extend(['U3','U3'])
        dtypeseq.extend(['f8','f8'])
        dtypeseq.extend(['U40','i4','U40','U40','U40'])
        dtypeseq.extend(['f8']*27)
        dtypeseq.extend(['U80']*4)
        dataBlock = np.genfromtxt(tmp, \
                                  dtype=dtypeseq, delimiter='|',skip_header=1)
        exoTic = dataBlock['f0']
        exoToi = dataBlock['f1']
        idx = np.where(exoTic == useTIC)[0]
        if len(idx) == 0:
            toihdr = 'No TOIs associated with TIC at ExoFop'
            toicheckstr = '<br>'
        else:
            toihdr = 'TIC Hosts TOIs'
            toilist = []
            for j in idx:
                toilist.append('{0:.2f}<br>'.format(exoToi[j]))     
            toicheckstr = ' '.join(toilist)
    else:
        toihdr = 'Using TOI {:.2f}'.format(args.toi)
        toicheckstr = ''
        
    # FORM THE URLS
    gaiaURLPart1 = 'https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery={"service":"GAIADR2","inputText":"'
    gaiaURLPart2 = '{0:f} {1:f} r{2:f}",'.format(starRa,starDec,1.0/60.0)
    gaiaURLPart3 = '"paramsService":"Mast.Catalogs.GaiaDR2.Cone","title":"Gaia (DR2)","columns":"*"}'
    gaiaURL = gaiaURLPart1 + gaiaURLPart2 + gaiaURLPart3


    # TESScut Target pixel file
    tcutURLPart1 = 'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'
    tcutURLPart2 = 'ra={0}&dec={1}'.format(starRa, starDec)
    tcutURLPart3 = '&y=15&x=15&units=px&sector=All'
    tcutURL = tcutURLPart1 + tcutURLPart2 + tcutURLPart3

        
    # ESO archive
    DATE_STR = datetime.date.today().strftime('%Y-%m-%d')
    esoURLPart1 = 'https://archive.eso.org/scienceportal/home?data_release_date=*:{0}&'.format(DATE_STR)
    esoURLPart2 = 'pos={0},{1}'.format(starRa,starDec)
#    esoURLPart3 = '&r=0.008333&sort=dist,-fov,-obs_date&s=P%2fDSS2%2fcolor&f=0.064387&fc=84.485552,-80.451816&cs=J2000&av=true&ac=false&c=8,9,10,11,12,13,14,15,16,17,18&mt=true&dts=true'
    esoURLPart3 = '&r=0.008333&sort=dist,-fov,-obs_date&s=P%2fDSS2%2fcolor&f=0.064387&cs=J2000&av=true&ac=false&c=8,9,10,11,12,13,14,15,16,17,18&mt=true&dts=true'

    esoURL = esoURLPart1 + esoURLPart2 + esoURLPart3

        
    # IRSA finderchart
    irsaURLPart1 = 'https://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?locstr='
    #irsaURLPart1 = 'https://irsa.ipac.caltech.edu/applications/finderchart/?__action=table.search&request=%7B%22startIdx%22%3A0%2C%22pageSize%22%3A100%2C%22id%22%3A%22QueryFinderChartWeb%22%2C%22tbl_id%22%3A%22results%22%2C%22UserTargetWorldPt%22%3A%22126.61572%3B10.08046%3BEQ_J2000%3B'
    irsaURLPart2 = urlencode('2MASS J{0}'.format(star2mass))
#    irsaURLPart2 = '2MASS J{0}'.format(star2mass)
    irsaURLPart3 = '&mode=getResult&subsetsize=4.0&searchCatalog=no&survey='
    irsaURLPart4 = 'DSS,SDSS,2MASS'
#    irsaURLPart4 = urlencode('DSS,SDSS,2MASS')
    #irsaURLPart3 = '%3Bned%22%2C%22imageSizeAndUnit%22%3A%220.042777777777777776%22%2C%22thumbnail_size%22%3A%22256%22%2C%22selectImage%22%3A%22dss%2Csdss%2C2mass%22%2C%22searchCatalog%22%3A%22no%22%2C%22ckgDSS%22%3A%22dss1Blue%2Cdss1Red%2Cdss2Blue%2Cdss2Red%2Cdss2IR%22%2C%22ckgSDSS%22%3A%22u%2Cg%2Cr%2Cz%22%2C%22ckg2MASS%22%3A%22j%2Ch%2Ck%22%2C%22imageSearchOptions%22%3A%22closed%22%2C%22META_INFO%22%3A%7B%22title%22%3A%22QueryFinderChartWeb%22%2C%22tbl_id%22%3A%22results%22%7D%7D&options=%7B%22tbl_group%22%3A%22results%22%2C%22removable%22%3Afalse%2C%22showTitle%22%3Afalse%2C%22pageSize%22%3A100%7D'
    irsaURL = irsaURLPart1 + irsaURLPart2 + irsaURLPart3 + irsaURLPart4
    
    
    # MAST Data portal
    #mstURLPart1 = 'https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery=%7B%22service%22%3A%22CAOMDB%22%2C%22inputText%22%3A%22'
    #mstURLPart2 = 'TIC%20{0:d}'.format(useTIC)
    #mstURLPart3 = '%22%2C%22paramsService%22%3A%22Mast.Caom.Cone%22%2C%22columns%22%3A%22*%22%2C%22caomVersion%22%3Anull%7D'
    mstURLPart1 = 'https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery={"service":"CAOMFILTERED","inputText":[{"paramName":"obs_collection","niceName":"obs_collection","values":[],"valString":"TESS","isDate":false,"separator":";","freeText":"TESS","displayString":"TESS"}],"position":"'
    mstURLPart2 = '{0:f}, {1:f}, {2:f}",'.format(starRa,starDec,2.0/60.0/60.0)
    mstURLPart3 = '"paramsService":"Mast.Caom.Filtered.Position","columns":"*"}'    
    mstURL = mstURLPart1 + mstURLPart2 + mstURLPart3
    #print(mstURL)

    #Vizier search
    vizPOSTDict = {'-c': '2MASS J{0}'.format(star2mass)}
    vizURLPart1 = 'http://vizier.u-strasbg.fr/viz-bin/VizieR-4?-out.max=50&-out.form=HTML+Table&-out.add=_r&-out.add=_RAJ+_DEJ&outaddvalue=default&-sort=_r&-order=I&-oc.form=sexa&'
    vizURLPart2 = dict_urlencode(vizPOSTDict)
    vizURLPart3 = '&-c.eq=J2000&-c.r=+30&-c.u=arcsec&-c.geom=r'
    vizURL = vizURLPart1 + vizURLPart2 + vizURLPart3
    #print(vizURL)

    # exofop 
    exofopURL = 'https://exofop.ipac.caltech.edu/tess/target.php?id={0}'.format(useTIC)
    
    # Simbad page
    simbadPOSTDict = {'Ident': '2MASS J{0}'.format(star2mass)}
    simbadURL = 'http://simbad.u-strasbg.fr/simbad/sim-basic?{0}'.format(dict_urlencode(simbadPOSTDict))
    #print(simbadURL)

    # Prepare strings for the HTML page
    DATE_STR = datetime.date.today().strftime('%Y%m%d')

    if args.toi is not None:
        path = os.path.abspath('twexo_temp_{0}_{1:.2f}.html'.format(DATE_STR,args.toi))
        ModeHeader = 'TOI: {0} is associated with TIC: {1} at ExoFop'.format(args.toi, useTIC)
    if args.ticId is not None:
        path = os.path.abspath('twexo_temp_{0}_{1:016d}.html'.format(DATE_STR,args.ticId))
        ModeHeader = 'Using TIC: {0}'.format(useTIC)
    if args.name is not None:
        # strip white space and only keep regular characters to protect against
        #  special characters that would be bad to put into filenames also clip lenght

        useName = safe_char(args.name[0])
        if (len(useName)>30):
            useName = useName[0:11]
        path = os.path.abspath('twexo_temp_{0}_{1}.html'.format(DATE_STR,useName))
        ModeHeader = 'From {0} Using TIC: {1}'.format(args.name[0], useTIC)
    if args.coord is not None:
        path = os.path.abspath('twexo_temp_{0}_{1:016d}.html'.format(DATE_STR,useTIC))
        ModeHeader = 'From input coordinates Using TIC: {0}'.format(useTIC)
    if args.outfile is not None:
        path = os.path.relpath(args.outfile[0]+'.html')
    print(path)
    ticPos = 'Using TIC catalog position {0:.6f} {1:.6f} [J2000.0; epoch 2000.0]'.format(starRa, starDec)
    twoMassId = '2MASS J{0} From TIC'.format(star2mass)
    closeTICN = '{0} TIC entries within {1} arcsec of target {2}'.format(len(ticList)-1, SEARCH_RAD, useTIC)
    neighTIC = []
    for i, curTIC in enumerate(ticList[1:21]):
        if np.isfinite(ticTeffs[i+1]):
            curTeff = ticTeffs[i+1]
        else:
            curTeff = 0.0
        if np.isfinite(ticLoggs[i+1]):
            curLogg = ticLoggs[i+1]
        else:
            curLogg = 0.0
        if np.isfinite(ticRads[i+1]):
            curRad = ticRads[i+1]
        else:
            curRad = 0.0
        neighTIC.append('{0:d} Sep [arcsec]: {1:8.3f} Tmag: {2:5.2f} Teff: {3:6.1f} Logg: {4:5.2f} Rs[Rsun]: {5:6.2f}<br>'.format(\
                        curTIC, seps[i+1].arcsecond, ticTmags[i+1], curTeff, \
                        curLogg, curRad))     
    neighTICStr = ' '.join(neighTIC)

    # Since we got our positions from the TIC (which are in epoch 2000.0)
    # propagate the proper motion to find coordinates in 2015.5 which is gaia
    # First sanity check pm data from TIC
    pmOK = True
    if (not np.isfinite(starPmTot)) or (not np.isfinite(starPmTotE)):
        pmOK = False
    if np.abs(starPmTot)/starPmTotE < 1.5:
        pmOK = False
    if pmOK:
        print('TIC has proper motion data.  Getting 2015.5 position for GAIA')
#        ctic = SkyCoord(ra=starRa * u.deg,
#                        dec=starDec * u.deg,
#                        pm_ra_cosdec=starPmRa * u.mas/u.yr,
#                        pm_dec=starPmDec * u.mas/u.yr,
#                        obstime=Time('J2000.0'))#,
#                        distance=1000.0*u.pc, # Use fake distance just for skycoord functionality
#                        radial_velocity = 0.0 * u.km/u.s) # Use fake rv "
        # convert position to GAIA DR2 2015.5
#        gaiaEpc = Time('J2015.5')
#        ctic_gaia_epoch = ctic.apply_space_motion(gaiaEpc)
#        gaiaPredRa = ctic_gaia_epoch.ra.degree
#        gaiaPredDec = ctic_gaia_epoch.dec.degree
        # Compare to manual way of adding 15.5 years of proper motion
        starPmRa_arcsecDYr = starPmRa / 1000.0
        starPmDec_arcsecDYr = starPmDec / 1000.0
        cosDec = np.cos(starDec * np.pi/180.0)
        starPmRa_degDYrWcosDec = starPmRa_arcsecDYr / 3600.0 / cosDec
        starPmDec_degDYr = starPmDec_arcsecDYr / 3600.0 
        starPmRaDel_degWcosDec = starPmRa_degDYrWcosDec * 15.5
        starPmDecDel_deg = starPmDec_degDYr * 15.5
        gaiaPredRa = starRa + starPmRaDel_degWcosDec
        gaiaPredDec = starDec + starPmDecDel_deg

    else:
        gaiaPredRa = starRa
        gaiaPredDec = starDec
    gaiaPos = 'Predicted GAIA position {0:.6f} {1:.6f} [J2000.0; epoch 2015.5]'.format(gaiaPredRa, gaiaPredDec)

    # We are go for GAIA Cone search
    ADSQL_Str = "SELECT \
       DISTANCE( POINT('ICRS', ra, dec),\
       POINT('ICRS', {0}, {1}) ) AS dist, phot_g_mean_mag, teff_val, teff_percentile_lower, \
       teff_percentile_upper, radius_val, radius_percentile_lower, radius_percentile_upper,\
       astrometric_gof_al, astrometric_excess_noise_sig, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, \
       a_g_val, e_bp_min_rp_val, parallax, \
        ra, dec , solution_id \
        from gaiadr2.gaia_source WHERE 1 = CONTAINS( POINT('ICRS', ra, dec), \
        CIRCLE('ICRS', {0}, {1}, 0.016666666666666666)) order by dist".format(gaiaPredRa, gaiaPredDec)
    #print(ADSQL_Str)
    job = Gaia.launch_job(ADSQL_Str)
    r = job.get_results()
    gaiaTeff = r['teff_val'][0]
    gaiaTeffE = ((r['teff_percentile_upper'][0] - gaiaTeff) + (gaiaTeff - r['teff_percentile_lower'][0]))/2.0
    gaiaRad = r['radius_val'][0]
    gaiaRadE = ((r['radius_percentile_upper'][0] - gaiaRad) + (gaiaRad - r['radius_percentile_lower'][0]))/2.0
    gaiaRpMag = r['phot_rp_mean_mag'][0]
    gaiaGMag = r['phot_g_mean_mag'][0]
    gaiaBpMag = r['phot_bp_mean_mag'][0]
    gaiapar = r['parallax'][0]
    gaiaag = float(r['a_g_val'][0])
    if not np.isfinite(gaiaag):
        gaiaag = 0.0        
    gaiaexclr = float(r['e_bp_min_rp_val'][0])
    if not np.isfinite(gaiaexclr):
        gaiaexclr = 0.0
    gaiaastrogof = r['astrometric_gof_al'][0]
    gaiaastrexsig = r['astrometric_excess_noise_sig'][0]
    #print(gaiaGMag, gaiapar, gaiaag, gaiaBpMag, gaiaRpMag, gaiaexclr)
    gaiaAbsGMag = gaiaGMag + 5.0 + 5.0*np.log10(gaiapar/1000.0)-gaiaag
    gaiaClro = gaiaBpMag - gaiaRpMag - gaiaexclr

    gaiaExStr = 'Other GAIA G: {0:5.2f} Bp: {1:5.2f} AbsG: {2:5.2f} (Bp-Rp)o: {3:5.3f} AstroGOF: {4:6.2f} AstroExNoiSig: {5:6.2f}<br>'.format(
            gaiaGMag, gaiaBpMag, gaiaAbsGMag, gaiaClro, gaiaastrogof, gaiaastrexsig)
#    print(r['teff_percentile_upper'])
#    print(r['teff_percentile_lower'])
    #print(r)

    #LOOK for DV reports for this target available at MAST
    # Setup mast query first step is to get obsid's for the target
    print('Query MAST for DV reports')
    dvStr = 'No DV Results for this target at MAST<br>'
    request = {'service':'Mast.Caom.Filtered.Position', \
       'params':{'columns':'*', \
                 'position':'{0:f}, {1:f}, {2:f}'.format(starRa,starDec,2.0/60.0/60.0), \
                 'filters':[
                         {'paramName':'project',
                         'values':['TESS']
                         },\
                         {'paramName':'dataproduct_type',
                          'values':['timeseries']
                          }]
                }, \
        'format':'json', 'removenullcolumns':True}
    headers, outString = mastQuery(request)
    outObject = json.loads(outString)
    if not outObject['status'] == 'EXECUTING':
        # Check if any objects returned
        if len(outObject['data']) > 0:
            # Now get list of obsids for the time series data
            #print(len(outObject['data']))
            obsidStr = []
            for i, oo in enumerate(outObject['data']):
                #print(i)
                #print(oo)
                #print(oo['obsid'])
                obsidStr.append('{0:d}'.format(oo['obsid']))
            # Have list of obs ids for this TIC
            # Request the data products associated with this obs id
            request = {'service':'Mast.Caom.Products', \
               'params':{'obsid':','.join(obsidStr), \
                         }, \
                'format':'json', 'removenullcolumns':True}
            headers, outString = mastQuery(request)
            outObject = json.loads(outString)
            #print(len(outObject['data']))
            dvLines = []
            for i, oo in enumerate(outObject['data']):
                #print(i)
                dataURI = oo['dataURI']
                suffix = dataURI[-3:]
                if suffix == 'pdf':
                    useURL = 'https://mast.stsci.edu/api/v0/download/file?uri=' + dataURI
                    dvLines.append('<a href="{0}" target="_blank">{0}</a><br>'.format(useURL))
            dvStr = ' '.join(dvLines)        
    
        else:
            dvStr = 'No DV Results for this target at MAST<br>'
    else:
        print('Timeout on the DV report search at MAST')
        dvStr = '***DV report search at MAST timed out.  Try Later***<br>'
    
    # Do an optional report that shows which sectors target is TESS observable
    #  One needs tess-point install for this to work
    tesspStr = '<br>'
    if foo_have_tp:
            outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, \
            outColPix, outRowPix, scinfo = tess_stars2px_function_entry(
                    useTIC, gaiaPredRa, gaiaPredDec)
            tpLines = ['<br>TESS Observe - Sec Cam Ccd  Col  Row<br>']
            for i, curTIC in enumerate(outID):
                tpLines.append('{0:d} {1:d} {2:d} {3:d} {4:.2f} {5:.2f}<br>'.format(outID[i], \
                               outSec[i],outCam[i], outCcd[i], outColPix[i], outRowPix[i]))
            tesspStr = ' '.join(tpLines)


    # Wkhtmltopdf is used to convert html to pdf
    #  On certain platforms URLs don't work unless they are one per line
    #  with args.noweb use a html line break
    URLBreak = '|'
    if args.noweb:
        URLBreak = '<br>'
    # We had all the HTML strings, put them in a dictionary that will
    #  will be used for replacement in HTML template below        
    htmlDict = {'neighTIC':neighTICStr, 'closeTICN':closeTICN, \
                'twoMassId':twoMassId, 'ticPos':ticPos, \
                'ModeHeader':ModeHeader, 'useTIC':useTIC, \
                'exofopURL':exofopURL, 'simbadURL':simbadURL, 'vizURL':vizURL, \
                'mstURL':mstURL, 'irsaURL':irsaURL, 'esoURL':esoURL, \
                'toicheckstr':toicheckstr, 'toihdr':toihdr, 'gaiaPos':gaiaPos,\
                'T_Tmag':starTmag, 'T_Teff':'{0:7.1f}&plusmn{1:5.1f}'.format(starTeff,starTeffE), \
                'T_Logg':'{0:4.2f}&plusmn{1:4.2f}'.format(starLogg, starLoggE), \
                'T_Rstar':'{0:6.2f}&plusmn{1:4.1f}'.format(starRad, starRadE), \
                'T_Mstar':'{0:5.2f}&plusmn{1:5.2f}'.format(starMass, starMassE), \
                'G_Rmag':'{0:6.3f}'.format(gaiaRpMag), \
                'G_Teff':'{0:7.1f}&plusmn{1:5.1f}'.format(gaiaTeff, gaiaTeffE),\
                'G_Rstar':'{0:6.2f}&plusmn{1:4.1f}'.format(gaiaRad, gaiaRadE), \
                'tesspStr':tesspStr, 'dvStr':dvStr, 'gaiaExStr':gaiaExStr, \
                'tcutURL':tcutURL, 'gaiaURL':gaiaURL, 'URLBreak':URLBreak}
    #print(htmlDict['mstURL'])
    #print('{mstURL}'.format(**htmlDict))
    #print('hello')
    #HTML TEMPLATE
    template = """
<html>
<head>
<title>TESS Web Exploder {useTIC}</title>
<style>
table, th, td {{
  border: 1px solid black;
  border-collapse: collapse;
}}
th, td {{
  padding: 5px;
}}
th {{
  text-align: left;
}}
</style>
</head>
<body>
<h1>{ModeHeader}</h1>
<h3>{ticPos}</h3>
<h3>{gaiaPos}</h3>
<h3>{twoMassId}</h3>
<h3>{closeTICN}</h3>
{neighTIC}
<h3>{toihdr}</h3>
{toicheckstr}
<h3>Target Parameters</h3>
<table style="width:100%">
  <tr>
    <th>Catalog</th>
    <th>Tmag/Rpmag</th> 
    <th>Teff</th>
    <th>Logg</th>
    <th>Rstar</th>
    <th>Mstar</th>
  </tr>
  <tr>
    <td>TIC</td>
    <td>{T_Tmag}</td> 
    <td>{T_Teff}</td>
    <td>{T_Logg}</td>
    <td>{T_Rstar}</td>
    <td>{T_Mstar}</td>
  </tr>
  <tr>
    <td>GAIA DR2</td>
    <td>{G_Rmag}</td> 
    <td>{G_Teff}</td>
    <td>...</td>
    <td>{G_Rstar}</td>
    <td>...</td>
  </tr>
</table>
{gaiaExStr}
<h2>Target Links</h2>
<a href="{exofopURL}" target="_blank">ExoFOP</a> {URLBreak}
<a href="{simbadURL}" target="_blank">Simbad</a> {URLBreak}
<a href="{vizURL}" target="_blank">Vizier</a> {URLBreak}
<a href='{mstURL}' target='_blank'>MAST TESS Data Holdings</a> {URLBreak}
<a href="{irsaURL}" target="_blank">IRSA Finderchart</a> {URLBreak}
<a href="{esoURL}" target="_blank">ESO Data Archive Holdings</a> {URLBreak}
<a href="{tcutURL}" target="_blank">TESScut TPF Download</a> {URLBreak}
<a href='{gaiaURL}' target='_blank'>GAIA DR2 60'' Cone Search @MAST</a> 
<h2>NASA Ames SPOC DV Results Available at MAST</h2>
{dvStr}
<br>
{tesspStr}
</body>
</html>
"""
    # Replace tags in HTML template with the strings in the htmlDict
    page_string = template.format(**htmlDict)
    url = 'file://' + path
    # Write the html to a local file
    with open(path, 'w') as f:
        f.write(page_string)

    if not args.noweb:
        # If the user requested a web explode load em all!
        if args.explode:
            webbrowser.open(gaiaURL, new=2)
            webbrowser.open(esoURL, new=2)
            webbrowser.open(irsaURL, new=2)
            webbrowser.open(mstURL, new=2)
            webbrowser.open(vizURL, new=2)
            webbrowser.open(simbadURL, new=2)
            webbrowser.open(exofopURL, new=2)  
    
    
        # Load the webpage
        webbrowser.open(url)
    
        
    