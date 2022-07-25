#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:35:04 2018
This gathers information for every detection you want to analyze
with TEC and makes a convenient pickle file with all this information 
@author: Christopher J. Burke (MIT)
"""
import numpy as np
import glob
import gzip
#from xml.dom.minidom import parse, parseString
import xml.etree.cElementTree as ET
import copy
import h5py



class tce_seed(object):
    """ Define a storage class that keeps all the info
        we want to keep for a TCE to propagate forward to vetting
    """
    # explictly state the name, datatype and shape for every
    #  class variable
    #  The names MUST exactly match the class variable names in the __init__
    store_names = ['epicId', 'planetNum', 'totPlanetNum', 'sourceId', \
                   'decDeg', 'raDeg', 'teff', 'teff_e', 'limbc', \
                   'feh', 'logg', 'rstar', 'rstar_e', \
                   'sector', 'ccd', 'camera', 'pixposvalid', \
                   'row', 'col', 'pmra', 'pmdec', \
                   'tmag', 'at_valid', 'at_snr', 'at_epochbtjd', \
                   'at_epochbtjd_e', 'at_rp', 'at_rp_e', 'at_imp', \
                   'at_dur', 'at_depth', 'at_depth_e', 'at_period', \
                   'at_period_e', 'at_rpDrstar', 'at_rpDrstar_e', 'at_aDrstar', \
                   'at_eqtemp', 'at_effflux', 'ghostcoreval', 'ghostcoresig', \
                   'ghosthalosig', 'ghosthaloval', 'modchi2', 'gof', \
                   'ntran', 'chi2', 'tce_epoch', 'mes', \
                   'maxsesinmes', 'tce_period', 'robstat', 'pulsedur', \
                   'trp_valid', 'trp_snr', 'trp_epochbtjd', 'trp_dur', \
                   'trp_depth', 'cent_oot_offset', 'cent_oot_offset_e', 'cent_tic_offset', \
                   'cent_tic_offset_e', 'oe_signif', 'data_start', 'data_end', \
                   'all_sectors', 'all_cadstart', 'all_cadend']
    store_types = ['i8', 'i4', 'i4', 'S80', \
                   'f8', 'f8', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'i4', 'i4', 'i4', 'i4', \
                   'f8', 'f8', 'f8', 'f8', \
                   'f8', 'i4', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'i4', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'i4', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'f8', 'f8', \
                   'f8', 'f8', 'i4', 'i4', \
                   'i4', 'i8', 'i8']
    store_shapes = [None, None, None, None, \
                    None, None, None, None, [4], \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    None, None, None, None, \
                    [100], [100], [100]]                    
    # Make the tuples that will define the numpy structured array
    # https://docs.scipy.org/doc/numpy/user/basics.rec.html
    sz = len(store_names)
    store_def_tuples = []
    for i in range(sz):
        if store_shapes[i] is not None:
            store_def_tuples.append((store_names[i], store_types[i], store_shapes[i]))
        else:
            store_def_tuples.append((store_names[i], store_types[i]))
    # Actually define the numpy structured/compound data type
    store_struct_numpy_dtype = np.dtype(store_def_tuples)

    def __init__(self):
        self.epicId = 0
        self.planetNum = 0
        self.totPlanetNum = 0
        self.sourceId = ''
        # Target information
        self.decDeg = 0.0
        self.raDeg = 0.0
        self.teff = 0.0
        self.teff_e = 0.0
        self.limbc = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)
        self.feh = 0.0
        self.logg = 0.0
        self.rstar = 0.0
        self.rstar_e = 0.0
        self.sector = 0
        self.ccd = 0
        self.camera = 0
        self.pixposvalid = 0
        self.row = 0.0
        self.col = 0.0
        self.pmra = 0.0
        self.pmdec = 0.0
        self.tmag = 0.0
        # All transit fit
        self.at_valid = 0
        self.at_snr = 0.0
        self.at_epochbtjd = 0.0
        self.at_epochbtjd_e = 0.0
        self.at_rp = 0.0
        self.at_rp_e = 0.0
        self.at_imp = 0.0
        self.at_dur = 0.0
        self.at_depth = 0.0
        self.at_depth_e = 0.0
        self.at_period = 0.0
        self.at_period_e = 0.0
        self.at_rpDrstar = 0.0
        self.at_rpDrstar_e = 0.0
        self.at_aDrstar = 0.0
        self.at_eqtemp = 0.0
        self.at_effflux = 0.0
        # Ghost diag
        self.ghostcoreval = 0.0
        self.ghostcoresig = 0.0
        self.ghosthalosig = 0.0
        self.ghosthaloval = 0.0
        # TCE Infor
        self.modchi2 = 0.0
        self.gof = 0.0
        self.ntran = 0
        self.chi2 = 0.0
        self.tce_epoch = 0.0
        self.mes = 0.0
        self.maxsesinmes = 0.0
        self.tce_period = 0.0
        self.robstat = 0.0
        self.pulsedur = 0.0
        # Trpzed fit
        self.trp_valid = 0
        self.trp_snr = 0.0
        self.trp_epochbtjd = 0.0
        self.trp_dur = 0.0
        self.trp_depth = 0.0
        # Centroid Fits
        self.cent_oot_offset = 0.0
        self.cent_oot_offset_e = 0.0
        self.cent_tic_offset = 0.0
        self.cent_tic_offset_e = 0.0
        # Odd/Even
        self.oe_signif = 0.0
        # Unit of Work data star and and
        self.data_start = 0
        self.data_end = 0
        # list of all sectors and start and end cadences numbers 
        self.all_sectors = np.array([-1]*100, dtype=int)
        self.all_cadstart = np.array([-1]*100, dtype=int)
        self.all_cadend = np.array([-1]*100, dtype=int)
        
    def store_objlist_as_hd5f(self, objlist, fileName):
        # save the class structure into hd5
        # objlist is a list of the tce_seed objects
        # First create the array of numpy structered arrays
        np_dset = np.ndarray(len(objlist), dtype=self.store_struct_numpy_dtype)
        # Convert the class variables into the numpy structured dtype
        for i, curobj in enumerate(objlist):
            for j in range(len(self.store_names)):
                np_dset[i][self.store_names[j]] = getattr(curobj, self.store_names[j])
        # Data set should be all loaded ready to write out
        fp = h5py.File(fileName, 'w')
        h5f_dset = fp.create_dataset('dset', shape=(len(objlist),), dtype=self.store_struct_numpy_dtype)
        h5f_dset[:] = np_dset
        fp.close()
        
    def fill_objlist_from_hd5f(self, fileName):
        fp = h5py.File(fileName, 'r')
        np_dset = np.array(fp['dset'])
        # Start with empty list
        all_objs = []
        # iterate through the numpy structured array and save to objects
        for i in range(len(np_dset)):
            tmp = tce_seed()
            for j in range(len(self.store_names)):
                setattr(tmp, self.store_names[j], np_dset[i][self.store_names[j]])
            # Append object to list
            all_objs.append(tmp)
        return all_objs
        
        
if __name__ == "__main__":
    tceSeedOutFile = 'sector-53_20220724_tce.h5'
    headXMLPath = '/pdo/spoc-data/sector-053/dv-results/'
    # Namespace there is extra junk prepended to tags
    #  This is supposed to make it easier to use 
    ns = {'ns': 'http://www.nasa.gov/2018/TESS/DV'}
    
    # Get list of XML files 
    fileList = glob.glob(headXMLPath + '*_dvr.xml*')
    # Gather data for each TCE
    all_tces = []
    for i in range(len(fileList)):
        if np.mod(i,20)==0:
            print("Parsed {0:d} of {1:d} Targs w/TCEs: {2:d}".format(i, len(fileList), len(all_tces)))

        
        # open and read gzipped xml file locally files are gzipped
        #  but they are not gzip from MAST deal with both situations
        # parse xml file content
        try:
            infile = gzip.open( fileList[i] )
            tree = ET.parse(infile)
        except:
            infile = open(fileList[i])
            tree = ET.parse(infile)

        #dom = parseString( content )
        root = tree.getroot()
        # Number of candidates
        try:
            ticId = int(root.get('ticId'))
            nCand = int(root.get('planetCandidateCount'))
        except:
            print('ERROR READING FILE {0}'.format(fileList[i]))
        # instantiate tce_seed class
        targdata = tce_seed()
        # Fill in the target specific information
        targdata.epicId = ticId
        targdata.totPlanetNum = nCand
        targdata.data_start = int(root.get('startCadence'))
        targdata.data_end = int(root.get('endCadence'))
        targdata.sourceId = 'SPOC'
        targdata.decDeg = float(root.find('ns:decDegrees', ns).attrib['value'])
        targdata.raDeg = float(root.find('ns:raDegrees', ns).attrib['value'])
        targdata.teff = float(root.find('ns:effectiveTemp', ns).attrib['value'])
        targdata.teff_e = float(root.find('ns:effectiveTemp', ns).attrib['uncertainty'])
        targdata.limbc[0] = float(root.find('ns:limbDarkeningModel', ns).attrib['coefficient1'])
        targdata.limbc[1] = float(root.find('ns:limbDarkeningModel', ns).attrib['coefficient2'])
        targdata.limbc[2] = float(root.find('ns:limbDarkeningModel', ns).attrib['coefficient3'])
        targdata.limbc[3] = float(root.find('ns:limbDarkeningModel', ns).attrib['coefficient4'])
        targdata.feh = float(root.find('ns:log10Metallicity', ns).attrib['value'])
        targdata.logg = float(root.find('ns:log10Metallicity', ns).attrib['value'])
        targdata.pmra = float(root.find('ns:pmRa', ns).attrib['value'])
        targdata.pmdec = float(root.find('ns:pmDec', ns).attrib['value'])
        targdata.rstar = float(root.find('ns:radius', ns).attrib['value'])
        targdata.rstar_e = float(root.find('ns:radius', ns).attrib['uncertainty'])
        targdata.tmag = float(root.find('ns:tessMag', ns).attrib['value'])
        # Get the camera position information from the first availalbe diff iage analysis
        tmp = root.find("ns:planetResults[@planetNumber='1']/ns:differenceImageResults", ns)
        targdata.sector = int(tmp.get('sector'))
        targdata.ccd = int(tmp[0].get('ccdNumber'))
        targdata.camera = int(tmp[0].get('cameraNumber'))
        tmp2 = tmp.find("ns:ticReferenceCentroid", ns)
        targdata.row = float(tmp2[0].get('value'))
        targdata.col = float(tmp2[1].get('value'))
        if np.isfinite(targdata.row) and np.isfinite(targdata.col) and (targdata.row > 0.0) and (targdata.col > 0.0):
            targdata.pixposvalid = 1
            
        # Get the cadence start and end for all sectors with data
        tmpall = root.findall("ns:planetResults[@planetNumber='1']/ns:differenceImageResults", ns)
        for idx, itmp in enumerate(tmpall):
            targdata.all_sectors[idx] =  int(itmp.get('sector'))
            targdata.all_cadstart[idx] = int(itmp.get('startCadence'))
            targdata.all_cadend[idx] = int(itmp.get('endCadence'))
        # Double check that the number of sectors agrees with sectorsObserved
        secobs = root.get('sectorsObserved')
        secsum = 0
        for x in secobs:
            secsum = secsum + int(x)
        if not len(np.where(targdata.all_sectors>=0)[0]) == secsum:
            print('Sector avail mismatch! {0:d} {1:d}'.format(i, targdata.epicId))
        
        
        for j in range(nCand):
            # Copy the target specific information to a new tce_seed class
            tcedata = copy.deepcopy(targdata)
            pres = root.find("ns:planetResults[@planetNumber='{0:d}']".format(j+1), ns)
            tcedata.planetNum = j+1
            # Get all transit fit information
            atfit = pres.find('ns:allTransitsFit', ns)
            tcedata.at_valid = int(atfit.get('fullConvergence') == 'true')
            tcedata.at_snr = float(atfit.get('modelFitSnr'))
            tcedata.at_epochbtjd = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='transitEpochBtjd']", ns).get('value'))
            tcedata.at_epochbtjd_e = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='transitEpochBtjd']", ns).get('uncertainty'))
            tcedata.at_rp = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='planetRadiusEarthRadii']", ns).get('value'))
            tcedata.at_rp_e = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='planetRadiusEarthRadii']", ns).get('uncertainty'))
            tcedata.at_imp = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='minImpactParameter']", ns).get('value'))
            tcedata.at_dur = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='transitDurationHours']", ns).get('value'))
            tcedata.at_depth = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='transitDepthPpm']", ns).get('value'))
            tcedata.at_period = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='orbitalPeriodDays']", ns).get('value'))
            tcedata.at_period_e = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='orbitalPeriodDays']", ns).get('uncertainty'))
            tcedata.at_rpDrstar = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='ratioPlanetRadiusToStarRadius']", ns).get('value'))
            tcedata.at_rpDrstar_e = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='ratioPlanetRadiusToStarRadius']", ns).get('uncertainty'))
            tcedata.at_aDrstar = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='ratioSemiMajorAxisToStarRadius']", ns).get('value'))
            tcedata.at_eqtemp = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='equilibriumTempKelvin']", ns).get('value'))
            tcedata.at_effflux = float(atfit.find("ns:modelParameters/ns:modelParameter[@name='effectiveStellarFlux']", ns).get('value'))
            ghostdat = pres.find('ns:ghostDiagnosticResults', ns)
            tcedata.ghostcoreval = float(ghostdat.find('ns:coreApertureCorrelationStatistic', ns).get('value'))
            tcedata.ghostcoresig = float(ghostdat.find('ns:coreApertureCorrelationStatistic', ns).get('significance'))
            tcedata.ghosthaloval = float(ghostdat.find('ns:haloApertureCorrelationStatistic', ns).get('value'))
            tcedata.ghosthalosig = float(ghostdat.find('ns:haloApertureCorrelationStatistic', ns).get('significance'))
            pcdat = pres.find('ns:planetCandidate', ns)
            tcedata.modchi2 = float(pcdat.get('modelChiSquare2'))/np.sqrt(float(pcdat.get('modelChiSquareDof2')))
            tcedata.ntran = int(pcdat.get('observedTransitCount'))
            tcedata.chi2 = float(pcdat.get('chiSquare2'))/np.sqrt(float(pcdat.get('chiSquareDof2')))
            tcedata.tce_epoch = float(pcdat.get('epochTjd'))
            tcedata.mes = float(pcdat.get('maxMultipleEventSigma'))
            tcedata.maxsesinmes = float(pcdat.get('maxSesInMes'))
            tcedata.tce_period = float(pcdat.get('orbitalPeriodInDays'))
            tcedata.robstat = float(pcdat.get('robustStatistic'))
            tcedata.pulsedur = float(pcdat.get('trialTransitPulseDurationInHours'))
            trpdat = pres.find('ns:trapezoidalFit', ns)
            tcedata.trp_valid = int(trpdat.get('fullConvergence') == 'true')
            tcedata.trp_snr = float(trpdat.get('modelFitSnr'))
            if tcedata.trp_valid == 1:                
                tcedata.trp_epochbtjd = float(trpdat.find("ns:modelParameters/ns:modelParameter[@name='transitEpochBtjd']", ns).get('value'))
                tcedata.trp_dur = float(trpdat.find("ns:modelParameters/ns:modelParameter[@name='transitDurationHours']", ns).get('value'))
                tcedata.trp_depth = float(trpdat.find("ns:modelParameters/ns:modelParameter[@name='transitDepthPpm']", ns).get('value'))

            centdat = pres.find('ns:centroidResults', ns)
            tcedata.cent_tic_offset = float(centdat.find("ns:differenceImageMotionResults/ns:msTicCentroidOffsets/ns:meanSkyOffset", ns).get('value'))
            tcedata.cent_tic_offset_e = float(centdat.find("ns:differenceImageMotionResults/ns:msTicCentroidOffsets/ns:meanSkyOffset", ns).get('uncertainty'))
            tcedata.cent_oot_offset = float(centdat.find("ns:differenceImageMotionResults/ns:msControlCentroidOffsets/ns:meanSkyOffset", ns).get('value'))
            tcedata.cent_oot_offset_e = float(centdat.find("ns:differenceImageMotionResults/ns:msControlCentroidOffsets/ns:meanSkyOffset", ns).get('uncertainty'))
            bindiscrim = pres.find('ns:binaryDiscriminationResults', ns)
            tcedata.oe_signif = float(bindiscrim.find("ns:oddEvenTransitDepthComparisonStatistic", ns).get('significance'))
            all_tces.append(tcedata)
#            for child in pcdat:
#                print(child.tag, child.attrib)
#            for child in pres:
#                print(child.tag, child.attrib)
#            print("hello World")

    print("Found {0:d} TCEs".format(len(all_tces)))
        # Write out hd5 file
    targdata.store_objlist_as_hd5f(all_tces, tceSeedOutFile)

    print("Wrote {0}".format(tceSeedOutFile))
