# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 14:46:59 2017
Some personal numpy array filtering and finding intersections with indices

@author: Christopher J. Burke
"""

import numpy as np
import glob
import os
import shutil

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list

def intersect(*arrays):
    """ This only works if arrays are sorted and unique"""
    
    matched = np.array(list(set(arrays[0]).intersection(*arrays[1:])))
    return np.array([np.where(np.in1d(array, matched))[0] for array in arrays])

def copy_dir_diff(dir1, dir2, dirout):
    """ Copy files in dir1 that are missingin dir2 into dirout """
    print(dir1)
    fileList = glob.glob(os.path.join(dir1, '*'))
    for curFullFile in fileList:
        curFile = os.path.basename(curFullFile)
        checkFullFile = os.path.join(dir2, curFile)
        if os.path.isfile(checkFullFile):
            print('{0} exist'.format(curFile))
        else:
            print('{0} miss'.format(curFile))
            newFullFile = os.path.join(dirout, curFile)
            shutil.copyfile(curFullFile, newFullFile)
            
