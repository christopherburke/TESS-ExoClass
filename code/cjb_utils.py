# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 14:46:59 2017
Some personal numpy array filtering and finding intersections with indices

@author: Christopher J. Burke
"""

import numpy as np

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list

def intersect(*arrays):
    """ This only works if arrays are sorted and unique"""
    
    matched = np.array(list(set(arrays[0]).intersection(*arrays[1:])))
    return np.array([np.where(np.in1d(array, matched))[0] for array in arrays])
