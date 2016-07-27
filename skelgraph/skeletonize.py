# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:44:56 2016

@author: yuliya

This module contains several functions for the initial data preparation and their skeletonization.

"""

#import skel  #I used this module initially, but now one can use skimage function
import numpy as np
import tables as tb
from skimage import measure, morphology
from scipy import ndimage


def op(data1=str):
    """
    Opens .h5 file.
    """
    return tb.openFile(data1, mode='r') 

def copy_seg(data1, x1, x2, y1, y2, w1, w2):
    """
    Creates an npy array as a segment of the input file.
    """
    return np.copy(data1.root.segment[x1:x2, y1:y2, w1:w2])

def fill_holes(data1):
    """
    Returns binary image without holes.
    """
    return ndimage.binary_fill_holes(data1)    



def remove_isolated_br(data1, n):
    """
    Leaves n largest connected regions of the input image.
    
    Parameters
    ---------
    data1 : input array
    n : int, the number of regions to return
    
    Returns
    -------
    mask : array
    """
    lab = measure.label(data1, connectivity=data1.ndim, background=0)
    props = measure.regionprops(lab)
    a=[]
    for i in props:
        a = np.append(a, i.area)
    a_s = np.argsort(a)[::-1]
    mask = np.zeros(data1.shape)
    for i in range(n):
        el = a_s[i] + 1
        m = lab==el
        mask = mask + m
    return mask.astype('uint8')

def skeleton(data1):
    """
    Returns the skeleton of the input binary image using
    the function skeletonize_3d (#initially skel.compute_thin_image).
    """
    return morphology.skeletonize_3d(data1)
     
