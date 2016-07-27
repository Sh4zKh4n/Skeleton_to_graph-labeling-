# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:01:52 2016
@author: yuliya

Allows to calculate the skeleton graph from the input .h5
array.

Creates first a skeleton and then a NetxorkX graph using label_method module.
"""

import label_method as lm
import networkx as nx
import glob
import numpy as np
import tables as tb
import skeletonize as sk
from scipy.ndimage import distance_transform_edt
from scipy import ndimage
from graph_vis import calc
   
    

data = glob.glob('path/*.h5') #read .h5 files
data.sort()

#Skeleton graph transfrom
for i in data:
    dat = tb.openFile(i, mode='r')
    data_seg = np.copy(dat.root.segment)
    mask = data_seg > -1   #remove background
    seg = data_seg * mask
    seg_sort = sk.remove_isolated_br(seg, 1) #we leave only largest domain
    seg_f = (sk.fill_holes(seg_sort)).astype('int8') #fill the holes in order to avoid bubles in the skeleton

    dat.close()
    
    skel = sk.skeleton(seg_f)   #skeleton
    try:
        distance = distance_transform_edt(seg_f)  
        dist = skel*distance    #local thicknesses mapped on the skeleton
    except MemoryError:    #for large data volumes I had MemoryError sometimes at this step  
        dist = skel
           
    #Here we create .h5 files containing the skeletons, distance transform, list of thicknesses,
    #array of labeled skeleton branches. Any other arrays to conserve can be added by analogy
    fileName = 'path/skeletons_largedom/' + i[19:-3] + '.h5' #indices in i differs depending on the length 
                                                            #of path; name can be arbitrary
    shape = seg.shape
    atom_s = tb.UInt8Atom()
    atom_d = tb.UInt16Atom()
    filters = tb.Filters(complevel=5, complib='zlib')

    d = tb.open_file(fileName, 'w')
    ca_s = d.create_carray(d.root, 'skel', atom_s, shape,
                       filters=filters)
    ca_dist = d.create_carray(d.root, 'dist', atom_d, shape,
                       filters=filters)    
    
    ca_s[:,:,:] = skel[:,:,:]
    ca_dist[:,:,:] = dist[:,:,:]     
    
    n = lm.numb(skel)
    br = lm.label_br(n)
    br1 = lm.rem_bound(br)  #exclude boundary effects by removing boundary branches
    
    #create a list of mean thicknesses of the branches
    find_br = ndimage.find_objects(br1)
    un = np.unique(br1)[1:]
    t = []
    for k in un:
        mask = br1[find_br[k-1]]==k
        v = np.unique(mask * dist[find_br[k-1]])
        n1 = float(np.count_nonzero(v.ravel()))
        t = np.append(t, (np.sum(v.ravel())/n1))
            
    atom = tb.UInt16Atom()
    filters = tb.Filters(complevel=5, complib='zlib')
    shape = t.shape
    ca_t = d.create_carray(d.root, 'mean_thick', atom, shape,
                       filters=filters)
    ca_t[:] = t[:]
    
    shape = skel.shape
    atom = tb.UInt32Atom()
    filters = tb.Filters(complevel=5, complib='zlib')

    ca_br = d.create_carray(d.root, 'branches1', atom, shape,
                       filters=filters)       
    ca_br[:,:,:] = br1[:,:,:]
    d.close()
    
    #Create the graph from the skeleton
    nodes = lm.label_nodes(n)
    nodes1 = lm.rem_bound(nodes)
    G = lm.neigh_br(nodes1, br1)
    G = lm.create_con(G, br1, t)
   
#    nx.write_gpickle(G, 'path' + i[19:-3] + '_1.gpickle')   #If it is necessary to write the graph.

#This function deals with the break-ups.
calc()
