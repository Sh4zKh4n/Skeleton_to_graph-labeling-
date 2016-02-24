# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:46:13 2016

Allows to calculate the graph ou of the skeleton
directly using the function from the label_method module.
"""
import numpy as np
import label_method as lb

data = np.load('/home/yuliya/Desktop/v34_0140_skel(200,600,400,700,500,700)filled.npy', mmap_mode='r')

dat1 = data[:200, :150, :100]

def calc_graph(dat):
    """
    Creates a graph from the input skeleton data.
    
    The skeleton is segmented to the regions depending on the
    neighborhood of each points: the nodes have more than 3 neighbors
    and the branches have 1 or 2 neighbors. The isolated points are
    not treated. The nodes are converted to the nodes of a networkx 
    graph with the attribute 'neigh', which is the array of the adjoint
    branches. Binary dilation is used then to establish the connections
    between the nodes and the branches and add the edges to the graph, 
    creating necessary end nodes. All branches have the attribute 'length'
    that is equal the the number of elements (points) in the branch.
    
    Parameters
    ---------
    dat : 3D binary array
    
    Returns
    -------
    G : a graph
    """
    dat_n = lb.numb(dat)
    dat_nodes = lb.label_nodes(dat_n)
    dat_br = lb.label_br(dat_n)
    
    dat_G = lb.neigh_br(dat_nodes, dat_br)
    G = lb.create_con(dat_G, dat_br)
    return G

def isolated_points(graph, dat):
    """
    Treats the isolated points and adds them to a graph as
    the nodes.
    
    Parameters
    ---------
    graph : a graph
    dat : 3D skeleton binary array
    
    Returns
    -------
    G : a copy of input graph with isolated nodes
    """
    G = graph.copy()
    dat_n = lb.numb(dat)
    dat_is = lb.label_iso(dat, dat_n)
    G = lb.add_isol_points(G, dat_is)
    return G
    

