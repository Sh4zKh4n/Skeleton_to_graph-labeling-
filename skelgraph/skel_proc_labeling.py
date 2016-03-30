# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:46:13 2016

Allows to calculate the graph ou of the skeleton
directly using the function from the label_method module.
"""
import numpy as np
import label_method as lb


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
    # this function could go to the label_method module instead
    dat_n = lb.numb(dat)
    dat_nodes = lb.label_nodes(dat_n)
    dat_br = lb.label_br(dat_n)
    dat_nodes = lb.rem_bound(dat_nodes)
    dat_br = lb.rem_bound(dat_br)
    
    dat_G = lb.neigh_br(dat_nodes, dat_br)
    G = lb.create_con(dat_G, dat_br)
    return G

    

