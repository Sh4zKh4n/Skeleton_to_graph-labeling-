# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:05:58 2016

Contains the functions necessary to transform the 3D skeleton
into a graph using the initial labeling of the input binary array as
nodes and branches on the base of the neighborhood of each point.
"""

from skimage import measure, morphology
import numpy as np
from skimage import util
from scipy import ndimage
import networkx as nx


def numb(c) :
    """
    Counts the number of neighboring 1 for every nonzero
    element in 3D.
    
    Parameters
    ------
    skel_mat : 3d binary array
    
    Returns
    ------
    arr : The array of uint8 of the same size as skel_mat 
    with the elements equal to the number of neighbors 
    for the current point.
    
    Examples
    ------
    >>> a = np.random.random_integers(0,1,(3,3,3))
    >>> a 
    array([[[0, 0, 1],
            [0, 1, 1],
            [1, 0, 1]],        
           [[1, 0, 1],
            [1, 0, 1],
            [1, 1, 0]],        
           [[1, 1, 1],
            [1, 0, 0],
            [1, 1, 0]]])
    >>> neigh = numb(a)
    >>> neigh 
    array([[[ 0,  0,  4],
            [ 0, 10,  6],
            [ 4,  0,  4]],            
           [[ 5,  0,  6],
            [10,  0,  9],
            [ 7, 10,  0]],            
           [[ 4,  7,  3],
            [ 8,  0,  0],
            [ 5,  6,  0]]], dtype=uint8)
    """
    c_pad = util.pad(c, 1, 'constant')
    mask = c_pad > 0
    fil = 3**3 * ndimage.uniform_filter(c_pad.astype('float'), 
                                        size=3) - 1
    return (fil * mask)[1:-1, 1:-1, 1:-1].astype('uint8')


def label_br(num):
    """
    Returns the array where each region of points with the number
    of neighbors 1 or 2 is labeled from 1 to N, where N is the
    number of regions.
    
    Parameters
    ---------
    num : a 3D array counting the number of neighbors at each point
    
    Returns
    -------
    arr : an array of the same shape as the input with labeled regions
    
    Examples
    --------
    >>>c = np.zeros((3,4,3))
    >>>c[1,:,0] = c[1, 1,:] = 1
    array([[[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]],
           [[ 1.,  0.,  0.],
            [ 1.,  1.,  1.],
            [ 1.,  0.,  0.],
            [ 1.,  0.,  0.]],
           [[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]]])
    >>>num = numb(c)
    >>>num
    array([[[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]],
           [[2, 0, 0],
            [3, 4, 1],
            [3, 0, 0],
            [1, 0, 0]],
           [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]], dtype=uint8)
    >>>label_br(num)
    array([[[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]],
           [[1, 0, 0],
            [0, 0, 2],
            [0, 0, 0],
            [3, 0, 0]],
           [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]], dtype=uint8)   
    """
    br2 = num == 2
    br1 = num == 1
    branches = br1 + br2
    label_branches = measure.label(branches, background=0)
    return (label_branches + np.ones_like(label_branches)).astype('uint8')
    

 
def label_nodes(num):
    """
    Returns the array where each region of points with the number
    of neighbors more than 3 is labeled from 1 to N, where N is the
    number of regions.
    
    Parameters
    ---------
    num : a 3D array counting the number of neighbors at each point
    
    Returns
    -------
    arr : an array of the same shape as the input with labeled regions
    
    Examples
    --------
    >>>c = np.zeros((3,4,3))
    >>>c[1,:,0] = c[1, 1,:] = 1
    array([[[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]],
           [[ 1.,  0.,  0.],
            [ 1.,  1.,  1.],
            [ 1.,  0.,  0.],
            [ 1.,  0.,  0.]],
           [[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]]])
    >>>num = numb(c)
    >>>num
    array([[[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]],
           [[2, 0, 0],
            [3, 4, 1],
            [3, 0, 0],
            [1, 0, 0]],
           [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]], dtype=uint8)
    >>>label_nodes(num)
    array([[[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]],
           [[0, 0, 0],
            [1, 1, 0],
            [1, 0, 0],
            [0, 0, 0]],
           [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]], dtype=uint8) 
    """
    no = num>2
    label_n = measure.label(no, background=0)
    # uint8 might be a problem for large arrays with more than 255 nodes...
    return (label_n+np.ones((label_n.shape))).astype('uint8')
    
    
def label_iso(skel, num):
    """
    Returns the array where each region labels isolated points
    (the number of neighbors 0) from 1 to N, where N is the
    number of regions.
    
    Parameters
    ---------
    skel : binary array of skeleton
    num : a 3D array counting the number of neighbors at each point
    
    Returns
    -------
    arr : an array of the same shape as the input with labeled regions
    """
    iso_z = num==0
    iso = skel*iso_z
    label_iso = measure.label(iso, background=0)
    return (label_iso+np.ones((label_iso.shape))).astype('uint8')    
    


def neigh_br(arr_nodes, arr_br):
    """
    Creates a graph with nodes from the array of nodes and branches
    by considering the branches that adjoin to the current node.
    The labels of these branches are set as an attribute of the
    node.
    
    Parameters
    ---------
    arr_nodes : array of labeled regions that are considered as nodes
    arr_br : array of labeled regions considered as branches
    
    Returns
    ------
    G : a graph of nodes with the attribute 'neigh' 
    
    Examples
    -------
    >>>arr_nodes = label_nodes(num)
    >>>arr_nodes
    array([[[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]],
           [[0, 0, 0],
            [1, 1, 0],
            [1, 0, 0],
            [0, 0, 0]],
           [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]], dtype=uint8) 
    >>>arr_br = label_br(num)
    >>>arr_br
    array([[[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]],
           [[1, 0, 0],
            [0, 0, 2],
            [0, 0, 0],
            [3, 0, 0]],
           [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]], dtype=uint8) 
    >>>G = neigh_br(arr_nodes, arr_br)
    >>>G.nodes()
    [1]
    >>>nx.get_node_attributes(G, 'neigh')
    {1: array([1, 2, 3])}
    
    """
    G = nx.Graph()
    el = morphology.cube(3)
    nodes = np.unique(arr_nodes)[1:]
    for j in nodes:
        mask = arr_nodes==j
        dil = morphology.binary_dilation(mask, el)
        vol_c = arr_br * (dil.astype('uint8'))
        G.add_node(j, neigh=np.unique(vol_c)[1:])
    return G


def create_con(graph, arr_br):
    """
    Creates the edges in the input graph with nodes that have
    the labels of adjoint branches as an attribute 'neigh'.
    End nodes are created while adding the edges to the graph and
    do not have the attribute 'neigh'. The edges have the attribute
    'length' that counts the number of elements in the branch and
    corresponds to the length of the branch.
    
    Parameters
    ---------
    graph : a graph of nodes
    arr_br : array of labeled regions considered as branches
    
    Returns
    -------
    G : a copy of the input graph with edges and end nodes added
    
    Examples
    --------
    >>>c = np.zeros((3,4,3))
    >>>c[1,:,0] = c[1, 1,:] = 1
    array([[[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]],
           [[ 1.,  0.,  0.],
            [ 1.,  1.,  1.],
            [ 1.,  0.,  0.],
            [ 1.,  0.,  0.]],
           [[ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.],
            [ 0.,  0.,  0.]]])
    >>>G = neigh_br(arr_nodes, arr_br)
    >>>G.nodes()
    [1]
    >>>nx.get_node_attributes(G, 'neigh')
    {1: array([1, 2, 3])}
    >>>G_con = create_con(G, arr_br)
    >>>G_con.nodes()
    [1, 2, 3, 4]
    >>>G_con.edges()
    [(1, 2), (1, 3), (1, 4)]
    >>>nx.get_edge_attributes(G_con, 'length')
    {(1, 2): 1, (1, 3): 1, (1, 4): 1}    
    """
    G = graph.copy()
    n = nx.get_node_attributes(G, 'neigh')
    nodes = G.nodes()
    list_br = [0]
    for i in nodes:
        for k in n[i]:
            flag = True            
            mask = arr_br==k
            # you could compute all the length outside of the i loop
            # and store them into a dictionary
            # here they are computed many times
            le = len(np.transpose(np.nonzero(mask)))
            for l in nodes:
                if ((k in n[l]) and (i != l)):
                    G.add_edge(i,l, length = le)
                    flag = False
                    break               
            # not sure I understand the use of flag
            if flag:
                G.add_edge(i, G.number_of_nodes()+1, length = le)
            list_br = np.append(list_br, k)
    for i in np.unique(arr_br):
        if (i not in list_br):
            mask = arr_br==i
            le = len(np.transpose(np.nonzero(mask)))
            G.add_edge(G.number_of_nodes()+1, G.number_of_nodes()+2, length = le) 
    return G


def add_isol_points(graph, arr_is):
    """
    Adds the points as nodes to the input graph using the
    array of labeled isolated regions.
    
    Parameters
    ---------
    graph : a graph
    arr_is : an array of labeled regions considered as isolated points
    
    Returns
    -------
    G : a copy of input graph with isolated nodes added
    """
    G = graph.copy()
    for i in np.unique(arr_is)[1:]:
        G.add_node(G.number_of_nodes()+1)
    return G
    
