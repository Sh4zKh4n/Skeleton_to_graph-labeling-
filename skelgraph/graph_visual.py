# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:36:33 2016

@author: yuliya
"""


import networkx as nx
import numpy as np
from mayavi import mlab

from numpy import array


def vis(g, f):
    """
    Visualization of the graph as the set of glyphes and tubes.
    
    Parameters
    ------------
    g : input graph with the node attributes 'posistion' storing the coordinates
        of nodes in space
    f : number of the current figure where to show the graph
    
    """
    mlab.figure(f, bgcolor=(1,1,1), size=(1200,1200))
    pos = nx.get_node_attributes(g, 'position')
    xp, yp, zp = [], [], []
    s = []
    for i in pos:
        xp = np.append(xp, pos[i][0])
        yp = np.append(yp, pos[i][1])
        zp = np.append(zp, pos[i][2])
        s = np.append(s, nx.degree(g, i)/10.)   
    
    glyph1 = mlab.points3d(xp, yp, zp, s)
    glyph1.glyph.glyph.range = array([ 0. ,  0.4])
    glyph1.glyph.glyph.range = array([ 0. ,  0.4])
    glyph1.glyph.glyph.range = array([ 0. ,  0.4])
    glyph1.glyph.glyph.scale_factor = 10.0
    th = nx.get_edge_attributes(g, 'thick')
    v1, v2 = [], []
    for i in g.edges():
        v1 = np.append(v1, i[0])
        v2 = np.append(v2, i[1])
    for i, j, k in zip(v1, v2, th):
        xl = [pos[i][0], pos[j][0]]
        yl = [pos[i][1], pos[j][1]]
        zl = [pos[i][2], pos[j][2]]
        lines = mlab.plot3d(xl, yl, zl, tube_radius=th[k]/100.)
     