# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:51:33 2016

@author: yuliya

The module allows to visualize the graph. It also contains
all necessary function for the break-ups detection.
"""

import networkx as nx
import numpy as np
from mayavi import mlab
from nx_merge_nodes import merge_nodes
import itertools
import glob
from scipy import ndimage, spatial
import tables as tb
from skimage import measure

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
    mlab.figure(f, bgcolor=(0,0,0), size=(1200,1200))
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
     
        
        
def merge(graph, min_d):
    """
    This function merges close nodes into a single node conserving the connections,
    but in fact it is not used. 
    
    Parameters
    ------------
    graph : input graph with the node attributes 'posistion' storing the coordinates
    min_d : radius of merging. All the nodes that are closer than this radius 
            are merged into a single node.
            
    Return
    --------
    gr : computed graph 
    """
    p = nx.get_node_attributes(graph, 'position')
    n = graph.number_of_nodes()
    gr = graph.copy()
    edges = list()
    ed = graph.edges()
    for comb in itertools.combinations(p, 2):  
        i = comb[0]
        j = comb[1]
        nodes = []
        dist = np.linalg.norm(p[i] - p[j])
        if (dist < min_d) and ((i,j) in ed):
            if (i not in edges) and (j not in edges):
                nodes = np.append(nodes, i)
                nodes = np.append(nodes, j)
                pos = (p[i] + p[j])/2
                gr.remove_edge(i, j)
                merge_nodes(gr, nodes, n+1, position = pos)
                edges.append(i)
                edges.append(j)
                n+=1
   
    return gr
 
def remove_short_branches(graph, min_l):
    """
    The function removes all branches in the graph that are shorter than a 
    given length. Is not used.
    
    Parameters
    ----------
    graph : input graph with the node attributes 'posistion' storing the coordinates
    min_l : minimum length of the branch
    
    Return
    --------
    gr : computed graph 
    """
    gr = graph.copy()
    p = nx.get_node_attributes(gr, 'position')
    bound = []
    for i in gr.nodes():
        if nx.degree(gr, i)==1:
            bound = np.append(bound, i)
    for i in bound:
        j = list(nx.neighbors(gr, i))[0]
        d = np.linalg.norm(p[i] - p[j])
        if d < min_l :
            gr.remove_edge(min(int(i), j), max(int(i), j), 0)
            gr.remove_node(int(i))
    return gr

def remove_terminal_br(graph):
    """
    Removes all terminal branches in a given graph ONCE.
    """
    g = graph.copy()
    for i in graph.nodes():
        if nx.degree(g, i)==1 :
            g.remove_node(i)
    return g 
  
  
def remove_nodes(g):
    """
    Modifies the INPUT graph by removing all nodes, whose degree is 2.
    Is not used.
    """
    g1 = g.copy()
    nodes = []
    le = nx.get_edge_attributes(g, 'length')
    for i in g1.nodes():
        if nx.degree(g1, i) == 2:
            nodes = np.append(nodes, i)
    for i in nodes:
        i = int(i)
        n = list(nx.neighbors(g, i))
        le = nx.get_edge_attributes(g, 'length')
        th = nx.get_edge_attributes(g, 'thick')
        num = nx.get_edge_attributes(g, 'number')
        if len(n)==1:
            continue
        else:
            n1 = n[0]
            n2 = n[1]
            k1 = num[min(n1, i), max(i, n1), 0]
            k2 = num[min(n2, i), max(i, n2), 0]
            l = le[min(n1, i), max(i, n1), 0] + le[min(i, n2), max(i, n2), 0]
            t = (th[min(n1, i), max(i, n1), 0] + th[min(i, n2), max(i, n2), 0])/2.
            le.update({(min(n1, n2), max(n2, n1), 0): l})
            th.update({(min(n1, n2), max(n2, n1), 0): t})
            num.update({(min(n1, n2), max(n2, n1), 0): k1+k2})
            g.remove_node(i)
            g.add_edge(n1, n2, length = l, thick = t, number = k1 + k2)
    return g
    

def dist(gr):  
    """
    Computes the euclidean distances between all the nodes.
    Is not used.
    """
    p = nx.get_node_attributes(gr, 'position')
    d = []
    for comb in itertools.combinations(p, 2):  
        i = comb[0]
        j = comb[1]
        dist = np.linalg.norm(p[i] - p[j])
        d = np.append(d, dist)
    return d

   



def branches_map(g1, skel2, bran1, thick):
    """
    Parameters
    --------
    g1 : graph of the first image
    skel2: full skeleton of the second image
    bran1: array of labeled branches of the first image
    thick: distance transfrom that stores local thicknesses of the skeleton,
           must have the same shape as skel2
    
    Return
    ------
    branches: dictionnary {k: sl} where k is the number of broken branch and
              sl is the slice of the skel2 array with this branch.
    pos: a list of coordinates of ruptures
    """
    
    number = nx.get_edge_attributes(g1, 'number') #dictionnary edges-numbers

    find_br = ndimage.find_objects(bran1)
    skel_pts = np.transpose(np.nonzero(skel2))
    skel_tree = spatial.KDTree(skel_pts) #We use Kd tree to work with the skeleton,
                                        #it stores the posistions of all nonzero points
    branches = {}
    pos = list()
    for j in number:
        br = number[j]
        flag = False  #We use flag here to break the cicle when the rupture event is found
        for k in br:
            sl = find_br[k-1]
            mask = (bran1[sl]==k).astype('uint8')
            start = []           
            for i in range(3):
                start = np.append(start, sl[i].start)
                points1 = np.transpose(np.nonzero(mask))
            for i in range(points1.shape[0]):
                point = points1[i] + start
                rad = thick[int(point[0]), int(point[1]), int(point[2])] #choose the radius of the ball 
                                                                            #as the local thickness of branch
                if (rad<5) and (points1.shape[0]>8): #this criteria is experimental, applied to very thin 
                    rad = 5                                #branches
                ball = skel_tree.query_ball_point(point, rad) #ball contains the points of the next skeleton
                if len(ball)==0:                #If it is empty than the branch has broken
                    branches.update({k: sl})
                    pos.append(point)
                    flag = True
                    break
            if flag:
                break
    return branches, pos



def vis_breakups():
    """
    This function allows to visualize break-ups as the highlighted branches in
    the original skeleton.
    """
    data_sk = glob.glob('/backup/yuliya/vsi05/skeletons_largdom/*.h5')
    data_sk.sort()
    
    data_br = glob.glob('/backup/yuliya/vsi05/breakups_con_correction/dict/*.npy')
    data_br.sort()
    
    for i,j in zip(data_sk[27:67][::-1], data_br[19:][::-1]):
        d = tb.openFile(i, mode='r')
        bran1 = np.copy(d.root.branches)
        skel1 = np.copy(d.root.skel)
        d.close() 

        br = np.load(j).item()
        
        mlab.figure(1, bgcolor=(1,1,1), size=(1200,1200))
        surf_skel = mlab.contour3d(skel1, colormap='hot')
        surf_skel.actor.property.representation = 'points'
        surf_skel.actor.property.line_width = 2.3
        mask = np.zeros_like(skel1)
        for k in br:
            sl = br[k]            
            mask[sl] = (bran1[sl] == k)
        surf_mask = mlab.contour3d(mask.astype('uint8'))
        surf_mask.actor.property.representation = 'wireframe'
        surf_mask.actor.property.line_width = 5
        
        mlab.savefig('/backup/yuliya/vsi05/breakups_con_correction/video_hot/' + i[39:-3] + '.png')
        mlab.close()
        

            
def calc():
    """
    The function that implements break-ups detection algorithm and
    write the data in a file. Input skeletons and graphs must correspond each
    other.
    """
    
    #input graphs
    data_gr = glob.glob('/path/graphs_largedom/*_1.gpickle')
    data_gr.sort()

    #input skeletons
    data_sk = glob.glob('/path/skeletons_largedom/*.h5')
    data_sk.sort()

    for i, j, w in zip(data_gr[:-1][::-1], data_sk[:-1][::-1], data_sk[1:][::-1]):
        g = nx.read_gpickle(i)

        #We simplify the graph
        g1 = remove_terminal_br(g)
        g1 = remove_nodes(g1)
        
        #We use labeled branches arrays, skeletons, distans transfrom. All are of the same shape.
        d = tb.openFile(j, mode='r')
        bran1 = np.copy(d.root.branches1)
        skel1 = np.copy(d.root.skel) #this 'current' skeletons are used only for visualization.
        thick = np.copy(d.root.dist)
        d.close()
        
        #Skeletons shifted in time by one step.
        d = tb.openFile(w, mode='r')
        skel2 = np.copy(d.root.skel)
        d.close()
        
        b, pos = branches_map(g1, skel2, bran1, thick)
        
        #Visualization of break-ups if necessary.
        mlab.figure(1, bgcolor=(1,1,1), size=(1500,1500))
        mlab.contour3d(skel1, colormap='hot')
        mask = np.zeros_like(skel1)
        for k in b:
            sl = b[k]            
            mask[sl] = (bran1[sl] == k)
        mlab.contour3d(mask.astype('uint8'))
        
        np.save('/path/breakups/' + j[41:-3] + '_1.npy', b)  #j[41:-3] must be replaced according to the path
        np.save('/path/breakups/' + j[41:-3] + '_pos.npy', pos)
        mlab.savefig('/path/breakups/' + j[41:-3] + '_1.png')
        mlab.close()
        
     
