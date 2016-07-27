# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 09:28:13 2016

@author: yuliya

This module contains all the functions for the skeleton graph
analysis.
"""

import numpy as np
import glob
import tables as tb
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx
#from graph_vis import remove_terminal_br, remove_nodes
import itertools
import random
from scipy.optimize import curve_fit
from mayavi import mlab
from scipy import ndimage

def len_hist():
    """
    length histogram
    """
    data = glob.glob('/backup/yuliya/v34/graphs/*_1.gpickle')
    data.sort()
    
    for i in data:
        G = nx.read_gpickle(i)
        leng = nx.get_edge_attributes(G, 'length')
        hist = list()
        for (u,v,key) in leng:
            hist.append(leng[u,v,key])
        np.save('/backup/yuliya/v34/length_branches/' + i[26:-10] + '.npy', hist)

#--------------------------
#input data
d1 = list()
d = list()

data1 = glob.glob('/backup/yuliya/v34/breakups_con/fixed/*.npy')
data1.sort()
data = glob.glob('/backup/yuliya/v34/graphs/fixed/*_no_term_br.gpickle')
data.sort()

d1.append(data1)
d.append(data)

data1 = glob.glob('/backup/yuliya/vsi05/breakups_con_correction/dict/fixed/*.npy')
data1.sort()
data = glob.glob('/backup/yuliya/vsi05/graphs_largdom/fixed/*.gpickle')
data.sort()

d1.append(data1)
d.append(data)

data1 = glob.glob('/backup/yuliya/v30/breakups_con/fixed/*.npy')
data1.sort()
data = glob.glob('/backup/yuliya/v30/graphs_largedomain/fixed/*.gpickle')
data.sort()

d1.append(data1)
d.append(data)

data1 = glob.glob('/backup/yuliya/v35/breakups_con/fixed/*.npy')
data1.sort()
data = glob.glob('/backup/yuliya/v35/graphs_largedomain/fixed/*_no_term_br.gpickle')
data.sort()

d1.append(data1)
d.append(data)

colors = cm.autumn(np.linspace(0, 1, 114))
colors1 = cm.winter(np.linspace(0, 1, 12))
t = np.load('/backup/yuliya/vsi01/vsi01_TimeInSec.npy')
u = 0
ell = np.load('/backup/yuliya/vsi05/dataellg_ellp.npy')[8:-10]
ellv34 = np.load('/backup/yuliya/v34/dataellg_ellp.npy')[13:-1][::1]
ellv34[7:-1] = ellv34[8:]
ellv34[10:-1] = ellv34[11:]
ellv34[45:-1] = ellv34[46:]
ellvsi05 = np.load('/backup/yuliya/vsi05/dataellg_ellp.npy')[8:-10]
ellvsi05[19:-1] = ellvsi05[20:]
ellv30 = np.load('/backup/yuliya/v30/dataellg_ellp.npy')[10:-1]
ellv30[61:-1] = ellv30[62:]
ellv30[63:-1] = ellv30[64:]
ellv35 = np.load('/backup/yuliya/v35/dataellg_ellp.npy')[22:-1]
ellv35[4:-1] = ellv35[5:]
ellv35[6:-4] = ellv35[10:]
ellv35[9:-2] = ellv35[11:]
ellv35[13:-1] = ellv35[14:]
ellv35[22:-1] = ellv35[23:]

ellv23 = np.load('/backup/yuliya/v23/dataellg_ellp.npy')[53:]

ell = list()
ell.append(ellv34)
ell.append(ellvsi05)
ell.append(ellv30)
ell.append(ellv35)
#eu = np.load('/backup/yuliya/v34/graphs/Euler.npy')[12:-1]
#vol = np.load('/backup/yuliya/v34/volumes_pd.npy')[13:-1]
#vol = np.load('/backup/yuliya/v34/volumes_pd.npy')[60:114]
eu = list()
c_t = list()
vert = list()
bran = list()
n_of_br = list()
n_tot = []

pa = ('v34', 'vsi05', 'v30', 'v35') #experiments
ti = [24., 32., 8., 5.2] #time intervals between the images

#-----------------------


def fit_d(x, a, b):
    return a * np.exp( - x * b)
    
def length_distr_total():
    """    
    Length distribution total
    """
    for data, data1, path, el, t in zip(d, d1, pa, ell, ti):
        l1 = list() #all branches
        l2 = list() #breaking branches
        u = 0
        for i,j,le in zip(data, data1, el):
            g = nx.read_gpickle(i)
            br = np.load(j).item()
            length = nx.get_edge_attributes(g, 'length')
            l_mean = np.mean(np.asarray(length.values()))
            number = nx.get_edge_attributes(g, 'number')
            num_dict = {}
            for k in number:
                for v in number[k]:
                    num_dict.setdefault(v, []).append(k)
                
            for k in length.values():
                l1.append(k/float(1))
            for k in br.keys():
                for l in num_dict[k]:
                    l2.append(length[l]/float(1))
            u+=1
        hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.06))
        hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.06))
        center1 = (bins1[:-1] + bins1[1:])/2
        center2 = (bins2[:-1] + bins2[1:])/2
        #save to file if necessary
        #np.save('/home/yuliya/codes/lengths/' + path + '/total_lno_all_len.npy', center1)    
        #np.save('/home/yuliya/codes/lengths/' + path + '/total_lno_break_len.npy', center2)
        #np.save('/home/yuliya/codes/lengths/' + path + '/total_lno_all.npy', hist1/float(len(l1)))
        #np.save('/home/yuliya/codes/lengths/' + path + '/total_lno_break.npy', hist2/float(len(l1)))
    #
    ##plt.figure(2)
    #hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.06))
    #hist2, bins2 = np.histogram(l2, np.arange(0, max(l2)+1, 0.06))
    #center1 = (bins1[:-1] + bins1[1:])/2
    #center2 = (bins2[:-1] + bins2[1:])/2
    #plt.plot(center1, hist1/float(len(l1)), '.', color='red', label = 'all branches')
    #plt.plot(center2, hist2/float(len(l1)), '.', color='blue', label = 'breaking branches')
    #



    """
    Probability as a function of length.
    """

    hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.1))
    hist2, bins2 = np.histogram(l2, np.arange(0, max(l2)+1, 0.1))
    center1 = (bins1[:-1] + bins1[1:])/2
    plt.plot(center1, hist2/hist1.astype('float32'), '.', color='green', label = 'v35')


def short_branches():
    """
    Visualization of short branches of the skeleton.
    
    """
    data1_sk = glob.glob('/backup/yuliya/vsi05/skeletons_largdom/*.h5')
    data1_sk.sort()

    for i,j, k in zip(d[1][37:47], data1_sk[46:56], ell[1][37:47]):
        g = nx.read_gpickle(i)
        dat = tb.openFile(j)
        skel = np.copy(dat.root.skel)
        bra = np.copy(dat.root.branches)
        mask = np.zeros_like(skel)    
        dat.close()
    
        length = nx.get_edge_attributes(g, 'length')
        number = nx.get_edge_attributes(g, 'number')
        num_dict = {}
        for m in number:
            for v in number[m]:
                num_dict.setdefault(v, []).append(m)
        find_br = ndimage.find_objects(bra)
        for l in list(length.keys()):
            if length[l]<0.5*k: #Criteria
                for b in number[l]:
                    mask[find_br[b-1]] = bra[find_br[b-1]]==b
        mlab.figure(bgcolor=(1,1,1), size=(1200,1200))
        mlab.contour3d(skel, colormap='hot')
        mlab.contour3d(mask)
        mlab.savefig('/backup/yuliya/vsi05/skeletons/short_bran/'+ i[42:-10] + '.png')
        mlab.close()
    
            


def length_distr_dyn():
    """
    Length distributions evolving in time.
    """
    av = [20, 12, 12, 13] #Number of images for avereging for each experiment. Here 
                            #we use 20 images for one plot for v34 experiment etc.

    for data, data1, path, el, a in zip(d, d1, pa, ell, av):
        l1 = list() #all branches
        l2 = list() #breaking branches
        u = 0
        count=0
        co=0
        for i,j,le in zip(data, data1, el):
                g = nx.read_gpickle(i)
                br = np.load(j).item()
                length = nx.get_edge_attributes(g, 'length')
                l_mean = np.mean(np.asarray(length.values()))
                number = nx.get_edge_attributes(g, 'number')
                num_dict = {}
                for k in number:
                    for v in number[k]:
                        num_dict.setdefault(v, []).append(k)
                
                for k in length.values():
                    l1.append(k/float(l_mean))
                for k in br.keys():
                    for l in num_dict[k]:
                        l2.append(length[l]/float(l_mean))
                u+=1
                if count>a:
                    count=-1
                    hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.06))
                    hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.06))
                    center1 = (bins1[:-1] + bins1[1:])/2
                    center2 = (bins2[:-1] + bins2[1:])/2
                    
                    #write to file if necessary
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_lm_all_len' + str(co)+'.npy', center1)    
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_lm_break_len' + str(co)+'.npy', center2)
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_lm_all' + str(co)+'.npy', hist1/float(len(l1)))
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_lm_break' + str(co)+'.npy', hist2/float(len(l1)))
                    l1 = list()
                    l2 = list()
                    co+=1
                count+=1

def thick_distr_dyn():
    """
    Thickness distributions evolving in time.
    """
    av = [20, 12, 12, 13]

    for data, data1, path, el, a in zip(d, d1, pa, ell, av):
        l1 = list() #all branches
        l2 = list() #breaking branches
        
        count=0
        co=0
        for i,j,le in zip(data, data1, el):
                g = nx.read_gpickle(i)
                br = np.load(j).item()
                thick = nx.get_edge_attributes(g, 'thick')
                th_mean = np.mean(np.asarray(thick.values()))
                number = nx.get_edge_attributes(g, 'number')
                num_dict = {}
                for k in number:
                    for v in number[k]:
                        num_dict.setdefault(v, []).append(k)
                for k in thick.values():
                    l1.append(k/float(th_mean))
                for k in br.keys():
                    for l in num_dict[k]:
                        l2.append(thick[l]/float(th_mean))
                if count>a:
                    count=-1
                    hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.06))
                    hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.06))
                    center1 = (bins1[:-1] + bins1[1:])/2
                    center2 = (bins2[:-1] + bins2[1:])/2
                    #write to file if necessary
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_rm_all_len' + str(co)+'.npy', center1)    
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_rm_break_len' + str(co)+'.npy', center2)
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_rm_all' + str(co)+'.npy', hist1/float(len(l1)))
                    np.save('/home/yuliya/codes/lengths/' + path + '/dyn_rm_break' + str(co)+'.npy', hist2/float(len(l1)))
                    l1 = list()
                    l2 = list()
                    co+=1
                count+=1

def degree_nodes():
    """
    Degree of nodes histogram.
    """
    for i in d[0][::5]:
        g = nx.read_gpickle(i)
        hist = nx.degree_histogram(g)
        n = nx.number_of_nodes(g)
        np.save('/home/yuliya/codes/lengths/v34/degree' + i[30:-10] + '.npy', np.asarray(hist)/float(n))


def thick_distr_total():
    """
    Thickness distributions.
    """
    for data, data1, path, el in zip(d, d1, pa, ell):
        l1 = list() #all branches
        l2 = list() #breaking branches
        u = 0
        for i,j,le in zip(data, data1, el):
            g = nx.read_gpickle(i)
            br = np.load(j).item()
            thick = nx.get_edge_attributes(g, 'thick')
            th_mean = np.mean(np.asarray(thick.values()))
            number = nx.get_edge_attributes(g, 'number')
            num_dict = {}
            for k in number:
                for v in number[k]:
                    num_dict.setdefault(v, []).append(k)
            for k in thick.values():
                l1.append(k/float(th_mean))
            for k in br.keys():
                for l in num_dict[k]:
                    l2.append(thick[l]/float(th_mean))
            u+=1
        hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.05))
        hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.05))
        center1 = (bins1[:-1] + bins1[1:])/2
        center2 = (bins2[:-1] + bins2[1:])/2
        #write to file if necessary
#        np.save('/home/yuliya/codes/lengths/' + path + '/total_rm_all_len.npy', center1)    
#        np.save('/home/yuliya/codes/lengths/' + path + '/total_rm_break_len.npy', center2)
#        np.save('/home/yuliya/codes/lengths/' + path + '/total_rm_all.npy', hist1/float(len(l1)))
#        np.save('/home/yuliya/codes/lengths/' + path + '/total_rm_break.npy', hist2/float(len(l1)))
    plt.plot(center1, hist1/float(len(l1)), '.', color='red', label = 'all branches')
    plt.plot(center2, hist2/float(len(l1)), '.', color='blue', label = 'breaking branches')
    plt.legend()
    plt.xlabel('d/l_typ', fontsize=18)
    plt.ylabel('P(d/l_typ)', fontsize = 18)


    """
    Probability as a function of thickness.
    """
    hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.02))
    hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.02))
    center1 = (bins1[:-1] + bins1[1:])/2
    plt.plot(center1, hist2/hist1.astype('float32'), '.', color='blue', label = 'v35')
    plt.legend()
    plt.xlabel('d/l_typ', fontsize=18)
    plt.ylabel('P(d/l_typ)', fontsize = 18)
    
    f = curve_fit(fit_d, center1[10:20], (hist2/hist1.astype('float32'))[10:20], p0=(1,0.2))[0]
    plt.plot(center1, fit_d(center1, f[0], f[1]), '.', color='blue', label = 'fitting vsi05')

def mean_length():
    """
    Time dependence of the mean length and mean thickness.
    """    

    l1 = list() #all branches
    l2 = list() #breaking branches
    d_all = glob.glob('/backup/yuliya/v35/graphs_largedomain/*_no_term_br.gpickle')
    d_all.sort()
    timev34 = np.load('/backup/yuliya/v35/v35_TimeInSec.npy')[12:]
    for i in d_all:
        g = nx.read_gpickle(i)
        length = nx.get_edge_attributes(g, 'length')
        thick = nx.get_edge_attributes(g, 'thick')
        m_le = np.mean(np.asarray(length.values()))
        m_th = np.mean(np.asarray(thick.values()))
        l1.append(m_le * 1.1)
        l2.append(m_th * 1.1)
    np.save('/home/yuliya/codes/lengths/v35/time_mean.npy', timev34)
    np.save('/home/yuliya/codes/lengths/v35/length_mean_t.npy', l1)
    np.save('/home/yuliya/codes/lengths/v35/thick_mean_t.npy', l2)


def aspect_rat_distr_total():
    """
    Aspect ratio distribution.
    """

    for data, data1, path, el, t in zip(d, d1, pa, ell, ti):
        l1 = list() #all branches
        l2 = list() #breaking branches
        u = 0
        for i,j in zip(data, data1):
            g = nx.read_gpickle(i)
            br = np.load(j).item()
            thick = nx.get_edge_attributes(g, 'thick')
            m_th = np.mean(np.asarray(thick.values()))
            length = nx.get_edge_attributes(g, 'length')
            m_le = np.mean(np.asarray(length.values()))
            number = nx.get_edge_attributes(g, 'number')
            num_dict = {}
            for k in number:
                for v in number[k]:
                    num_dict.setdefault(v, []).append(k)
            for k,l in zip(thick.values(), length.values()):
                l1.append(float(l)*m_th/(m_le * k))
            for k in br.keys():
                for l in num_dict[k]:
                    l2.append(float(length[l] * m_th)/(m_le * thick[l]))
            u+=1

        hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.3))
        hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.3))
        center1 = (bins1[:-1] + bins1[1:])/2
        center2 = (bins2[:-1] + bins2[1:])/2
    #    np.save('/home/yuliya/codes/lengths/' + path + '/total_arm_all_len.npy', center1)    
    #    np.save('/home/yuliya/codes/lengths/' + path + '/total_arm_break_len.npy', center2)
    #    np.save('/home/yuliya/codes/lengths/' + path + '/total_arm_all.npy', hist1/float(len(l1)))
    #    np.save('/home/yuliya/codes/lengths/' + path + '/total_arm_break.npy', hist2/float(len(l1)))
    plt.plot(center1, hist1/float(len(l1)), '.', color='red', label = 'all branches')
    plt.plot(center2, hist2/float(len(l1)), '.', color='blue', label = 'breaking branches')
    plt.legend()
    plt.xlabel('l / d', fontsize=18)
    plt.ylabel('P(l/d)', fontsize = 18)

    """
    Probability as a function of aspect ratio.
    """
    hist1, bins1 = np.histogram(l1, np.arange(0, max(l1)+1, 0.3))
    hist2, bins2 = np.histogram(l2, np.arange(0, max(l1)+1, 0.3))
    center1 = (bins1[:-1] + bins1[1:])/2
    plt.plot(center1, hist2/hist1.astype('float32'), '.', color='red', label = 'v34')
    plt.legend()
    plt.xlabel('d / l', fontsize = 18)
    plt.ylabel('Breakup Probability', fontsize = 18)


def correl():
###length-thickness correlation
#l1, l2 =list(), list()
#t1, t2 = list(), list()
#for i,j in zip(d[0], d1[0]):
#    g = nx.read_gpickle(i)
#    br = np.load(j).item()
#    thick = nx.get_edge_attributes(g, 'thick')
#    length = nx.get_edge_attributes(g, 'length')
#    number = nx.get_edge_attributes(g, 'number')
#    num_dict = {}
#    for k in number:
#        for v in number[k]:
#            num_dict.setdefault(v, []).append(k)
#    for k,l in zip(thick.values(), length.values()):
#        l1.append(l/float(ellv34[u]))
#        t1.append(k/float(ellv34[u]))
#    for k in br.keys():
#        for l in num_dict[k]:
#            l2.append(length[l]/float(ellv34[u]))
#            t2.append(thick[l]/float(ellv34[u]))
#    u+=1
#
#a1 = np.append([np.asarray(l1)], [np.asarray(t1)], axis=0)
#a2 = np.append([np.asarray(l2)], [np.asarray(t2)], axis=0)    
#cor_all = np.corrcoef(a1)
#cor_br = np.corrcoef(a2)


######Number of breakups
#for i,j in zip(data1, data):
#    br = np.load(i).item()
#    n_of_br.append(len(br.keys()))
#    g = nx.read_gpickle(j)
#    n_tot = np.append(n_tot, nx.number_of_edges(g))
    
########length histogram (loops_brannches - input graphs) 
#u = 0
#for i, j, c, c1 in zip(data[9:-1][::10], data1[1:][::10], colors, colors1):
#    g = nx.read_gpickle(i)
#    br = np.load(j).item()
#    length = nx.get_edge_attributes(g, 'length')
#    number = nx.get_edge_attributes(g, 'number')
#    num_dict = {}
#    for k in number:
#        for v in number[k]:
#            num_dict.setdefault(v, []).append(k)
#    l1= list()
#    l2 = list()
#    l2.append(0)
#    for k in length.values():
#        l1.append(k/float(ellv34[u]))
#    for k in br.keys():
#        for l in num_dict[k]:
#            l2.append(length[l]/float(ellv34[u]))
#
#    hist1,bins1 = np.histogram(l1, np.arange(0,max(l1)+1, 0.6))
#    center1 = (bins1[:-1] + bins1[1:])/2
#    cum1 = np.cumsum(hist1)
#    hist2,bins2 = np.histogram(l2[1:], np.arange(0,max(l2)+1, 0.6))
#    center2 = (bins2[1:] + bins2[:-1])/2
#    cum2 = np.cumsum(hist2)
#    plt.plot(center1, cum1/float(len(l1)), '--', color=c, label = j[32:-6])
##    plt.plot(center2, cum2/float(len(l2[1:])), '--', color=c1, label = j[32:-6])
#    plt.legend()
#    plt.xlabel('l/l_typ', fontsize=18)
#    plt.ylabel('P(l/l_typ)', fontsize=18)
#    u+=1


#plt.plot(ellv30[:-2], n_of_br/n_tot.astype('float32'), '.', color='blue', label='#breakups/vol')
#plt.plot(ellv34, -eu, '.', color='red', label = '- euler char./vol')    
#plt.legend()
#plt.xlabel('length', fontsize=18)


########lthickness histogram (loops_brannches - input graphs) 
#u = 0
#for i, j, c, c1 in zip(data[9:-1][::10], data1[1:][::10], colors, colors1):
#    g = nx.read_gpickle(i)
#    br = np.load(j).item()
#    thick = nx.get_edge_attributes(g, 'thick')
#    number = nx.get_edge_attributes(g, 'number')
#    num_dict = {}
#    for k in number:
#        for v in number[k]:
#            num_dict.setdefault(v, []).append(k)
#    t1= list()
#    t2 = list()
#    t2.append(0)
#    for k in thick.values():
#        t1.append(k/float(ellv34[u]))
#    for k in br.keys():
#        for l in num_dict[k]:
#            t2.append(thick[l]/float(ellv34[u]))
#
#    hist1,bins1 = np.histogram(t1, np.arange(0,max(t1)+1, 0.06))
#    center1 = (bins1[:-1] + bins1[1:])/2
#    cum1 = np.cumsum(hist1)
#    hist2,bins2 = np.histogram(t2[1:], np.arange(0,max(t2)+1, 0.06))
#    center2 = (bins2[1:] + bins2[:-1])/2
#    cum2 = np.cumsum(hist2)
#    plt.plot(center1, cum1/float(len(t1)), '--', color=c, label = j[32:-6])
#    plt.plot(center2, cum2/float(len(t1)), '--', color=c1, label = j[32:-6])
##    plt.legend()
#    plt.xlabel('d/l_typ', fontsize=18)
#    plt.ylabel('P(d/l_typ)', fontsize=18)
#    u+=1

#######aspect ratio histogram (loops_brannches - input graphs) 
#u = 0
#for i, j, c, c1 in zip(data[9:-1][::10], data1[1:][::10], colors, colors1):
#    g = nx.read_gpickle(i)
#    br = np.load(j).item()
#    thick = nx.get_edge_attributes(g, 'thick')
#    length = nx.get_edge_attributes(g, 'length')
#    number = nx.get_edge_attributes(g, 'number')
#    num_dict = {}
#    for k in number:
#        for v in number[k]:
#            num_dict.setdefault(v, []).append(k)
#    a1= list()
#    a2 = list()
#    a2.append(0)
#    for k, l in zip(thick.values(), length.values()):
#        a1.append(l/float(k))
#    for k in br.keys():
#        for l in num_dict[k]:
#            a2.append(length[l]/float(thick[l]))
#
#    hist1,bins1 = np.histogram(a1, np.arange(0,max(a1)+1, 0.6))
#    center1 = (bins1[:-1] + bins1[1:])/2
#    cum1 = np.cumsum(hist1)
#    hist2,bins2 = np.histogram(a2[1:], np.arange(0,max(a2)+1, 1))
#    center2 = (bins2[1:] + bins2[:-1])/2
#    cum2 = np.cumsum(hist2)
##    plt.plot(center1, cum1/float(len(a1)), '--', color=c, label = j[32:-6])
#    plt.plot(center2, cum2/float(len(a2[1:])), '--', color=c1, label = j[32:-6])
#    plt.legend()
#    plt.xlabel('l / d', fontsize=18)
#    plt.ylabel('P(l/d)', fontsize=18)
#    u+=1



######distance measure - euclid
#plt.figure(2)
#for g, c in zip(data2[10:][::3], colors):
#    pos = np.load(g)
#    d = list()
#    d.append(0)
#    for comb in itertools.combinations(pos, 2):  
#        i = comb[0]
#        j = comb[1]
#        dist = np.linalg.norm(i - j)
#        d.append(dist)
#    hist, bins = np.histogram(d[1:], np.arange(0, round(max(d)/100) * 100 + 1, 100))
#    center = (bins[1:] + bins[:-1])/2
#    plt.plot(center/float(ell[u]), hist/float(len(d[1:])), '-', color=c)
#    u+=1    

######distance measure - skeleton
#@profile
#def dist_skel():
##    plt.figure(1)
#    d_skel = list()
#    u=0
#    for b, g, c in zip(data1[52:53], data[60:61], colors):
#        br = np.load(b).item().keys()
#        gr = nx.read_gpickle(g)
#        number = nx.get_edge_attributes(gr, 'number')
#        num_dict = {}
#        for k in number:
#            for v in number[k]:
#                num_dict.setdefault(v, []).append(k)
#        for comb in itertools.combinations(br, 2):
#            n1 = num_dict[comb[0]][0][0]
#            n2 = num_dict[comb[0]][0][1]
#            m1 = num_dict[comb[1]][0][0]
#            m2 = num_dict[comb[1]][0][1]
#            p_min = p_max = nx.shortest_path_length(gr, n1, m1, 'length')
#            for i in itertools.product((n1,n2), (m1,m2)):
#                p = nx.shortest_path_length(gr, i[0], i[1], 'length')
#                if p < p_min:
#                    p_min = p
#                if p > p_max:
#                    p_max = p
#            d_skel.append((p_min+p_max)/(2. * ellv34[u]))
#    hist, bins = np.histogram(d_skel, np.arange(0, max(d_skel)+ 1), 0.01)
#    center = (bins[1:] + bins[:-1])/2
#    plt.plot(center/float(ellv34), hist/float(len(d_skel)), '-', label='breakup distances')
#        u+=1
#    return d_skel
    




    
###picking random breakups
def dist_skel(data1, data, el):
    """
    Chooses random pairs of break-ups and computes the distances
    between them.
    
    Parameters
    ------------
    data1 : list of pathes for the break-ups dictionaries
    data : list of pathes for the graphs
    el : list of length scales
    
    Return
    -------
    d_skel : list of distances along the skeleton.
    """
#    plt.figure(1)
    d_skel = list()
    br_peak = list()
    u=0
    for b, g in zip(data1[:], data[:]):
        br = np.load(b).item().keys()
        if br==[]:
            continue
        gr = nx.read_gpickle(g)
        number = nx.get_edge_attributes(gr, 'number')
        num_dict = {}
        for k in number:
            for v in number[k]:
                num_dict.setdefault(v, []).append(k)
        a = len(list(br))
        if a > 1000:
            a = 1000
        for j in range(a):
            e1 = random.choice(br)
            e2 = random.choice(br)
            if (e1==e2) or (num_dict[e1]==num_dict[e2]):
                continue
            comb = (e1, e2)
            n1 = num_dict[comb[0]][0][0]
            n2 = num_dict[comb[0]][0][1]
            m1 = num_dict[comb[1]][0][0]
            m2 = num_dict[comb[1]][0][1]
            try:            
                p_min = p_max = nx.shortest_path_length(gr, n1, m1, 'length')
            except nx.NetworkXNoPath:
                continue
            for i in itertools.product((n1,n2), (m1,m2)):
                p = nx.shortest_path_length(gr, i[0], i[1], 'length')
                if p < p_min:
                    p_min = p
                if p > p_max:
                    p_max = p
            d = (p_min+p_max)/(2. * el[u])
            if d<2:
                br_peak.append((b, e1, e2))
            d_skel.append(d)
#    hist, bins = np.histogram(d_skel, np.arange(0, max(d_skel)+ 1), 0.01)
#    center = (bins[1:] + bins[:-1])/2
#    plt.plot(center/float(ellv34), hist/float(len(d_skel)), '-', label='breakup distances')
        u+=1
    return d_skel

def random_pairs_dist(data1, data, el):
    """
    Chooses random pairs of branches in the skeleton and computes the distances
    between them.
    
    Parameters
    ------------
    data1 : list of pathes for the break-ups dictionaries. Used here for the 
            estimation of the necessary number of pairs.
    data : list of pathes for the graphs
    el : list of length scales
    
    Return
    -------
    d_ran : list of distances along the skeleton between random branches.
    """
    d_ran = list()
    u=0
    for b, g in zip(data1, data):
        br = np.load(b).item().keys()
        gr = nx.read_gpickle(g)
        edg = list(nx.edges(gr))
        number = nx.get_edge_attributes(gr, 'number')
        num_dict = {}
        for k in number:
            for v in number[k]:
                num_dict.setdefault(v, []).append(k)
        a = len(list(br))
        if a > 1500: #We cannot use here infinitly many pairs because of the computing time.
            a = 1500
        for j in range(a):
            e1 = random.choice(edg)
            e2 = random.choice(edg)
            if (e1==e2):
                continue
            n1 = e1[0]
            n2 = e1[1]
            m1 = e2[0]
            m2 = e2[1]
            try: 
                p_min = p_max = nx.shortest_path_length(gr, n1, m1, 'length')
            except nx.NetworkXNoPath:
                continue    
            for i in itertools.product((n1,n2), (m1,m2)):
                p = nx.shortest_path_length(gr, i[0], i[1], 'length')
                if p < p_min:
                    p_min = p
                if p > p_max:
                    p_max = p
            d = (p_min+p_max)/(2. * el[u])
            d_ran.append(d)
        u+=1
    return d_ran

def distance_correl():
    """
    Computation of distribution of distances between random pairs of branches
    and random break-ups.
    
    """
    d_skel1 = list()
    d_ran1 = list()
    #d_skel_peak = dist_skel()
    for i, j, k in zip(d1, d, ell):
        d_skel1.append(np.asarray(dist_skel(i, j, k)))
        d_ran1.append(np.asarray(random_pairs_dist(i, j, k)))

    pa = ('v34', 'vsi05', 'v30', 'v35')
    for i, j, path in zip(d_ran1, d_skel1, pa):
        hist1, bins1 = np.histogram(j, np.arange(0, max(j)+ 1, 0.5))
        center1 = (bins1[1:] + bins1[:-1])/2
        hist2, bins2 = np.histogram(i, np.arange(0, max(i)+ 1, 0.5))
        center2 = (bins2[1:] + bins2[:-1])/2
    #    np.save('/home/yuliya/codes/lengths/' + path + '/dist_distr_array_all.npy', i)    
    #    np.save('/home/yuliya/codes/lengths/' + path + '/dist_distr_array_break.npy', j)
    #    np.save('/home/yuliya/codes/lengths/' + path + '/dist_distr_hist_all.npy', hist2/float(len(i)))
    #    np.save('/home/yuliya/codes/lengths/' + path + '/dist_distr_bins_all.npy', center2)
    #    np.save('/home/yuliya/codes/lengths/' + path + '/dist_distr_hist_break.npy', hist1/float(len(j)))
    #    np.save('/home/yuliya/codes/lengths/' + path + '/dist_distr_bins_break.npy', center1)
        
    plt.figure(figsize=(9, 7))
    plt.plot(center1, hist1/float(len(j)), '--', label='breakup distances')
    plt.plot(center2, hist2/float(len(i)), '-', label='random distances')
    plt.legend()
    plt.xlabel('dist/l_typ', fontsize=18)
    plt.ylabel('P(dist/l_typ)', fontsize=18)

#####probability to break
def prob_break():
    u = 0
    ratio1 = []
    for i, j in zip(data, data1):

        gr = nx.read_gpickle(i)
        br = np.load(j).item().keys()
    
        number = nx.get_edge_attributes(gr, 'number')
        num_dict = {}    
        for k in number:
            for v in number[k]:
                num_dict.setdefault(v, []).append(k)
        phys_br = list()
        for k in br:
            phys_br.append(number[num_dict[k][0]][0])
        
        if len(br)<>len(phys_br):
            azjkhrg

        events = len(phys_br)    

        branches = nx.number_of_edges(gr)
        
        ratio1 = np.append(ratio1, float(events)/branches)
        u+=1
    return ratio1

def topology(data, ell):
    """
    Computation of topological characteristics.
    
    Parameters
    ------------
    data : array of pathes to the graphs
    ell : list of length scales
    """    
    for i in data:
        G = nx.read_gpickle(i)
        B = nx.number_of_edges(G)
        V = nx.number_of_nodes(G)
        Euler = V - B
        C = (B-V)/float(V)
        eu.append(Euler)
        c_t.append(C)
        vert.append(V)
        bran.append(B)
        plt.figure(2)
        plt.plot(t[k], C, '.', color='red')
        plt.xlabel('t, s', fontsize=18)
        plt.ylabel('log C_t, (B - V)/V', fontsize=18)
        plt.yscale('log')
        plt.ylabel('Euler char., V - B', fontsize=18)
        k += 1

    plt.plot(ell, c_t, '.', label='v23')
    #
    #np.save('/backup/yuliya/v23/graphs_largedom/Euler.npy', eu)
    #np.save('/backup/yuliya/v23/graphs_largedom/C_t.npy', c_t)
    #np.save('/backup/yuliya/v23/graphs_largedom/V.npy', vert)
    #np.save('/backup/yuliya/v23/graphs_largedom/B.npy', bran)
    #np.save('/backup/yuliya/vsi01/graphs_largdom/time.npv23/graphs_largedom/y', t)
    plt.yscale('log')
    plt.xscale('log')   
    #plt.savefig('/backup/yuliya/v34/graphs/v34_logC_t.png', dpi=500)
    plt.close()    
    
    
    



    