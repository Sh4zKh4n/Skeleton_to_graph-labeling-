from skelgraph.label_method import numb, label_nodes, label_br, neigh_br
import numpy as np
from numpy.testing import assert_equal, assert_raises
import networkx as nx

def test_numb():
    a = np.ones((2, 2, 2))
    assert np.all(numb(a) == 7) # only corners of the cube
    a = np.zeros((3, 3, 3))
    assert np.all(numb(a) == 0)
    a[1, 1, 1] = 1
    assert np.all(numb(a) == 0)


def test_numb_dtype():
    a = np.ones((2, 2, 2), dtype=np.bool)
    assert np.all(numb(a) == 7) # only corners of the cube


def test_neigh_br():
    num = np.array([[[0, 0, 0],
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
                     [0, 0, 0]]])
    arr_nodes = label_nodes(num)
    arr_br = label_br(num)
    G = neigh_br(arr_nodes, arr_br)
    assert np.all(nx.get_node_attributes(G, 'neigh')[1] == [1, 2, 3])
