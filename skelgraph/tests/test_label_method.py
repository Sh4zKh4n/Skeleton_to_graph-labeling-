from skelgraph.label_method import numb
import numpy as np
from numpy.testing import assert_equal, assert_raises


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
