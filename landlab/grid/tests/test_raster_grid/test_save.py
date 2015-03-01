#! /usr/bin/env python
import os

import numpy as np
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_true, assert_equal, assert_raises
try:
    from nose.tools import assert_list_equal
except ImportError:
    from landlab.testing.tools import assert_list_equal

from landlab.testing.tools import cdtemp
from landlab.io import write_esri_ascii, read_esri_ascii
from landlab import RasterModelGrid

def test_save():
    grid = RasterModelGrid(4, 5, 2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))
    with cdtemp() as _:
        grid.save('test.asc')
