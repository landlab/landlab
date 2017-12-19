#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:39:32 2017

@author: KRB
"""

from landlab import RasterModelGrid
from landlab.components import (TaylorNonLinearDiffuser)
import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_raises

import warnings

def test_raise_stability_error():
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.node_x.copy()**2
    Cdiff = TaylorNonLinearDiffuser(mg)
    assert_raises(RuntimeError, Cdiff.soilflux, 10, if_unstable='raise')

def test_raise_kwargs_error():
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.node_x.copy()**2
    assert_raises(TypeError, TaylorNonLinearDiffuser, mg, bad_name='true')

def test_infinite_taylor_error():
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.node_x.copy()**4
    Cdiff = TaylorNonLinearDiffuser(mg, nterms=400)
    assert_raises(RuntimeError, Cdiff.soilflux, 10)

#def test_warn():
#    mg = RasterModelGrid((5, 5))
#    z = mg.add_zeros('node', 'topographic__elevation')
#    z += mg.node_x.copy()**2
#    Cdiff = TaylorNonLinearDiffuser(mg)
#    
#    with warnings.catch_warnings(record=True) as w:
#        # Cause all warnings to always be triggered.
#        warnings.simplefilter("always")
#        # Trigger a warning.
#        Cdiff.soilflux(dt=10, if_unstable='warn')
#        # Verify some things
#        assert len(w) == 1
#        assert issubclass(w[-1].category, RuntimeWarning)

