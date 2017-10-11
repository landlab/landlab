#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:39:32 2017

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab.components import (DepthDependentDiffuser,
                                ExponentialWeatherer)
import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_raises

import warnings

def test_raise_kwargs_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros('node', 'soil__depth')
    z = mg.add_zeros('node', 'topographic__elevation')
    BRz = mg.add_zeros('node', 'bedrock__elevation')
    z += mg.node_x.copy()
    BRz += mg.node_x/2.
    soilTh[:] = z - BRz
    expweath = ExponentialWeatherer(mg)
    assert_raises(TypeError, DepthDependentDiffuser, mg, diffusivity=1)
    
    DDdiff = DepthDependentDiffuser(mg)
    
    assert_raises(TypeError, DDdiff.soilflux, 2., bad_var=1)

