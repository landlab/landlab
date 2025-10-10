#!/usr/bin/env python2
"""
Created on Fri Mar  3 10:39:32 2017

@author: gtucker
"""
import pytest
from numpy.testing import assert_equal

from landlab import RasterModelGrid
from landlab.components import DepthDependentDiffuser
from landlab.components import ExponentialWeatherer


def test_raise_kwargs_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("soil__depth", at="node")
    z = mg.add_zeros("topographic__elevation", at="node")
    BRz = mg.add_zeros("bedrock__elevation", at="node")
    z += mg.node_x.copy()
    BRz += mg.node_x / 2.0
    soilTh[:] = z - BRz
    ExponentialWeatherer(mg)
    with pytest.raises(TypeError):
        DepthDependentDiffuser(mg, diffusivity=1)

    DDdiff = DepthDependentDiffuser(mg)

    with pytest.raises(TypeError):
        DDdiff.soilflux(2.0, bad_var=1)


def test_handle_creep_coefficient():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("bedrock__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    # linear_diffusivity is deprecated, but should still assign
    # the correct value
    dd = DepthDependentDiffuser(mg, linear_diffusivity=0.03)
    assert_equal(dd._K, 0.03)

    # if soil_transport_velocity and linear_diffusivity are both
    # given, the latter should be ignored
    dd = DepthDependentDiffuser(
        mg, soil_transport_velocity=0.02, linear_diffusivity=0.03
    )
    assert_equal(dd._K, 0.02)

    # if neither are given, the default should be used
    dd = DepthDependentDiffuser(mg)
    assert_equal(dd._K, 1.0)
