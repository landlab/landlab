#!/usr/bin/env python2
"""
Created on Fri Mar  3 10:39:32 2017

@author: gtucker
"""
import pytest

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
