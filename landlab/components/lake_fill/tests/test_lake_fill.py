#!/usr/env/python

import pytest
import numpy as np

from landlab.components import LakeMapperBarnes
from landlab import RasterModelGrid, HexModelGrid, CLOSED_BOUNDARY

"""
These tests test specific aspects of LakeMapperBarnes not picked up in the
various docstrings.
"""

def test_bad_init_method1():
    rmg = RasterModelGrid((5, 5), dx=2.)
    rmg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(rmg, method='Nope')


def test_bad_init_method1():
    rmg = RasterModelGrid((5, 5), dx=2.)
    rmg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(rmg, method='d8')


def test_bad_init_gridmethod():
    hmg = HexModelGrid(30, 29, dx=3.)
    hmg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(hmg, method='D8')

def closed_up_grid():
    mg = RasterModelGrid((5, 5), dx=1.)
    for edge in ('left', 'right', 'top', 'bottom'):
        mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        with pytest.raises(ValueError):
            lmb = LakeMapperBarnes(mg)
