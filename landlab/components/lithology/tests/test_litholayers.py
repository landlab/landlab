#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 09:17:36 2018

@author: barnhark
"""

import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.bmi import wrap_as_bmi
from landlab.components import LithoLayers


def test_litholayers_as_bmi():
    """Test Litholayers can be wrapped with a BMI."""
    wrap_as_bmi(LithoLayers)


def test_z0s_ids_different_shape():
    """Test that providing z0s and ids of different shapes raises an error."""
    mg = RasterModelGrid((3, 3))
    z0s = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        LithoLayers(mg, z0s, ids, attrs)


def test_z0s_bad_order():
    """Test that providing z0s in a bad order raises an error."""
    mg = RasterModelGrid((3, 3))
    z0s = [-4, -3, -2, -1, 0, 1, 2, 6, 4]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        LithoLayers(mg, z0s, ids, attrs)


def test_bad_function():
    """Test that providing a function of three variables."""
    mg = RasterModelGrid((3, 3))
    z0s = [-4, -3, -2, -1, 0, 1, 2, 4, 6]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        LithoLayers(mg, z0s, ids, attrs, function=lambda x, y, z: 0 * x + 0 * y + z)


def test_function_returns_scalar():
    mg = RasterModelGrid((3, 3))
    z0s = [-4, -3, -2, -1, 0, 1, 2, 4, 6]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        LithoLayers(mg, z0s, ids, attrs, function=lambda x, y: 1.0)


def test_function_returns_wrong_number_of_values():
    mg = RasterModelGrid((3, 3))
    z0s = [-4, -3, -2, -1, 0, 1, 2, 4, 6]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        LithoLayers(mg, z0s, ids, attrs, function=lambda x, y: np.array([1.0]))
