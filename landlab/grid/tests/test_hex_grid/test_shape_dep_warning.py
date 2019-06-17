# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 10:36:03 2015

@author: gtucker
"""
import warnings

from landlab import HexModelGrid


def test_shape_dep_warning():
    """Test dep warning on use of shape keyword instead of node_layout."""
    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger the deprecation warning.
        HexModelGrid(3, 2, shape="rect")
        # Verify some things
        assert len(w) == 1
        assert issubclass(w[-1].category, DeprecationWarning)
        assert "node_layout" in str(w[-1].message)
