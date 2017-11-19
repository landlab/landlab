#! /usr/bin/env python
"""Calculate steepest descent on a raster grid.

Steepest descent functions for raster grids
+++++++++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

"""
import numpy as np

from landlab.core.utils import make_optional_arg_into_id_array
from landlab.utils.decorators import deprecated
from landlab.grid.raster_gradients import (
    calc_grad_across_cell_faces,
    calc_grad_across_cell_corners,
    calc_grad_along_node_links)


_VALID_ROUTING_METHODS = set(['d8', 'd4'])


def _assert_valid_routing_method(method):
    """Check if name is a valid routing method.

    Parameters
    ----------
    method : str
        Name of a routing method.

    Raises
    ------
    ValueError
        If the method name is invalid.
    """
    if method not in _VALID_ROUTING_METHODS:
        raise ValueError(
            '%s: routing method not understood. should be one of %s' %
            (method, ', '.join(_VALID_ROUTING_METHODS)))
