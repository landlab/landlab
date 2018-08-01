#!/usr/env/python

"""
lake_filler_with_flux.py: Component to fill depressions in a landscape while
honouring mass balance.

Similar to the DepressionFinderAndRouter, but will not fill a lake to the brim
if there is not enough incoming flux to do so. Designed to "play nice" with
the FlowAccumulator.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message

from landlab import BAD_INDEX_VALUE
import six
import numpy as np


class LakeEvaporator(Component):
    """
    Component to take a water surface overlying a bedrock surface as a series
    of (specified) existing lakes, then reduce the volume of water in each
    lake (i.e., drop the water surface) according to a loss term specified
    at each lake node.
    """
    def __init__(self):
        """
        Initialise the component.
        """
        pass


class LakeFillerWithFlux(Component):
    """
    """

    def run_one_step(dt):
        """
        """
        # First, we need a conventional lake fill. This involves calling the
        # FillSinksBarnes component. Nice, FAST special cases where lakes fill
        # up or only have one pit, which we can use to accelerate things,
        # so this is very much worth it.
        pass
