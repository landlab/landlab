# coding: utf8
#! /usr/env/python
"""
"""
import numpy as np

from landlab.components.profiler.base_profiler import _NetworkProfiler


class RidgeProfiler(_NetworkProfiler):
    """
    """
    def __init__(self, grid, threshold=None):
        """
        """
        super(RidgeProfiler, self).__init__(grid)
        self._grid = grid

        if threshold is None:
            threshold = 4. * np.amin(grid.area_of_cell)
        self.threshold = threshold

        self._ridge_mask = self._drainage_area < self.threshold


    def something(self):
        """
        """
        pass
