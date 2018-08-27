# coding: utf8
#! /usr/env/python
"""
"""

from landlab.components.profiler.base_profiler import Profiler


class RidgeProfiler(Profiler):
    """
    """
    def __init__(self, grid):
        """
        """
        super(RidgeProfiler, self).__init__(grid)
        self._grid = grid
    
    def something(self):
        """
        """
        pass
