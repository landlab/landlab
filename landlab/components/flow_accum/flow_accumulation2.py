#! /usr/env/python
"""The landlab flow accumulation module.

A python flow accumulation module. It is designed to be general, and to
operate across multiple grids and multiple flow direction patterns. However,
at the moment, only a steepest descent (single path) routing scheme is
implemented.

Notes
-----
There remain some outstanding issues with the handling of boundary cells,
which this component has inherited from `flow_routing_D81.
"""
# Created DEJH, 8/2013
from __future__ import print_function

from six.moves import range

import warnings

import numpy as np

from landlab import CLOSED_BOUNDARY


class AccumFlow(object):

    """Flow accumulation module.

    This class allows the routing of flow around a landscape according to a
    previously calculated flow direction vector. It is not sensitive to grid
    type. It will eventually be able to work with discharges which are split
    across more than one node, but at the moment, assumes a single line of
    descent for a given node.
    """

    def __init__(self, grid):
        self._grid = grid

        # Prefilled with zeros, size of WHOLE grid + 1, to allow -1 ids
        self._flow_accum_by_area = np.empty(grid.number_of_nodes + 1,
                                            dtype=float)
        self._flow_accum_by_area[-1] = 0.

    @property
    def grid(self):
        return self._grid

    def calc_flowacc(self, z, flowdirs):
        """Calculate flow accumulation.

        Parameters
        ----------
        z : ndarray
            Elevations at nodes.
        flowdirs : ndarray
            Downstream neighbor node.

        Returns
        -------
        ndarray
            Total upstream area at each node.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> grid = RasterModelGrid((4, 5), spacing=(3, 4))
        >>> z = np.array([3., 3., 0., 3., 3.,
        ...               3., 2., 1., 2., 3.,
        ...               3., 2., 2., 2., 3.,
        ...               3., 3., 3., 3., 3.])

        First calculate the flow directions.

        >>> from landlab.components.flow_routing import RouteFlowD8
        >>> flow_router = RouteFlowD8(len(z))
        >>> flowdirs, _ = flow_router.calc_flowdirs(grid ,z)
        >>> flowdirs
        array([ 6,  7, -1,  7,  8,
                6,  2,  2,  2,  8,
               11,  7,  7,  7, 13,
               11, 11, 12, 13, 13])

        Then calculate the drainage areas.

        >>> from landlab.components.flow_accum import AccumFlow
        >>> accumulator = AccumFlow(grid)
        >>> accumulator.calc_flowacc(z, flowdirs)
        array([  0.,   0.,  72.,   0.,   0.,
                 0.,  12.,  48.,  12.,   0.,
                 0.,  12.,  12.,  12.,   0.,
                 0.,   0.,   0.,   0.,   0.])
        """
        nodes = np.where(self.grid.status_at_node != CLOSED_BOUNDARY)[0]

        # Cell areas for nodes without cells is 0.
        self._flow_accum_by_area[nodes] = self.grid.cell_area_at_node[nodes]

        # Perform test to see if the flowdir data is a single vector, or
        # multidimensional, here. Several ways possible: 1. Is the vector
        # multidimensional?, e.g., try: data.flowdirs.shape[1] 2. set a flag
        # in flowdir.

        try:
            # Descending order.
            height_order_nodes = np.argsort(z[nodes])[::-1]
        except:
            warnings.warn("Cells could not be sorted by elevation. Does the "
                          "data object contain the elevation vector?")

        try:
            sorted_flowdirs = flowdirs[nodes][height_order_nodes]
        except:
            warnings.warn("Flow directions could not be sorted by elevation. "
                          "Does the data object contain the flow direction "
                          "vector?")

        ## Inefficient Python code.
        for i in range(len(sorted_flowdirs)):
            iter_height_order = height_order_nodes[i]
            iter_sorted_fldirs = sorted_flowdirs[i]
            self._flow_accum_by_area[iter_sorted_fldirs] += (
                self._flow_accum_by_area[nodes][iter_height_order])

        return self._flow_accum_by_area[:-1]
