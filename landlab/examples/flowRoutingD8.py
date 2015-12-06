#! /usr/env/python
"""

Python implementation of d8 routing scheme for
Landlab with a rectangular, uniform, mesh.

Some issues that might come up ...
This class relies on a method in model_grid called
find_node_in_direction_of_max_slope.  That method does not do
any boundary checking.  It needs to be done somewhere, maybe here,
maybe there.  Alternatively, if boundary elevations are always set
in a consistent way depending on the type of bounary, maybe no boundary
checking is needed?

Last updated NG 6/2013

"""

#from landlab.model_grid import RasterModelGrid
from numpy import *


class RouteFlowD8(object):
    """
    This class finds the steepest path among 8 possible directions, so
    diagonals are considered.
    The class assumes that the model is using a rectangular, uniform (raster)
    grid.
    """

    def __init__(self, num_cells):
        """
        This sets the num_cells parameter.
        This class assumes that the number of cells does not change after a
        class item has been instantiated.
        """

        self.num_cells = num_cells
        self.initialize()

        # print 'RouteFlowD8.__init__'

    def initialize(self):
        """
        This sets up the flow direction vector.
        It is initialized to -1 for all nodes.
        A -1 flowdirs value indicates a boundary node.
        Is this method really needed?  Can we just put this in __init__?
        """

        self.flowdirs = -ones(self.num_cells, dtype=int)

    def calc_flowdirs(self, mg, z):
        """
        This assigns the flow directions using the function
        find_node_in_direction_of_max_slope in the model_grid.
        The flowdirs vector contains the node id that a node flows to.
        If the node is a boundary node, the flowdirs vector has a value of -1.

        Method inputs: the model grid and elevation vector
        Method returns: the flow direction vector
        """

        for i in range(0, self.num_cells):
            if mg.is_interior(i):
                self.flowdirs[i] = mg.find_node_in_direction_of_max_slope(z, i)

        return self.flowdirs
