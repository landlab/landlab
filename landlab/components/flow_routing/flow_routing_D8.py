#! /usr/env/python
"""

Python implementation of d8 routing scheme for
Landlab with a rectangular, uniform, mesh.

Some issues that might come up ...
This class relies on a method in model_grid called
find_node_in_direction_of_max_slope.  That method does not do
any boundary checking.  It needs to be done somewhere, maybe here,
maybe there?  Alternatively, if boundary elevations are always set
in a consistent way depending on the type of bounary, maybe no boundary
checking is needed?

Last updated NG 8/2013

"""

#from landlab.model_grid import RasterModelGrid
from numpy import *


class RouteFlowD8(object):
    """
    This class finds the steepest path among 8 possible directions, so
    diagonals are considered.
    The class assumes that the model is using a rectangular, uniform (raster)
    grid.

    This sets the num_cells parameter.
    This class assumes that the number of cells does not change after a
    class item has been instantiated.
    """

    def __init__(self, num_nodes):
        """
        This sets the num_cells parameter.
        This class assumes that the number of cells does not change after a
        class item has been instantiated.
        """

        self.num_nodes = num_nodes
        self.initialize()

        #print 'RouteFlowD8.__init__'

    def initialize(self):
        """
        This sets up the flow direction vector.
        It is initialized to -1 for all nodes.
        A -1 flowdirs value indicates a boundary node.
        """

        self.flowdirs = -ones(self.num_nodes, dtype=int)


    def calc_flowdirs(self, mg, z):
        """
        This assigns the flow directions using the function
        find_node_in_direction_of_max_slope in the model_grid.
        The flowdirs vector contains the node id that a node flows to.
        If the node is a boundary node, the flowdirs vector has a value of -1.

        This now only applies to interior nodes.
        However, should a fixed value node, which counts as an interior node,
        have a flow direction?  Still not sure about this.

        Method inputs: the model grid and elevation vector
        Method returns: the flow direction vector, and as of DEJH modifications
        Sept 2013, the maximum (most steeply downhill) slope leaving each node.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing import RouteFlowD8

        >>> grid = RasterModelGrid((4, 5), spacing=(1, 1))
        >>> z = np.array([3., 3., 0., 3., 3.,
        ...               3., 2., 1., 2., 3.,
        ...               3., 2., 2., 2., 3.,
        ...               3., 3., 3., 3., 3.])

        First calculate the flow directions.

        >>> flow_router = RouteFlowD8(len(z))
        >>> flowdirs, _ = flow_router.calc_flowdirs(grid ,z)
        >>> flowdirs
        array([ 6,  7, -1,  7,  8,
                6,  2,  2,  2,  8,
               11,  7,  7,  7, 13,
               11, 11, 12, 13, 13])
        """
        #for i in range(0, self.num_nodes):
        #    if mg.is_interior(i):
        #        self.flowdirs[i] = mg.find_node_in_direction_of_max_slope(z, i)

        #Alt. method, DEJH:
        self.gradients_on_active_links = mg.calc_grad_at_active_link(z)
        max_slopes, self.flowdirs = mg._calculate_steepest_descent_on_nodes(z, self.gradients_on_active_links, dstr_node_ids=self.flowdirs)

        return self.flowdirs, max_slopes
