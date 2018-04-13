#! /usr/bin/env python
"""

"""
import numpy as np

from ..graph import Graph
from ..field import GraphFields

from ..utils.decorators import cache_result_in_object
from .decorators import return_readonly_id_array
from landlab.utils.decorators import make_return_array_immutable

from .nodestatus import (CORE_NODE, FIXED_VALUE_BOUNDARY,
                         FIXED_GRADIENT_BOUNDARY, LOOPED_BOUNDARY,
                         CLOSED_BOUNDARY)
from .linkstatus import ACTIVE_LINK, FIXED_LINK, INACTIVE_LINK
from .linkstatus import set_status_at_link
from ..core.utils import add_module_functions_to_class

class NetworkModelGrid(Graph, GraphFields):
    """A ModelGrid of just nodes and links.

    Parameters
    ----------
    yx_of_node : tuple of ndarray
        Node y and x coordinates.
    links : array of tuple of int
        Nodes at link tail and head.

    Examples
    --------
    >>> from landlab.grid.network import NetworkModelGrid
    >>> y_of_node = (0, 1, 2, 2)
    >>> x_of_node = (0, 0, -1, 1)
    >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.x_of_node
    array([ 0.,  0., -1.,  1.])
    >>> grid.y_of_node
    array([ 0.,  1.,  2.,  2.])
    >>> grid.nodes_at_link
    array([[0, 1],
           [2, 1],
           [1, 3]])
    """

    def __init__(self, yx_of_node, links, **kwds):
        Graph.__init__(self, yx_of_node, links=links)
        GraphFields.__init__(self,
                             {'node': self.number_of_nodes,
                              'link': self.number_of_links,
                              'grid': 1},
                             default_group='node')

        self._node_status = np.zeros(self.number_of_nodes, dtype=np.uint8)
        self.bc_set_code = 0


    @property
    #@override_array_setitem_and_reset('reset_status_at_node') # this is in BASE, not sure if we need it.
    def status_at_node(self):
        """Get array of the boundary status for each node.

        Examples
        --------
        >>> TODO
        """
        return self._node_status

    @status_at_node.setter
    def status_at_node(self, new_status):
        """Set the array of node boundary statuses."""
        self._node_status[:] = new_status[:]
        self.reset_status_at_node()

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def status_at_link(self):
        """Get array of the status of all links.

        Examples
        --------
        Examples
        --------
        >>> # TODO
        """
        return set_status_at_link(self.status_at_node[self.nodes_at_link])

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def active_links(self):
        """Get array of active links.

        Examples
        --------
        >>> # TODO
        """
        return np.where(self.status_at_link == ACTIVE_LINK)[0]

add_module_functions_to_class(NetworkModelGrid, 'mappers.py', pattern='map_*')
add_module_functions_to_class(NetworkModelGrid, 'gradients.py', pattern='calc_*')
add_module_functions_to_class(NetworkModelGrid, 'divergence.py', pattern='calc_*')
