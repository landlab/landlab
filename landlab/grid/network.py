#! /usr/bin/env python
"""
A class used to create and manage network models in 2D.
"""
import numpy as np

from landlab.utils.decorators import make_return_array_immutable

from ..core.utils import add_module_functions_to_class
from ..field import GraphFields
from ..graph import NetworkGraph
from ..utils.decorators import cache_result_in_object
from .decorators import override_array_setitem_and_reset, return_readonly_id_array
from .linkstatus import ACTIVE_LINK, set_status_at_link


class NetworkModelGrid(NetworkGraph, GraphFields):
    """Create a ModelGrid of just nodes and links.

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
        NetworkGraph.__init__(self, yx_of_node, links=links)
        GraphFields.__init__(
            self,
            {"node": self.number_of_nodes, "link": self.number_of_links, "grid": 1},
            default_group="node",
        )
        self._node_status = np.zeros(self.number_of_nodes, dtype=np.uint8)
        self.bc_set_code = 0

    @property
    @override_array_setitem_and_reset("reset_status_at_node")
    def status_at_node(self):
        """Get array of the boundary status for each node.

        Examples
        --------
        >>> from landlab.grid.network import NetworkModelGrid
        >>> from landlab import CLOSED_BOUNDARY, CORE_NODE

        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
        >>> grid.status_at_node
        array([0, 0, 0, 0], dtype=uint8)
        >>> grid.status_at_link
        array([0, 0, 0], dtype=uint8)

        Now we change the status at node 0 to a closed boundary. This will
        result in changing the link status as well.

        >>> grid.status_at_node = [CLOSED_BOUNDARY, CORE_NODE, CORE_NODE, CORE_NODE]
        >>> grid.status_at_node
        array([4, 0, 0, 0], dtype=uint8)
        >>> grid.status_at_link
        array([4, 0, 0], dtype=uint8)

        LLCATS: NINF BC
        """
        return self._node_status

    @status_at_node.setter
    def status_at_node(self, new_status):
        """Set the array of node boundary statuses."""
        self._node_status[:] = new_status[:]
        self.reset_status_at_node()

    def reset_status_at_node(self):
        attrs = [
            "_active_link_dirs_at_node",
            "_status_at_link",
            "_active_links",
            "_fixed_links",
            "_activelink_fromnode",
            "_activelink_tonode",
            "_fixed_links",
            "_active_adjacent_nodes_at_node",
            "_fixed_value_boundary_nodes",
            "_link_status_at_node",
        ]
        for attr in attrs:
            try:
                del self.__dict__[attr]
            except KeyError:
                pass

        self.bc_set_code += 1

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def status_at_link(self):
        """Get array of the status of all links.

        Examples
        --------
        >>> from landlab.grid.network import NetworkModelGrid
        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
        >>> grid.status_at_link
        array([0, 0, 0], dtype=uint8)

        LLCATS: LINF BC
        """
        return set_status_at_link(self.status_at_node[self.nodes_at_link])

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def active_links(self):
        """Get array of active links.

        Examples
        --------
        >>> from landlab.grid.network import NetworkModelGrid
        >>> from landlab import FIXED_LINK

        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
        >>> grid.active_links
        array([0, 1, 2])

        LLCATS: NINF BC SUBSET
        """
        return np.where(self.status_at_link == ACTIVE_LINK)[0]


# add only the correct functions
add_module_functions_to_class(
    NetworkModelGrid, "mappers.py", pattern="map_*", exclude="cell|patch"
)
add_module_functions_to_class(
    NetworkModelGrid, "gradients.py", pattern="calc_grad_at_link"
)
