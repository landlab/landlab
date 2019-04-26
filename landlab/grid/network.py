#! /usr/bin/env python
"""
A class used to create and manage network models in 2D.
"""
import numpy as np

from landlab.grid.base import _default_axis_names, _default_axis_units
from landlab.utils.decorators import make_return_array_immutable

from ..core import load_params
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
    xy_of_reference : tuple, optional
        Coordinate value in projected space of (0., 0.)
        Default is (0., 0.)

    Examples
    --------
    >>> from landlab import NetworkModelGrid
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

        self.axis_name = kwds.get("axis_name", _default_axis_names(self.ndim))

        self.axis_units = kwds.get("axis_units", _default_axis_units(self.ndim))

        self._ref_coord = tuple(kwds.get("xy_of_reference", (0.0, 0.0)))

    @classmethod
    def from_file(cls, file_like):
        params = load_params(file_like)
        return cls.from_dict(params)

    @classmethod
    def from_dict(cls, params):
        return cls(**params)

    @property
    def xy_of_reference(self):
        """Return the coordinates (x, y) of the reference point.

        For RasterModelGrid and HexModelGrid the reference point is the
        minimum of x_of_node and of y_of_node. By default it is (0, 0). For
        VoronoiDelaunayGrid the reference point is (0, 0). For RadialModelGrid
        it is the (x, y) of the center point.

        The intention of these coordinates is to provide a method to store
        the large float values of projected coordinates.

        Example
        -------

        >>> from landlab import NetworkModelGrid
        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node),
        ...                         nodes_at_link,
        ...                         xy_of_reference = (12345, 678910))
        >>> grid.xy_of_reference
        (12345, 678910)
        >>> grid.xy_of_reference = (98765, 43210)
        >>> grid.xy_of_reference
        (98765, 43210)
        """
        return self._ref_coord

    @xy_of_reference.setter
    def xy_of_reference(self, new_xy_of_reference):
        """Set a new value for the model grid xy_of_reference."""
        self._ref_coord = (new_xy_of_reference[0], new_xy_of_reference[1])

    @property
    def axis_units(self):
        """Get units for each axis.

        Returns
        -------
        tuple of str
            The units (as a string) for each of a grid's coordinates.

        Examples
        --------
        >>> from landlab import NetworkModelGrid
        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node),
        ...                         nodes_at_link)
        >>> grid.axis_units
        ('-', '-')
        >>> grid.axis_units = ('km', 'km')
        >>> grid.axis_units
        ('km', 'km')

        LLCATS: GINF
        """
        return self._axis_units

    @axis_units.setter
    def axis_units(self, new_units):
        """Set the units for each coordinate axis."""
        if len(new_units) != self.ndim:
            raise ValueError("length of units does not match grid dimension")
        self._axis_units = tuple(new_units)

    @property
    def axis_name(self):
        """Get the name of each coordinate axis.

        Returns
        -------
        tuple of str
            The names of each axis.

        Examples
        --------
        >>> from landlab import NetworkModelGrid
        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node),
        ...                         nodes_at_link)
        >>> grid.axis_name
        ('y', 'x')
        >>> grid.axis_name = ('lon', 'lat')
        >>> grid.axis_name
        ('lon', 'lat')

        LLCATS: GINF
        """
        return self._axis_name

    @axis_name.setter
    def axis_name(self, new_names):
        """Set the names of a grid's coordinate axes.

        Raises
        ------
        ValueError
            If the number of dimension do not match.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.axis_name = ('lon', 'lat')
        >>> grid.axis_name
        ('lon', 'lat')
        """
        if len(new_names) != self.ndim:
            raise ValueError("length of names does not match grid dimension")
        self._axis_name = tuple(new_names)

    @property
    @override_array_setitem_and_reset("reset_status_at_node")
    def status_at_node(self):
        """Get array of the boundary status for each node.

        Examples
        --------
        >>> from landlab import NetworkModelGrid
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
        >>> from landlab import NetworkModelGrid
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
        >>> from landlab import NetworkModelGrid
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

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def x_of_link(self):
        """Get array of the x-coordinates of link midpoints.

        Examples
        --------
        >>> from landlab import NetworkModelGrid
        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
        >>> grid.x_of_link
        array([ 0. , -0.5,  0.5])

        LLCATS: LINF MEAS
        """
        return np.mean(self.x_of_node[self.nodes_at_link], axis=1)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def y_of_link(self):
        """Get array of the y-coordinates of link midpoints.

        Examples
        --------
        >>> from landlab import NetworkModelGrid
        >>> y_of_node = (0, 1, 2, 2)
        >>> x_of_node = (0, 0, -1, 1)
        >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
        >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
        >>> grid.y_of_link
        array([ 0.5,  1.5,  1.5])

        LLCATS: LINF MEAS
        """
        return np.mean(self.y_of_node[self.nodes_at_link], axis=1)


# add only the correct functions
add_module_functions_to_class(
    NetworkModelGrid, "mappers.py", pattern="map_*", exclude="cell|patch"
)
add_module_functions_to_class(
    NetworkModelGrid, "gradients.py", pattern="calc_grad_at_link"
)
