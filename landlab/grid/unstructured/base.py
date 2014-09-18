import numpy as np

from .status import StatusGrid
from .links import link_is_active, find_active_links, LinkGrid, _split_link_ends
from .cells import CellGrid
from .nodes import NodeGrid


def _default_axis_names(n_dims):
    """Returns a tuple of the default axis names."""
    _DEFAULT_NAMES = ('z', 'y', 'x')
    return _DEFAULT_NAMES[- n_dims:]


def _default_axis_units(n_dims):
    """Returns a tuple of the default axis units."""
    return ('-', ) * n_dims


class BaseGrid(object):
    def __init__(self, nodes, axis_name=None, axis_units=None, node_status=None,
                 links=None, cells=None):
        """__init__([coord0, coord1, ...], axis_name=None, axis_units=None)

        Parameters
        ----------
        coord0, coord1, ... : sequence of array-like
            Coordinates of grid nodes
        axis_name : sequence of strings, optional
            Names of coordinate axes
        axis_units : sequence of strings, optional
            Units of coordinate axes

        Returns
        -------
        BaseGrid :
            A newly-created BaseGrid

        Examples
        --------
        >>> ngrid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
        >>> ngrid.number_of_nodes
        4
        >>> ngrid.x_at_node
        array([ 0.,  1.,  0.,  1.])
        >>> ngrid.x_at_node[2]
        0.0
        >>> ngrid.point_at_node[2]
        array([ 1.,  0.])
        >>> ngrid.coord_at_node[:, [2, 3]]
        array([[ 1.,  1.],
               [ 0.,  1.]])
        >>> ngrid.distance_between_nodes([0, 1, 2, 3], 0)
        array([ 0.        ,  1.        ,  1.        ,  1.41421356])

        >>> cells = ([0, 1, 2, 1, 3, 2], [3, 3], [0, 1])
        >>> ngrid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]), cells=cells)
        >>> ngrid.number_of_cells
        2
        >>> ngrid.node_at_cell
        array([0, 1])

        >>> links = [(0, 2), (1, 3), (0, 1), (1, 2), (0, 3)]
        >>> ngrid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]), links=zip(*links))
        >>> ngrid.number_of_links
        5
        >>> ngrid.links_leaving_at_node(0)
        array([0, 2, 4])
        >>> ngrid.links_entering_at_node(0)
        array([], dtype=int64)

        >>> grid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]),
        ...     node_status=[0, 0, 0, 4], links=zip(*links))
        >>> grid.status_at_node
        array([0, 0, 0, 4])
        >>> grid.active_links_entering_at_node(0)
        array([], dtype=int64)
        >>> grid.active_links_leaving_at_node(0)
        array([0, 2])
        """
        self._node_grid = NodeGrid(nodes)

        self._axis_name = tuple(axis_name or _default_axis_names(self.ndim))
        self._axis_units = tuple(axis_units or _default_axis_units(self.ndim))

        if cells is not None:
            self._cell_grid = CellGrid(*cells)

        if links is not None:
            self._link_grid = LinkGrid(links, self.number_of_nodes)

        if node_status is not None:
            self._status_grid = StatusGrid(node_status)

        if links is not None and node_status is not None:
            links = _split_link_ends(links)
            self._active_link_grid = BaseGrid.create_active_link_grid(
                self.status_at_node, links, self.number_of_nodes)

    @staticmethod
    def create_active_link_grid(node_status, links, number_of_nodes):
        active_link_ids = find_active_links(node_status, links)
        return LinkGrid((links[0][active_link_ids], links[1][active_link_ids]),
                        number_of_nodes, link_ids=active_link_ids)

    @property
    def ndim(self):
        return self._node_grid.ndim

    @property
    def axis_units(self):
        """Coordinate units of each axis.

        Returns
        -------
        tuple of strings :
            Coordinate units of each axis.

        Examples
        --------
        >>> ngrid = BaseGrid(([0, 1, 0], [1, 1, 0]))
        >>> ngrid.axis_units
        ('-', '-')

        >>> ngrid = BaseGrid(([0, 1, 0], [1, 1, 0]),
        ...     axis_units=['degrees_north', 'degrees_east'])
        >>> ngrid.axis_units
        ('degrees_north', 'degrees_east')
        """
        return self._axis_units

    @property
    def axis_name(self):
        """Name of each axis.

        Returns
        -------
        tuple of strings :
            Names of each axis.

        Examples
        --------
        >>> ngrid = BaseGrid(([0, 1, 0], [1, 1, 0]))
        >>> ngrid.axis_name
        ('y', 'x')

        >>> ngrid = BaseGrid(([0, 1, 0], [1, 1, 0]), axis_name=['lat', 'lon'])
        >>> ngrid.axis_name
        ('lat', 'lon')
        """
        return self._axis_name

    @property
    def number_of_links(self):
        """Number of links.
        """
        return self._link_grid.number_of_links

    @property
    def number_of_cells(self):
        """Number of cells.
        """
        return self._cell_grid.number_of_cells

    @property
    def number_of_nodes(self):
        """Number of nodes.
        """
        return self._node_grid.number_of_nodes

    @property
    def coord_at_node(self):
        return self._node_grid.coord

    @property
    def x_at_node(self):
        return self._node_grid.x

    @property
    def y_at_node(self):
        return self._node_grid.y

    @property
    def point_at_node(self):
        return self._node_grid.point

    def links_leaving_at_node(self, node):
        return self._link_grid.out_link_at_node(node)

    def links_entering_at_node(self, node):
        return self._link_grid.in_link_at_node(node)

    def active_links_leaving_at_node(self, node):
        return self._active_link_grid.out_link_at_node(node)

    def active_links_entering_at_node(self, node):
        return self._active_link_grid.in_link_at_node(node)

    @property
    def node_at_link_start(self):
        return self._link_grid.node_at_link_start

    @property
    def node_at_link_end(self):
        return self._link_grid.node_at_link_end

    @property
    def node_at_cell(self):
        return self._cell_grid.node_at_cell

    @property
    def cell_at_node(self):
        return self._cell_grid.cell_at_node

    def core_cells(self):
        return self.cell_at_node[self.core_nodes]

    @property
    def status_at_node(self):
        return self._status_grid.node_status

    @status_at_node.setter
    def status_at_node(self, status):
        self._status_grid.node_status = status
        self._active_link_grid = BaseGrid.create_active_link_grid(
            self.status_at_node, (self.node_at_link_start,
                                  self.node_at_link_end), self.number_of_nodes)

    def active_nodes(self):
        return self._status_grid.active_nodes()

    def core_nodes(self):
        return self._status_grid.core_nodes()

    def boundary_nodes(self):
        return self._status_grid.boundary_nodes()

    def closed_boundary_nodes(self):
        return self._status_grid.closed_boundary_nodes()

    def fixed_gradient_boundary_nodes(self):
        return self._status_grid.fixed_gradient_boundary_nodes()

    def fixed_value_boundary_nodes(self):
        return self._status_grid.fixed_value_boundary_nodes()

    def active_links(self):
        return self._active_link_grid.link_id

    def link_length(self, link):
        return self.distance_between_nodes(self.node_at_link_start[link],
                                           self.node_at_link_end[link])

    def distance_between_nodes(self, node0, node1):
        node0, node1 = np.broadcast_arrays(node0, node1)
        return np.sqrt(np.sum((self.coord_at_node[:, node1] -
                               self.coord_at_node[:, node0]) ** 2, axis=0))
