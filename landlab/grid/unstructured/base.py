import numpy as np

from .status import StatusGrid
from .links import link_is_active, LinkGrid
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

        >>> cells = ([0, 1, 2, 1, 3, 2], [3, 3])
        >>> ngrid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]), cells=cells)
        >>> ngrid.number_of_cells
        2
        >>> ngrid.nodes_at_cell(0)
        array([0, 1, 2])

        >>> links = [(0, 2), (1, 3), (0, 1), (1, 2), (0, 3)]
        >>> ngrid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]), links=zip(*links))
        >>> ngrid.number_of_links
        5
        >>> ngrid.links_leaving_at_node(0)
        array([0, 2, 4])
        >>> ngrid.links_entering_at_node(0)
        array([], dtype=int64)
        """
        self._axis_name = tuple(axis_name or _default_axis_names(self.ndim))
        self._axis_units = tuple(axis_units or _default_axis_units(self.ndim))

        self._node_grid = NodeGrid(nodes)

        if cells is not None:
            self._cell_grid = CellGrid(*cells)
        if links is not None:
            self._link_grid = LinkGrid(links, self.number_of_nodes)
        if node_status is not None:
            self._status_grid = StatusGrid(node_status)

    @property
    def ndim(self):
        return 2

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

    def nodes_at_cell(self, cell):
        return self._cell_grid.nodes_at_cell(cell)

    def links_leaving_at_node(self, node):
        return self._link_grid.out_link_at_node(node)

    def links_entering_at_node(self, node):
        return self._link_grid.in_link_at_node(node)

    @property
    def node_id_at_cells(self):
        return self._node_id_at_cells

    def _reset_list_of_active_links(self):
        """
        Creates or resets a list of active links. We do this by sweeping
        through the given lists of from and to nodes, and checking the status
        of these as given in the node_status list. A link is active if both its
        nodes are active interior points, or if one is an active interior and
        the other is an active boundary.
        """
        status_at_link_ends = (self.node_status[self._node_id_at_link_start],
                               self.node_status[self._node_id_at_link_end])

        (self._active_link_ids, ) = np.where(
            link_is_active(status_at_link_ends))

        self._number_of_active_links = len(self._active_link_ids)
        self._number_of_active_faces = len(self._active_link_ids)
        self._node_id_at_active_link_start = self._node_id_at_link_start[self._active_link_ids]
        self._node_id_at_active_link_end = self._node_id_at_link_end[self._active_link_ids]

        # Set up active inlink and outlink matrices
        self._setup_active_inlink_and_outlink_matrices()
