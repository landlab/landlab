import numpy as np

from .status import StatusGrid
from .links import link_is_active, find_active_links, LinkGrid
from .links import _split_link_ends
from .cells import CellGrid
from .nodes import NodeGrid
from landlab.utils.decorators import deprecated


def _default_axis_names(n_dims):
    """Returns a tuple of the default axis names."""
    _DEFAULT_NAMES = ('z', 'y', 'x')
    return _DEFAULT_NAMES[- n_dims:]


def _default_axis_units(n_dims):
    """Returns a tuple of the default axis units."""
    return ('-', ) * n_dims


class BaseGrid(object):
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
    >>> from landlab.grid.unstructured.base import BaseGrid
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
    >>> len(ngrid.links_entering_at_node(0)) == 0
    True

    >>> tails, heads = zip(*links)
    >>> grid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]),
    ...     node_status=[0, 0, 0, 4], links=[tails, heads])
    >>> grid.status_at_node
    array([0, 0, 0, 4])
    >>> len(grid.active_links_entering_at_node(0)) == 0
    True
    >>> grid.active_links_leaving_at_node(0)
    array([0, 2])
    """

    def __init__(self, nodes, axis_name=None, axis_units=None,
                 node_status=None, links=None, cells=None):
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
        >>> from landlab.grid.unstructured.base import BaseGrid
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
        >>> len(ngrid.links_entering_at_node(0)) == 0
        True

        >>> tails, heads = zip(*links)
        >>> grid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]),
        ...     node_status=[0, 0, 0, 4], links=[tails, heads])
        >>> grid.status_at_node
        array([0, 0, 0, 4])
        >>> len(grid.active_links_entering_at_node(0)) == 0
        True
        >>> grid.active_links_leaving_at_node(0)
        array([0, 2])
        """
        self._node_grid = NodeGrid(nodes)

        self._axis_name = tuple(axis_name or _default_axis_names(self.ndim))
        self._axis_units = tuple(axis_units or _default_axis_units(self.ndim))

        if cells is not None:
            try:
                self._cell_grid = CellGrid(*cells)
            except TypeError:
                self._cell_grid = cells

        if links is not None:
            try:
                self._link_grid = LinkGrid(links, self.number_of_nodes)
            except TypeError:
                self._link_grid = links

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
        >>> from landlab.grid.unstructured.base import BaseGrid
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
        >>> from landlab.grid.unstructured.base import BaseGrid
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

    @deprecated(use='length_of_link', version=1.0)
    def link_length(self, link=None):
        return self.length_of_link(link=link)

    def length_of_link(self, link=None):
        """Length of grid links.

        Parameters
        ----------
        link : array-like, optional
            Link IDs

        Examples
        --------
        >>> from landlab.grid.unstructured.base import BaseGrid
        >>> links = [(0, 2), (1, 3), (0, 1), (2, 3), (0, 3)]
        >>> grid = BaseGrid(([0, 0, 4, 4], [0, 3, 0, 3]), links=links)
        >>> grid.length_of_link()
        array([ 4.,  4.,  3.,  3.,  5.])
        >>> grid.length_of_link(0)
        array([ 4.])

        >>> grid.length_of_link().min()
        3.0
        >>> grid.length_of_link().max()
        5.0
        """
        if link is None:
            node0, node1 = (self.node_at_link_start, self.node_at_link_end)
        else:
            node0, node1 = (self.node_at_link_start[link],
                            self.node_at_link_end[link])

        return self.node_to_node_distance(node0, node1)

    def node_to_node_distance(self, node0, node1, out=None):
        """Distance between nodes.

        Parameters
        ----------
        node0 : array-like
            Node ID of start
        node1 : array-like
            Node ID of end

        Returns
        -------
        array :
            Distances between nodes.

        Examples
        --------
        >>> from landlab.grid.unstructured.base import BaseGrid
        >>> grid = BaseGrid(([0, 0, 4, 4], [0, 3, 0, 3]))
        >>> grid.node_to_node_distance(0, 3)
        array([ 5.])
        >>> grid.node_to_node_distance(0, [0, 1, 2, 3])
        array([ 0.,  3.,  4.,  5.])
        """
        return point_to_point_distance(
            self._get_coord_at_node(node0), self._get_coord_at_node(node1),
            out=out)

        node0, node1 = np.broadcast_arrays(node0, node1)
        return np.sqrt(np.sum((self.coord_at_node[:, node1] -
                               self.coord_at_node[:, node0]) ** 2, axis=0))

    def point_to_node_distance(self, point, node=None, out=None):
        """Distance from a point to a node.

        Parameters
        ----------
        point : tuple
            Coordinates of point
        node : array-like
            Node IDs

        Returns
        -------
        array :
            Distances from point to node.

        Examples
        --------
        >>> from landlab.grid.unstructured.base import BaseGrid
        >>> grid = BaseGrid(([0, 0, 4, 4], [0, 3, 0, 3]))
        >>> grid.point_to_node_distance((0., 0.), [1, 2, 3])
        array([ 3.,  4.,  5.])
        >>> grid.point_to_node_distance((0., 0.))
        array([ 0.,  3.,  4.,  5.])
        >>> out = np.empty(4)
        >>> out is grid.point_to_node_distance((0., 0.), out=out)
        True
        >>> out
        array([ 0.,  3.,  4.,  5.])
        """
        return point_to_point_distance(point, self._get_coord_at_node(node),
                                       out=out)

    def point_to_node_angle(self, point, node=None, out=None):
        """Angle from a point to a node.

        Parameters
        ----------
        point : tuple
            Coordinates of point
        node : array-like
            Node IDs

        Returns
        -------
        array :
            Angles from point to node as radians.

        Examples
        --------
        >>> from landlab.grid.unstructured.base import BaseGrid
        >>> grid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
        >>> grid.point_to_node_angle((0., 0.), [1, 2, 3]) / np.pi
        array([ 0.  ,  0.5 ,  0.25])
        >>> grid.point_to_node_angle((0., 0.)) / np.pi
        array([ 0.  ,  0.  ,  0.5 ,  0.25])
        >>> out = np.empty(4)
        >>> out is grid.point_to_node_angle((0., 0.), out=out)
        True
        >>> out / np.pi
        array([ 0.  ,  0.  ,  0.5 ,  0.25])
        """
        return point_to_point_angle(point, self._get_coord_at_node(node),
                                    out=out)

    def point_to_node_azimuth(self, point, node=None, out=None):
        """Azimuth from a point to a node.

        Parameters
        ----------
        point : tuple
            Coordinates of point
        node : array-like
            Node IDs

        Returns
        -------
        array :
            Azimuths from point to node.

        Examples
        --------
        >>> from landlab.grid.unstructured.base import BaseGrid
        >>> grid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
        >>> grid.point_to_node_azimuth((0., 0.), [1, 2, 3])
        array([ 90.,   0.,  45.])
        >>> grid.point_to_node_azimuth((0., 0.))
        array([ 90.,  90.,   0.,  45.])
        >>> grid.point_to_node_azimuth((0., 0.), 1)
        array([ 90.])
        >>> out = np.empty(4)
        >>> out is grid.point_to_node_azimuth((0., 0.), out=out)
        True
        >>> out
        array([ 90.,  90.,   0.,  45.])
        """
        return point_to_point_azimuth(point, self._get_coord_at_node(node),
                                      out=out)

    def point_to_node_vector(self, point, node=None, out=None):
        """Azimuth from a point to a node.

        Parameters
        ----------
        point : tuple
            Coordinates of point
        node : array-like
            Node IDs

        Returns
        -------
        array :
            Vector from point to node.

        Examples
        --------
        >>> from landlab.grid.unstructured.base import BaseGrid
        >>> grid = BaseGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
        >>> grid.point_to_node_vector((0., 0.), [1, 2, 3])
        array([[ 0.,  1.,  1.],
               [ 1.,  0.,  1.]])
        >>> grid.point_to_node_vector((0., 0.))
        array([[ 0.,  0.,  1.,  1.],
               [ 0.,  1.,  0.,  1.]])
        >>> grid.point_to_node_vector((0., 0.), 1)
        array([[ 0.],
               [ 1.]])
        >>> out = np.empty((2, 1))
        >>> out is grid.point_to_node_vector((0., 0.), 1, out=out)
        True
        >>> out
        array([[ 0.],
               [ 1.]])
        """
        return point_to_point_vector(point, self._get_coord_at_node(node),
                                     out=out)

    def _get_coord_at_node(self, node=None):
        if node is None:
            return self.coord_at_node
        else:
            return self.coord_at_node[:, node].reshape((2, -1))


def point_to_point_distance(point0, point1, out=None):
    """Length of vector that joins two points.

    Parameters
    ----------
    (y0, x0) : tuple of array_like
    (y1, x1) : tuple of array_like
    out : array_like, optional
        An array to store the output. Must be the same shape as the output
        would have.

    Returns
    -------
    l : array_like
        Length of vector joining points; if *out* is provided, *v* will be
        equal to *out*.

    Examples
    --------
    >>> from landlab.grid.unstructured.base import point_to_point_distance
    >>> point_to_point_distance((0, 0), (3, 4))
    array([ 5.])
    >>> point_to_point_distance((0, 0), ([3, 6], [4, 8]))
    array([  5.,  10.])
    """
    point0 = np.reshape(point0, (2, -1))
    point1 = np.reshape(point1, (2, -1))
    if out is None:
        sum_of_squares = np.sum((point1 - point0) ** 2., axis=0)
        return np.sqrt(sum_of_squares)
    else:
        sum_of_squares = np.sum((point1 - point0) ** 2., axis=0, out=out)
        return np.sqrt(sum_of_squares, out=out)


def point_to_point_angle(point0, point1, out=None):
    """Angle of vector that joins two points.

    Parameters
    ----------
    (y0, x0) : tuple of array_like
    (y1, x1) : tuple of array_like
    out : array_like, optional
        An array to store the output. Must be the same shape as the output
        would have.

    Returns
    -------
    a : array_like
        Angle of vector joining points; if *out* is provided, *v* will be equal
        to *out*.
    """
    point0 = np.reshape(point0, (-1, 1))
    diff = point_to_point_vector(point0, point1)
    if out is None:
        return np.arctan2(diff[0], diff[1])
    else:
        return np.arctan2(diff[0], diff[1], out=out)


def point_to_point_azimuth(point0, point1, out=None):
    """Azimuth of vector that joins two points.

    Parameters
    ----------
    (y0, x0) : tuple of array_like
    (y1, x1) : tuple of array_like
    out : array_like, optional
        An array to store the output. Must be the same shape as the output
        would have.

    Returns
    -------
    azimuth : array_like
        Azimuth of vector joining points; if *out* is provided, *v* will be
        equal to *out*.

    Examples
    --------
    >>> from landlab.grid.unstructured.base import point_to_point_azimuth
    >>> point_to_point_azimuth((0, 0), (1, 0))
    array([ 0.])
    >>> point_to_point_azimuth([(0, 1), (0, 1)], (1, 0))
    array([  0., -90.])
    >>> point_to_point_azimuth([(0, 1, 0), (0, 1, 2)], [(1, 1, 2), (0, 0, 4)])
    array([  0., -90.,  45.])
    """
    azimuth_in_rads = point_to_point_angle(point0, point1, out=out)
    if out is None:
        return (np.pi * .5 - azimuth_in_rads) * 180. / np.pi
    else:
        np.subtract(np.pi * .5, azimuth_in_rads, out=out)
        return np.multiply(out, 180. / np.pi, out=out)


def point_to_point_vector(point0, point1, out=None):
    """Vector that joins two points.

    Parameters
    ----------
    (y0, x0) : tuple of array_like
    (y1, x1) : tuple of array_like
    out : array_like, optional
        An array to store the output. Must be the same shape as the output
        would have.

    Returns
    -------
    (dy, dx) : tuple of array_like
        Vectors between points; if *out* is provided, *v* will be equal to
        *out*.

    Examples
    --------
    >>> from landlab.grid.unstructured.base import point_to_point_vector
    >>> point_to_point_vector((0, 0), (1, 2))
    array([[1],
           [2]])
    >>> point_to_point_vector([(0, 1), (0, 1)], (1, 2))
    array([[1, 0],
           [2, 1]])
    >>> point_to_point_vector([(0, 0, 0), (0, 1, 2)], [(1, 2, 2), (2, 4, 4)])
    array([[1, 2, 2],
           [2, 3, 2]])
    """
    point0 = np.reshape(point0, (2, -1))
    point1 = np.reshape(point1, (2, -1))

    if out is None:
        return np.subtract(point1, point0)
    else:
        return np.subtract(point1, point0, out=out)
