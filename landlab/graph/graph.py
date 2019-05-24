"""Define a graph of nodes-links-patches.

Nodes and links are required. If no patches are provided, no patches will
be created.

Examples
--------

>>> from landlab.graph import NetworkGraph, Graph

>>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
>>> graph = NetworkGraph((node_y, node_x))
>>> graph.x_of_node
array([ 0.,  1.,  2.,  0.,  1.,  2.,  0.,  1.,  2.])
>>> graph.y_of_node
array([ 0.,  0.,  0.,  1.,  1.,  1.,  2.,  2.,  2.])
>>> graph.ndim
2

>>> links = ((0, 1), (1, 2),
...          (0, 3), (1, 4), (2, 5),
...          (3, 4), (4, 5),
...          (3, 6), (4, 7), (5, 8),
...          (6, 7), (7, 8))
>>> graph = Graph((node_y, node_x), links=links)
>>> graph.nodes_at_link # doctest: +NORMALIZE_WHITESPACE
array([[0, 1], [1, 2],
       [0, 3], [1, 4], [2, 5],
       [3, 4], [4, 5],
       [3, 6], [4, 7], [5, 8],
       [6, 7], [7, 8]])
>>> graph.node_at_link_head
array([1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8])
>>> graph.node_at_link_tail
array([0, 1, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7])

>>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
       [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
       [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])

>>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
       [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
       [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]])

>>> patches = ((5, 3, 0, 2), (6, 4, 1, 3), (10, 8, 5, 7), (11, 9, 6, 8))
>>> graph = Graph((node_y, node_x), links=links, patches=patches)
>>> graph.links_at_patch
array([[ 3,  5,  2,  0],
       [ 4,  6,  3,  1],
       [ 8, 10,  7,  5],
       [ 9, 11,  8,  6]])
>>> graph.nodes_at_patch
array([[4, 3, 0, 1],
       [5, 4, 1, 2],
       [7, 6, 3, 4],
       [8, 7, 4, 5]])
"""
import json

import numpy as np
import six
import xarray as xr

from ..core.utils import as_id_array
from ..utils.decorators import read_only_array, store_result_in_grid
from .object.at_node import get_links_at_node
from .object.at_patch import get_nodes_at_patch
from .quantity.of_link import (
    get_angle_of_link,
    get_length_of_link,
    get_midpoint_of_link,
)
from .quantity.of_patch import get_area_of_patch, get_centroid_of_patch
from .sort import reindex_by_xy, reorder_links_at_patch
from .sort.sort import reorient_link_dirs, reverse_one_to_many
from .ugrid import ugrid_from_unstructured


def _parse_sorting_opt(sorting):
    SORTING_OPTS = ("xy", "ccw", "ne")

    as_dict = None

    if isinstance(sorting, bool):
        as_dict = dict([(opt, format) for opt in SORTING_OPTS])
    elif isinstance(sorting, dict):
        as_dict = dict(sorting.items())
        for opt in SORTING_OPTS:
            sorting.setdefault(opt, True)

    return as_dict


def find_perimeter_nodes(graph):
    """Find nodes on the perimeter of a graph.

    Uses a convex hull to locate the perimeter nodes of a graph.

    Parameters
    ----------
    graph : graph_like
        A Graph of nodes (just requires *xy_of_node*).

    Returns
    -------
    ndarray of int
        Identifiers of the perimeter nodes.
    """
    from scipy.spatial import ConvexHull

    hull = ConvexHull(graph.xy_of_node, qhull_options="Qt")
    return as_id_array(hull.vertices)


class thawed(object):
    def __init__(self, graph):
        self._graph = graph
        self._initially_frozen = graph.frozen

    def __enter__(self):
        self._graph.thaw()

    def __exit__(self, ex_type, ex_value, traceback):
        if self._initially_frozen:
            self._graph.freeze()


class NetworkGraph(object):
    """Define the connectivity of a graph of nodes and links.

    Unlike Graph, NetworkGraph does not have patches.
    """

    def __init__(self, mesh, **kwds):
        """Define a graph of connected nodes.


        Parameters
        ----------
        mesh : Dataset
            xarray Dataset that defines the topology in ugrid format.
        """
        if not isinstance(mesh, xr.Dataset):
            node_y_and_x = mesh
            links = kwds.get("links", None)
            mesh = ugrid_from_unstructured(node_y_and_x, links=links)

        self._ds = mesh

        self._frozen = False
        self.freeze()

        if kwds.get("sort", True):
            NetworkGraph.sort(self)

        self._origin = (0.0, 0.0)

    @property
    def frozen(self):
        return self._frozen

    def thawed(self):
        return thawed(self)

    def sort(self):
        with self.thawed():
            reorient_link_dirs(self)
            sorted_nodes, sorted_links, sorted_patches = reindex_by_xy(self)

        return sorted_nodes, sorted_links, sorted_patches

    def freeze(self):
        """Freeze the graph by making arrays read-only."""
        for var in self.ds.variables:
            self.ds[var].values.flags.writeable = False
        self._frozen = True

    def thaw(self):
        """Thaw the graph by making arrays writable."""
        for var in self.ds.variables:
            self.ds[var].values.flags.writeable = True
        self._frozen = False

    def _add_variable(self, name, var, dims=None, attrs=None):
        kwds = dict(data=var, dims=dims, attrs=attrs)
        self.ds.update({name: xr.DataArray(**kwds)})
        if self._frozen:
            self.freeze()

    @property
    def ds(self):
        return self._ds

    def to_dict(self):
        return self.ds.to_dict()

    def to_json(self):
        return json.dumps(self.ds.to_dict())

    def to_netcdf(self, *args, **kwds):
        """Write graph contents to a netCDF file.

        See xarray.Dataset.to_netcdf for a complete list of parameters.
        Below are only the most common.

        Parameters
        ----------
        path : str, optional
            Path to which to save this graph.
        mode : {'w', 'a'}, optional
            Write ('w') or append ('a') mode. If mode='w', any
            existing file at this location will be overwritten.
        format : {'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'}, optional
            File format for the resulting netCDF file:

            * NETCDF4: Data is stored in an HDF5 file, using netCDF4 API
              features.
            * NETCDF4_CLASSIC: Data is stored in an HDF5 file, using only
              netCDF 3 compatible API features.
            * NETCDF3_64BIT: 64-bit offset version of the netCDF 3 file format,
              which fully supports 2+ GB files, but is only compatible with
              clients linked against netCDF version 3.6.0 or later.
            * NETCDF3_CLASSIC: The classic netCDF 3 file format. It does not
              handle 2+ GB files very well.

            All formats are supported by the netCDF4-python library.
            scipy.io.netcdf only supports the last two formats.

            The default format is NETCDF4 if you are saving a file to disk and
            have the netCDF4-python library available. Otherwise, xarray falls
            back to using scipy to write netCDF files and defaults to the
            NETCDF3_64BIT format (scipy does not support netCDF4).
        """
        self.ds.to_netcdf(*args, **kwds)

    @classmethod
    def from_netcdf(cls, fname):
        return cls.from_dataset(xr.open_dataset(fname))

    @classmethod
    def from_dict(cls, meta):
        return cls(
            (meta["y_of_node"], meta["x_of_node"]),
            links=meta.get("nodes_at_link", None),
            patches=meta.get("links_at_patch", None),
        )

    @classmethod
    def load(cls, source):
        if isinstance(source, six.string_types):
            return cls.from_netcdf(source)
        elif isinstance(source, (dict, xr.Dataset)):
            return cls.from_dict(source)
        else:
            raise ValueError(
                "source must be dict-like or NetCDF ({type})".format(type=type(source))
            )

    def __str__(self):
        return str(self.ds)

    def __repr__(self):
        return repr(self.ds)

    @property
    def ndim(self):
        return 2

    @property
    @store_result_in_grid()
    @read_only_array
    def xy_of_node(self):
        """Get x and y-coordinates of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.xy_of_node[:, 0]
        array([ 0.,  1.,  2.,  0.,  1.,  2.])
        >>> graph.xy_of_node[:, 1]
        array([ 0.,  0.,  0.,  1.,  1.,  1.])

        LLCATS: NINF
        """
        return np.stack((self.x_of_node, self.y_of_node)).T.copy()

    @property
    def x_of_node(self):
        """Get x-coordinate of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.x_of_node
        array([ 0.,  1.,  2.,  0.,  1.,  2.])

        LLCATS: NINF
        """
        return self.ds["x_of_node"].values

    @property
    def y_of_node(self):
        """Get y-coordinate of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.y_of_node
        array([ 0.,  0.,  0.,  1.,  1.,  1.])

        LLCATS: NINF
        """
        return self.ds["y_of_node"].values

    @property
    def nodes(self):
        """Get identifier for each node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.nodes
        array([0, 1, 2, 3, 4, 5])

        LLCATS: NINF
        """
        return self.ds["node"].values

    @property
    @store_result_in_grid()
    @read_only_array
    def perimeter_nodes(self):
        """Get nodes on the convex hull of a Graph.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> np.sort(graph.perimeter_nodes)
        array([0, 2, 3, 5])

        LLCATS: NINF BC SUBSET
        """
        return find_perimeter_nodes(self)

    @property
    def number_of_nodes(self):
        """Get total number of nodes.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.number_of_nodes == 6
        True

        LLCATS: NINF
        """
        return self.ds.dims["node"]

    @property
    def nodes_at_link(self):
        """Get nodes at either end of links.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.nodes_at_link # doctest: +NORMALIZE_WHITESPACE
        array([[0, 1], [1, 2],
               [0, 3], [1, 4], [2, 5],
               [3, 4], [4, 5],
               [3, 6], [4, 7], [5, 8],
               [6, 7], [7, 8]])

        LLCATS: NINF
        """
        return self.ds["nodes_at_link"].values

    @property
    def node_at_link_tail(self):
        """Get nodes at link tail.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.node_at_link_tail
        array([0, 1, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7])

        LLCATS: NINF
        """
        return self.nodes_at_link[:, 0]

    @property
    def node_at_link_head(self):
        """Get nodes at link head.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.node_at_link_head
        array([1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8])

        LLCATS: NINF
        """
        return self.nodes_at_link[:, 1]

    @property
    def number_of_links(self):
        """Get nodes at link head.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.number_of_links == 12
        True

        LLCATS: LINF
        """
        try:
            return self.ds.dims["link"]
        except KeyError:
            return 0

    @property
    def links_at_node(self):
        """Get links touching a node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        >>> node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
               [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
               [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])

        LLCATS: LINF
        """
        try:
            return self._links_at_node
        except AttributeError:
            (
                self._links_at_node,
                self._link_dirs_at_node,
            ) = self._create_links_and_dirs_at_node()
            return self._links_at_node

    def _create_links_and_dirs_at_node(self):
        return get_links_at_node(self, sort=True)
        # return get_links_at_node(self, sort=self._sorting['ccw'])

    @property
    def link_dirs_at_node(self):
        """Get directions of links touching a node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
               [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
               [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]])

        LLCATS: LINF
        """
        try:
            return self._link_dirs_at_node
        except AttributeError:
            (
                self._links_at_node,
                self._link_dirs_at_node,
            ) = self._create_links_and_dirs_at_node()
            return self._link_dirs_at_node

    @property
    @store_result_in_grid()
    @read_only_array
    def angle_of_link(self):
        """Get the angle of each link.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import Graph

        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.angle_of_link * 180. / np.pi
        array([  0.,   0.,  90.,  90.,  90.,   0.,   0.])

        LLCATS: LINF
        """
        return get_angle_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def length_of_link(self):
        """Get the length of links.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import UniformRectilinearGraph

        >>> graph = UniformRectilinearGraph((2, 3), spacing=(1, 2))
        >>> graph.length_of_link
        array([ 2.,  2.,  1.,  1.,  1.,  2.,  2.])

        LLCATS: LINF
        """
        return get_length_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def midpoint_of_link(self):
        """Get the middle of links.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import UniformRectilinearGraph

        >>> graph = UniformRectilinearGraph((2, 3), spacing=(1, 2))
        >>> graph.midpoint_of_link # doctest: +NORMALIZE_WHITESPACE
        array([[ 1. ,  0. ], [ 3. ,  0. ],
               [ 0. ,  0.5], [ 2. ,  0.5], [ 4. ,  0.5],
               [ 1. ,  1. ], [ 3. ,  1. ]])

        LLCATS: LINF
        """
        return get_midpoint_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def xy_of_link(self):
        return get_midpoint_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def adjacent_nodes_at_node(self):
        """Get adjacent nodes.

        Examples
        --------
        >>> from landlab.graph import Graph

        First, a simple example with no diagonals.

        >>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.adjacent_nodes_at_node
        array([[ 1,  3, -1, -1],
               [ 2,  4,  0, -1],
               [ 5,  1, -1, -1],
               [ 4,  6,  0, -1],
               [ 5,  7,  3,  1],
               [ 8,  4,  2, -1],
               [ 7,  3, -1, -1],
               [ 8,  6,  4, -1],
               [ 7,  5, -1, -1]])

        Next, we add the diagonal from node 0 to node 4.

        >>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8),
        ...          (0, 4))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.adjacent_nodes_at_node
        array([[ 1,  4,  3, -1, -1],
               [ 2,  4,  0, -1, -1],
               [ 5,  1, -1, -1, -1],
               [ 4,  6,  0, -1, -1],
               [ 5,  7,  3,  0,  1],
               [ 8,  4,  2, -1, -1],
               [ 7,  3, -1, -1, -1],
               [ 8,  6,  4, -1, -1],
               [ 7,  5, -1, -1, -1]])

        LLCATS: NINF
        """
        node_is_at_tail = np.choose(
            self.link_dirs_at_node + 1, np.array((1, -1, 0), dtype=np.int8)
        )
        out = self.nodes_at_link[self.links_at_node, node_is_at_tail]
        out[node_is_at_tail == -1] = -1

        return out


class Graph(NetworkGraph):

    """Define the connectivity of a graph of nodes, links, and patches."""

    def __init__(self, mesh, **kwds):
        """Define a graph of connected nodes.

        Parameters
        ----------
        mesh : Dataset
            xarray Dataset that defines the topology in ugrid format.
        """
        if not isinstance(mesh, xr.Dataset):
            node_y_and_x = mesh
            links = kwds.get("links", None)
            patches = kwds.get("patches", None)
            mesh = ugrid_from_unstructured(node_y_and_x, links=links, patches=patches)

        self._ds = mesh

        self._frozen = False
        self.freeze()

        if kwds.get("sort", True):
            Graph.sort(self)

        self._origin = (0.0, 0.0)

    def sort(self):
        with self.thawed():
            reorient_link_dirs(self)
            sorted_nodes, sorted_links, sorted_patches = reindex_by_xy(self)
            reorder_links_at_patch(self)

        return sorted_nodes, sorted_links, sorted_patches

    @property
    @store_result_in_grid()
    @read_only_array
    def xy_of_patch(self):
        """Get the centroid of each patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.xy_of_patch
        array([[ 0.5,  0.5],
              [ 1.5,  0.5]])

        LLCATS: PINF
        """
        return get_centroid_of_patch(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def area_of_patch(self):
        """Get the area of each patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.area_of_patch
        array([ 1.,  1.])

        LLCATS: PINF
        """
        return get_area_of_patch(self)

    @property
    def number_of_patches(self):
        """Get the number of patches.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.number_of_patches == 2
        True

        LLCATS: PINF
        """
        try:
            return self.ds.dims["patch"]
        except KeyError:
            return 0

    @property
    def links_at_patch(self):
        """Get the links that define a patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.links_at_patch
        array([[3, 5, 2, 0],
               [4, 6, 3, 1]])

        LLCATS: LINF
        """
        return self.ds["links_at_patch"].values

    @property
    # @store_result_in_grid()
    @read_only_array
    def nodes_at_patch(self):
        """Get the nodes that define a patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1, 2, 2, 2])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.nodes_at_patch
        array([[4, 3, 0, 1],
               [5, 4, 1, 2]])

        LLCATS: NINF
        """
        return get_nodes_at_patch(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def patches_at_node(self):
        """Get the patches that touch each node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.patches_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0, -1], [ 0,  1], [ 1, -1],
               [ 0, -1], [ 0,  1], [ 1, -1]])

        LLCATS: PINF
        """
        return reverse_one_to_many(self.nodes_at_patch)

    @property
    @store_result_in_grid()
    @read_only_array
    def patches_at_link(self):
        """Get the patches on either side of each link.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.patches_at_link # doctest: +NORMALIZE_WHITESPACE
        array([[ 0, -1], [ 1, -1],
               [ 0, -1], [ 0,  1], [ 1, -1],
               [ 0, -1], [ 1, -1]])

        LLCATS: PINF
        """
        return reverse_one_to_many(self.links_at_patch, min_counts=2)
        try:
            return self.ds["patches_at_link"].values
        except KeyError:
            patches_at_link = xr.DataArray(
                data=reverse_one_to_many(self.links_at_patch, min_counts=2),
                dims=("link", "Two"),
                attrs={
                    "cf_role": "edge_node_connectivity",
                    "long_name": "patches on either side of a link",
                    "start_index": 0,
                },
            )
            self.ds.update({"patches_at_link": patches_at_link})
            return self.ds["patches_at_link"].values
