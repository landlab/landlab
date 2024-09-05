"""Define a graph of nodes-links-patches.

Nodes and links are required. If no patches are provided, no patches will
be created.

Examples
--------

>>> from landlab.graph import NetworkGraph, Graph

>>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
>>> graph = NetworkGraph((node_y, node_x), sort=True)
>>> graph.x_of_node
array([0., 1., 2., 0., 1., 2., 0., 1., 2.])
>>> graph.y_of_node
array([0., 0., 0., 1., 1., 1., 2., 2., 2.])
>>> graph.ndim
2

>>> links = [
...     (0, 1),
...     (1, 2),
...     (0, 3),
...     (1, 4),
...     (2, 5),
...     (3, 4),
...     (4, 5),
...     (3, 6),
...     (4, 7),
...     (5, 8),
...     (6, 7),
...     (7, 8),
... ]
>>> graph = Graph((node_y, node_x), links=links, sort=True)
>>> graph.nodes_at_link
array([[0, 1], [1, 2],
       [0, 3], [1, 4], [2, 5],
       [3, 4], [4, 5],
       [3, 6], [4, 7], [5, 8],
       [6, 7], [7, 8]])
>>> graph.node_at_link_head
array([1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8])
>>> graph.node_at_link_tail
array([0, 1, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7])

>>> graph.links_at_node
array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
       [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
       [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])

>>> graph.link_dirs_at_node
array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
       [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
       [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]],
      dtype=int8)

>>> patches = ((5, 3, 0, 2), (6, 4, 1, 3), (10, 8, 5, 7), (11, 9, 6, 8))
>>> graph = Graph((node_y, node_x), links=links, patches=patches, sort=True)
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
from functools import cached_property

import numpy as np
import xarray as xr

from ..core.utils import as_id_array
from ..utils.decorators import read_only_array
from .object.at_node import get_links_at_node
from .object.at_patch import get_nodes_at_patch
from .quantity.of_link import get_angle_of_link
from .quantity.of_link import get_length_of_link
from .quantity.of_link import get_midpoint_of_link
from .quantity.of_patch import get_area_of_patch
from .quantity.of_patch import get_centroid_of_patch
from .sort import reindex_by_xy
from .sort.sort import reorient_link_dirs
from .sort.sort import reverse_one_to_many
from .sort.sort import sort_spokes_at_hub
from .ugrid import ugrid_from_unstructured


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


class thawed:
    def __init__(self, graph):
        self._graph = graph
        self._initially_frozen = graph.frozen

    def __enter__(self):
        self._graph.thaw()

    def __exit__(self, ex_type, ex_value, traceback):
        if self._initially_frozen:
            self._graph.freeze()


def _update_node_at_cell(ugrid, node_at_cell):
    node_at_cell = xr.DataArray(
        data=as_id_array(node_at_cell),
        dims=("cell",),
        attrs={
            "cf_role": "cell_node_connectivity",
            "long_name": "nodes centered at cells",
            "start_index": 0,
        },
    )
    ugrid.update({"node_at_cell": node_at_cell})


def _update_nodes_at_face(ugrid, nodes_at_face):
    nodes_at_face = xr.DataArray(
        data=as_id_array(nodes_at_face),
        dims=("face", "Two"),
        attrs={
            "cf_role": "face_node_connectivity",
            "long_name": "nodes on either side of a face",
            "start_index": 0,
        },
    )
    ugrid.update({"nodes_at_face": nodes_at_face})


class NetworkGraph:
    """Define the connectivity of a graph of nodes and links.

    Unlike Graph, NetworkGraph does not have patches.
    """

    def __init__(self, node_y_and_x, links=None, sort=False):
        """Define a graph of connected nodes.

        Parameters
        ----------
        mesh : Dataset
            xarray Dataset that defines the topology in ugrid format.
        """
        self._ds = ugrid_from_unstructured(node_y_and_x, links=links)

        self._frozen = False
        self.freeze()

        if sort:
            (
                self._sorted_nodes,
                self._sorted_links,
                self._sorted_patches,
            ) = NetworkGraph.sort(self)

        self._origin = (0.0, 0.0)

    @property
    def frozen(self):
        return self._frozen

    def thawed(self):
        return thawed(self)

    def sort(self):
        """Sort graph elements."""
        with self.thawed():
            reorient_link_dirs(self)
            sorted_nodes, sorted_links, sorted_patches = reindex_by_xy(self)
            if "links_at_patch" in self.ds:
                sort_spokes_at_hub(
                    self.links_at_patch,
                    np.round(self.xy_of_patch, decimals=4),
                    np.round(self.xy_of_link, decimals=4),
                    inplace=True,
                )

        return sorted_nodes, sorted_links, sorted_patches

    def freeze(self):
        """Freeze the graph by making arrays read-only."""
        for var in self.ds.variables:
            array = self.ds[var].values
            while array is not None:
                array.flags.writeable = False
                array = array.base
        self._frozen = True

    def thaw(self):
        """Thaw the graph by making arrays writable."""
        for var in self.ds.variables:
            arrays = []
            array = self.ds[var].values
            while array is not None:
                arrays.append(array)
                array = array.base
            for array in arrays[::-1]:
                array.flags.writeable = True
        self._frozen = False

    def _add_variable(self, name, var, dims=None, attrs=None):
        kwds = {"data": var, "dims": dims, "attrs": attrs}
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
        if isinstance(source, str):
            return cls.from_netcdf(source)
        elif isinstance(source, (dict, xr.Dataset)):
            return cls.from_dict(source)
        else:
            raise ValueError(f"source must be dict-like or NetCDF ({type(source)})")

    def __str__(self):
        return str(self.ds)

    def __repr__(self):
        return repr(self.ds)

    @property
    def ndim(self):
        return 2

    @cached_property
    @read_only_array
    def xy_of_node(self):
        """Get x and y-coordinates of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.xy_of_node[:, 0]
        array([0., 1., 2., 0., 1., 2.])
        >>> graph.xy_of_node[:, 1]
        array([0., 0., 0., 1., 1., 1.])

        :meta landlab: info-node
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
        array([0., 1., 2., 0., 1., 2.])

        :meta landlab: info-node
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
        array([0., 0., 0., 1., 1., 1.])

        :meta landlab: info-node
        """
        return self.ds["y_of_node"].values

    @cached_property
    def node_x(self):
        return self.x_of_node

    @cached_property
    def node_y(self):
        return self.y_of_node

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

        :meta landlab: info-node
        """
        return self.ds["node"].values

    @cached_property
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

        :meta landlab: info-node, boundary-condition, subset
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

        :meta landlab: info-node
        """
        return self.ds.sizes["node"]

    @property
    def nodes_at_link(self):
        """Get nodes at either end of links.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.nodes_at_link
        array([[0, 1], [1, 2],
               [0, 3], [1, 4], [2, 5],
               [3, 4], [4, 5],
               [3, 6], [4, 7], [5, 8],
               [6, 7], [7, 8]])

        :meta landlab: info-node
        """
        return self.ds["nodes_at_link"].values

    @property
    def node_at_link_tail(self):
        """Get nodes at link tail.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.node_at_link_tail
        array([0, 1, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7])

        :meta landlab: info-node
        """
        return self.nodes_at_link[:, 0]

    @property
    def node_at_link_head(self):
        """Get nodes at link head.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.node_at_link_head
        array([1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8])

        :meta landlab: info-node
        """
        return self.nodes_at_link[:, 1]

    @property
    def number_of_links(self):
        """Get nodes at link head.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.number_of_links == 12
        True
        """
        try:
            return self.ds.sizes["link"]
        except KeyError:
            return 0

    @cached_property
    @read_only_array
    def links_at_node(self):
        """Get links touching a node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        >>> node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.links_at_node
        array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
               [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
               [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])

        :meta landlab: info-link
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

    @cached_property
    @read_only_array
    def link_dirs_at_node(self):
        """Return link directions into each node.

        A value of 1 indicates a link points toward a given node, while a value
        of -1 indicates a link points away from a node.

        Returns
        -------
        (n_nodes, max_links_per_node) ndarray of int
            Link directions relative to the nodes of a grid. The shape of the
            matrix will be number of nodes by the maximum number of links per
            node. A zero indicates no link at this position.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.link_dirs_at_node
        array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
               [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
               [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]],
              dtype=int8)
        """
        try:
            return self._link_dirs_at_node
        except AttributeError:
            (
                self._links_at_node,
                self._link_dirs_at_node,
            ) = self._create_links_and_dirs_at_node()
            return self._link_dirs_at_node

    @cached_property
    @read_only_array
    def angle_of_link(self):
        """Get the angle of each link.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import Graph

        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.angle_of_link * 180.0 / np.pi
        array([ 0.,  0., 90., 90., 90.,  0.,  0.])

        :meta landlab: info-link
        """
        return get_angle_of_link(self)

    @cached_property
    @read_only_array
    def length_of_link(self):
        """Get the length of links.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import UniformRectilinearGraph

        >>> graph = UniformRectilinearGraph((2, 3), spacing=(1, 2))
        >>> graph.length_of_link
        array([2., 2., 1., 1., 1., 2., 2.])

        :meta landlab: info-link
        """
        return get_length_of_link(self)

    @cached_property
    @read_only_array
    def midpoint_of_link(self):
        """Get the middle of links.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import UniformRectilinearGraph

        >>> graph = UniformRectilinearGraph((2, 3), spacing=(1, 2))
        >>> graph.midpoint_of_link
        array([[1. , 0. ], [3. , 0. ],
               [0. , 0.5], [2. , 0.5], [4. , 0.5],
               [1. , 1. ], [3. , 1. ]])

        :meta landlab: info-link
        """
        return get_midpoint_of_link(self)

    @cached_property
    @read_only_array
    def xy_of_link(self):
        return get_midpoint_of_link(self)

    @cached_property
    @read_only_array
    def adjacent_nodes_at_node(self):
        """Get adjacent nodes.

        Examples
        --------
        >>> from landlab.graph import Graph

        First, a simple example with no diagonals.

        >>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> graph = Graph((node_y, node_x), links=links, sort=True)
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
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ...     (0, 4),
        ... )
        >>> graph = Graph((node_y, node_x), links=links, sort=True)
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

        :meta landlab: info-node
        """
        node_is_at_tail = np.choose(
            self.link_dirs_at_node + 1, np.array((1, -1, 0), dtype=np.int8)
        )
        out = self.nodes_at_link[self.links_at_node, node_is_at_tail]
        out[node_is_at_tail == -1] = -1

        return out

    @cached_property
    @read_only_array
    def adjacent_links_at_link(self):
        from .object.ext.at_link import find_adjacent_links_at_link

        adjacent_links_at_link = np.empty((self.number_of_links, 2), dtype=int)

        find_adjacent_links_at_link(
            self.nodes_at_link, self.links_at_node, adjacent_links_at_link
        )

        return adjacent_links_at_link

    @cached_property
    @read_only_array
    def unit_vector_at_link(self):
        """Make arrays to store the unit vectors associated with each link.

        For each link, the x and y components of the link's unit vector (that
        is, the link's x and y dimensions if it were shrunk to unit length but
        retained its orientation).

        Examples
        --------

        The example below is a seven-node hexagonal grid, with six nodes around
        the perimeter and one node (#3) in the interior. There are four
        horizontal links with unit vector (1,0), and 8 diagonal links with
        unit vector (+/-0.5, +/-sqrt(3)/2) (note: sqrt(3)/2 ~ 0.866).

        >>> from landlab.graph import TriGraph
        >>> graph = TriGraph((3, 2), spacing=2.0, node_layout="hex", sort=True)

        >>> np.round(graph.unit_vector_at_link[:, 0], decimals=5)
        array([ 1. , -0.5,  0.5, -0.5,  0.5,  1. ,  1. ,  0.5, -0.5,  0.5, -0.5,
                1. ])
        >>> np.round(graph.unit_vector_at_link[:, 1], decimals=5)
        array([0.     , 0.86603, 0.86603, 0.86603, 0.86603, 0.     ,
               0.     , 0.86603, 0.86603, 0.86603, 0.86603, 0.     ])
        """
        u = np.diff(self.xy_of_node[self.nodes_at_link], axis=1).reshape((-1, 2))
        return u / np.linalg.norm(u, axis=1).reshape((-1, 1))

    @cached_property
    @read_only_array
    def unit_vector_at_node(self):
        """Get a unit vector for each node.

        Examples
        --------
        >>> from landlab.graph import UniformRectilinearGraph
        >>> graph = UniformRectilinearGraph((3, 3))
        >>> graph.unit_vector_at_node
        array([[1., 1.],
               [2., 1.],
               [1., 1.],
               [1., 2.],
               [2., 2.],
               [1., 2.],
               [1., 1.],
               [2., 1.],
               [1., 1.]])

        >>> from landlab.graph import TriGraph
        >>> graph = TriGraph((3, 2), spacing=2.0, node_layout="hex", sort=True)

        >>> unit_vector_at_node = np.round(graph.unit_vector_at_node, decimals=5)
        >>> unit_vector_at_node[:, 0]
        array([2., 2., 2., 4., 2., 2., 2.])
        >>> unit_vector_at_node[:, 1]
        array([1.73205, 1.73205, 1.73205, 3.4641 , 1.73205, 1.73205, 1.73205])
        """
        unit_vector_at_link = np.vstack((self.unit_vector_at_link, [0.0, 0.0]))
        return np.abs(unit_vector_at_link[self.links_at_node]).sum(axis=1)


class Graph(NetworkGraph):
    """Define the connectivity of a graph of nodes, links, and patches."""

    def __init__(self, node_y_and_x, links=None, patches=None, sort=False):
        if patches is not None and len(patches) == 0:
            patches = None
        self._ds = ugrid_from_unstructured(node_y_and_x, links=links, patches=patches)

        self._frozen = False
        self.freeze()

        if sort:
            Graph.sort(self)

        self._origin = (0.0, 0.0)

    def merge(self, dual, node_at_cell=None, nodes_at_face=None):
        self._dual = dual

        if node_at_cell is not None:
            _update_node_at_cell(self.ds, node_at_cell)
        if nodes_at_face is not None:
            _update_nodes_at_face(self.ds, nodes_at_face)

    def sort(self):
        with self.thawed():
            reorient_link_dirs(self)
            sorted_nodes, sorted_links, sorted_patches = reindex_by_xy(self)
            # reorder_links_at_patch(self)
            if "links_at_patch" in self.ds:
                sort_spokes_at_hub(
                    self.links_at_patch,
                    np.round(self.xy_of_patch, decimals=4),
                    np.round(self.xy_of_link, decimals=4),
                    inplace=True,
                )

        return sorted_nodes, sorted_links, sorted_patches

    @cached_property
    @read_only_array
    def xy_of_patch(self):
        """Get the centroid of each patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.xy_of_patch
        array([[0.5, 0.5],
               [1.5, 0.5]])

        :meta landlab: info-patch
        """
        return get_centroid_of_patch(self)

    @cached_property
    @read_only_array
    def area_of_patch(self):
        """Get the area of each patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.area_of_patch
        array([1.,  1.])

        :meta landlab: info-patch
        """
        return get_area_of_patch(self)

    @property
    def number_of_patches(self):
        """Get the number of patches.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.number_of_patches == 2
        True

        :meta landlab: info-patch
        """
        try:
            return self.ds.sizes["patch"]
        except KeyError:
            return 0

    @property
    def links_at_patch(self):
        """Get the links that define a patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches, sort=True)
        >>> graph.links_at_patch
        array([[3, 5, 2, 0],
               [4, 6, 3, 1]])

        :meta landlab: info-link
        """
        return self.ds["links_at_patch"].values

    @cached_property
    @read_only_array
    def nodes_at_patch(self):
        """Get the nodes that define a patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2])
        >>> links = (
        ...     (0, 1),
        ...     (1, 2),
        ...     (0, 3),
        ...     (1, 4),
        ...     (2, 5),
        ...     (3, 4),
        ...     (4, 5),
        ...     (3, 6),
        ...     (4, 7),
        ...     (5, 8),
        ...     (6, 7),
        ...     (7, 8),
        ... )
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches, sort=True)
        >>> graph.nodes_at_patch
        array([[4, 3, 0, 1],
               [5, 4, 1, 2]])

        :meta landlab: info-node
        """
        nodes_at_patch = get_nodes_at_patch(self)
        sort_spokes_at_hub(
            nodes_at_patch, self.xy_of_patch, self.xy_of_node, inplace=True
        )
        return nodes_at_patch

    @cached_property
    @read_only_array
    def patches_at_node(self):
        """Get the patches that touch each node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches, sort=True)
        >>> graph.patches_at_node
        array([[ 0, -1], [ 1,  0], [ 1, -1],
               [ 0, -1], [ 0,  1], [ 1, -1]])

        :meta landlab: info-patch
        """
        patches_at_node = reverse_one_to_many(self.nodes_at_patch)
        sort_spokes_at_hub(
            patches_at_node, self.xy_of_node, self.xy_of_patch, inplace=True
        )
        return patches_at_node

    @cached_property
    @read_only_array
    def patches_at_link(self):
        """Get the patches on either side of each link.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.patches_at_link
        array([[ 0, -1], [ 1, -1],
               [ 0, -1], [ 0,  1], [ 1, -1],
               [ 0, -1], [ 1, -1]])

        :meta landlab: info-patch
        """
        return reverse_one_to_many(self.links_at_patch, min_counts=2)
