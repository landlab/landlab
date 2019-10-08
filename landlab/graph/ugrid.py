import numpy as np
import xarray as xr

from ..utils.jaggedarray import flatten_jagged_array

_MESH_ATTRS = {
    "cf_role": "mesh_topology",
    "long_name": "Topology data of 2D unstructured mesh",
    "topology_dimension": 2,
    "node_coordinates": "x_of_node y_of_node",
    "face_node_connectivity": "nodes_at_patch",
    "face_dimension": "patch",
    "face_edge_connectivity": "links_at_patch",
    "edge_node_connectivity": "nodes_at_link",
    "edge_dimension": "link",
}


def ugrid_from_unstructured(node_y_and_x, links=None, patches=None):
    ugrid = xr.Dataset({"mesh": xr.DataArray(data="a", attrs=_MESH_ATTRS)})

    _update_node_coords(ugrid, node_y_and_x)

    if links is not None:
        _update_nodes_at_link(ugrid, links)

    if patches is not None and "nodes_at_link" in ugrid:
        _update_links_at_patch(ugrid, patches)

    return ugrid


def _update_node_coords(ugrid, node_y_and_x):
    node_y, node_x = (
        np.asarray(node_y_and_x[0], dtype=float),
        np.asarray(node_y_and_x[1], dtype=float),
    )
    y_of_node = xr.DataArray(
        data=node_y.reshape((-1,)),
        coords={"node": np.arange(node_y.size)},
        dims=("node",),
        attrs={"long_name": "y", "units": "m"},
    )
    x_of_node = xr.DataArray(
        data=node_x.reshape((-1,)),
        coords={"node": np.arange(node_x.size)},
        dims=("node",),
        attrs={"long_name": "x", "units": "m"},
    )
    ugrid.update({"y_of_node": y_of_node, "x_of_node": x_of_node})

    return ugrid


def _update_nodes_at_link(ugrid, node_links):
    node_links = np.asarray(node_links, dtype=np.int).reshape((-1, 2))
    nodes_at_link = xr.DataArray(
        data=node_links,
        dims=("link", "Two"),
        attrs={
            "cf_role": "edge_node_connectivity",
            "long_name": "nodes a link tail and head",
            "start_index": 0,
        },
    )
    ugrid.update({"nodes_at_link": nodes_at_link})


def _update_links_at_patch(ugrid, patches):
    from .matrix.at_patch import links_at_patch

    if len(patches) > 0:
        patches = flatten_jagged_array(patches, dtype=int)
    patch_links = links_at_patch(patches)
    links_at_patch = xr.DataArray(
        data=patch_links,
        dims=("patch", "max_patch_links"),
        attrs={
            "cf_role": "face_edge_connectivity",
            "long_name": "Maps every face to its edges",
            "start_index": 0,
        },
    )
    ugrid.update({"links_at_patch": links_at_patch})
