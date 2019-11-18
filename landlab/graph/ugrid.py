import numpy as np
import xarray as xr

from ..utils.jaggedarray import flatten_jagged_array

MESH_ATTRS = {
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

DUAL_MESH_ATTRS = {
    "cf_role": "mesh_topology",
    "long_name": "Topology data of 2D unstructured mesh",
    "topology_dimension": 2,
    "node_coordinates": "x_of_corner y_of_corner",
    "face_node_connectivity": "corners_at_cell",
    "face_dimension": "cell",
    "face_edge_connectivity": "faces_at_cell",
    "edge_node_connectivity": "corners_at_face",
    "edge_dimension": "face",
}


def ugrid_from_structured_quad(coords, shape=None):
    from .structured_quad.structured_quad import (
        setup_nodes_at_link,
        setup_links_at_patch,
        setup_node_coords_structured,
    )

    node_y_and_x = setup_node_coords_structured(coords, shape=shape)
    shape = node_y_and_x[0].shape

    return ugrid_from_unstructured(
        node_y_and_x, setup_nodes_at_link(shape), setup_links_at_patch(shape)
    )


def ugrid_from_rectilinear(coords):
    from .structured_quad.structured_quad import setup_node_coords_rectilinear

    node_y_and_x = setup_node_coords_rectilinear(coords)
    return ugrid_from_structured_quad(node_y_and_x)


def ugrid_from_uniform_rectilinear(shape, spacing=1.0, origin=0.0):
    from .structured_quad.structured_quad import setup_node_coords

    node_y_and_x = setup_node_coords(shape, spacing, origin)
    return ugrid_from_structured_quad(node_y_and_x)


def load_ugrid(ugrid):
    if not isinstance(ugrid, xr.Dataset):
        raise AssertionError("not an instance of xarray.Dataset")

    ds = xr.Dataset({"mesh": xr.DataArray(data="a", attrs=MESH_ATTRS)})

    meta = ugrid.a.attrs

    node_coords = meta["node_coordinates"].split()
    ds.update({"x_of_node": ugrid[node_coords[0]], "y_of_node": ugrid[node_coords[1]]})
    ds = ds.rename({ds.x_of_node.dims[0]: "node"})

    ds.update({"nodes_at_link": ugrid[meta["edge_node_connectivity"]]})
    ds = ds.rename({ds.nodes_at_link.dims[0]: "link", ds.nodes_at_link.dims[1]: "Two"})

    ds.update({"links_at_patch": ugrid[meta["face_edge_connectivity"]]})
    ds = ds.rename(
        {
            ds.links_at_patch.dims[0]: "patch",
            ds.links_at_patch.dims[1]: "max_patch_links",
        }
    )

    return ds


def ugrid_as_dual(ugrid):
    rename = {
        "node": "corner",
        "link": "face",
        "patch": "cell",
        "x_of_node": "x_of_corner",
        "y_of_node": "y_of_corner",
        "nodes_at_link": "corners_at_face",
        "links_at_patch": "faces_at_cell",
        "max_patch_links": "max_cell_faces",
    }
    return ugrid.rename(rename)


def ugrid_from_unstructured(node_y_and_x, links=None, patches=None):
    ugrid = xr.Dataset({"mesh": xr.DataArray(data="a", attrs=MESH_ATTRS)})

    update_node_coords(ugrid, node_y_and_x)

    if links is not None:
        update_nodes_at_link(ugrid, links)

    if patches is not None and "nodes_at_link" in ugrid:
        update_links_at_patch(ugrid, patches)

    return ugrid


def update_node_coords(ugrid, node_y_and_x):
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


def update_nodes_at_link(ugrid, node_links):
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


def update_links_at_patch(ugrid, patches):
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


def update_nodes_at_patch(ugrid, nodes_at_patch):
    nodes_at_patch = xr.DataArray(
        data=nodes_at_patch,
        dims=("patch", "max_patch_links"),
        attrs={
            "cf_role": "patch_node_connectivity",
            "long_name": "nodes defining patches",
            "start_index": 0,
        },
    )
    ugrid.update({"nodes_at_patch": nodes_at_patch})


def update_node_at_cell(ugrid, node_at_cell):
    node_at_cell = xr.DataArray(
        data=node_at_cell,
        dims=("cell",),
        attrs={
            "cf_role": "node_cell_connectivity",
            "long_name": "node at cell",
            "start_index": 0,
        },
    )
    ugrid.update({"node_at_cell": node_at_cell})


def update_nodes_at_face(ugrid, nodes_at_face):
    nodes_at_face = xr.DataArray(
        data=nodes_at_face,
        dims=("face", "Two"),
        attrs={
            "cf_role": "node_face_connectivity",
            "long_name": "nodes at face",
            "start_index": 0,
        },
    )
    ugrid.update({"nodes_at_face": nodes_at_face})
