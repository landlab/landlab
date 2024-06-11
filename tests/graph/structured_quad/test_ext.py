import numpy as np
import pytest

from landlab.graph.structured_quad.ext.at_cell import fill_node_at_cell
from landlab.graph.structured_quad.ext.at_face import fill_nodes_at_face
from landlab.graph.structured_quad.ext.at_link import fill_nodes_at_link
from landlab.graph.structured_quad.ext.at_link import fill_patches_at_link
from landlab.graph.structured_quad.ext.at_node import fill_links_at_node
from landlab.graph.structured_quad.ext.at_node import fill_patches_at_node
from landlab.graph.structured_quad.ext.at_patch import fill_links_at_patch


def number_of_nodes(shape):
    return shape[0] * shape[1]


def number_of_links(shape):
    return shape[0] * (shape[1] - 1) + (shape[0] - 1) * shape[1]


def number_of_patches(shape):
    return (shape[0] - 1) * (shape[1] - 1)


def number_of_corners(shape):
    return number_of_nodes((shape[0] - 1, shape[1] - 1))


def number_of_faces(shape):
    return number_of_links((shape[0] - 1, shape[1] - 1))


def number_of_cells(shape):
    return number_of_patches((shape[0] - 1, shape[1] - 1))


@pytest.mark.parametrize("n_rows", (3, 4, 8, 512))
@pytest.mark.parametrize("n_cols", (3, 4, 8, 512))
@pytest.mark.parametrize("dtype", (np.int64, np.longlong))
@pytest.mark.parametrize(
    "func,shape",
    (
        (fill_node_at_cell, lambda shape: (number_of_cells(shape), 1)),
        (fill_nodes_at_face, lambda shape: (number_of_faces(shape), 2)),
        (fill_patches_at_link, lambda shape: (number_of_links(shape), 2)),
        (fill_nodes_at_link, lambda shape: (number_of_links(shape), 2)),
        (fill_patches_at_node, lambda shape: (number_of_nodes(shape), 4)),
        (fill_links_at_node, lambda shape: (number_of_nodes(shape), 4)),
        (fill_links_at_patch, lambda shape: (number_of_patches(shape), 4)),
    ),
)
def test_fill_over_under_flow(n_rows, n_cols, dtype, func, shape):
    """Test buffer over/under flow."""
    buffer_size = 1024

    shape = shape((n_rows, n_cols))

    actual = np.full((shape[0] + 2 * buffer_size, shape[1]), -2, dtype=dtype).squeeze()

    shape = (n_rows, n_cols)

    func(shape, actual[buffer_size:-buffer_size])

    assert np.all(actual[:buffer_size] == -2)
    assert np.all(actual[-buffer_size:] == -2)
    assert np.all(actual[buffer_size:-buffer_size] >= -1)
