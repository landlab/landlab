import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid import raster_funcs as rfuncs


def test_face_with_scalar_cell_id():
    """Test using a scalar arg for cell."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0]), 0)
    assert_array_equal(node_ids, np.array([[1]]))


def test_face_with_iterable_cell_id():
    """Test using iterable arg for cell."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0]), (5,))
    assert_array_equal(node_ids, np.array([[8]]))


def test_face_with_no_cell_id():
    """Test without an arg for cells to mean all cells."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0]))
    assert_array_equal(node_ids, np.array([[1, 2, 3, 6, 7, 8]]).T)


def test_face_multiple_faces():
    """Test getting nodes for more than one corner."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0, 1]), (4,))
    assert_array_equal(node_ids, np.array([[7, 11]]))


def test_face_multiple_corners_and_cells():
    """Test getting nodes for more than one corner and more than one cell."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0, 1]), (4, 5))
    assert_array_equal(node_ids, np.array([[7, 11], [8, 12]]))


def test_face_type_tuple():
    """Test using tuple as face arg."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, (0, 1), (4, 5))
    assert_array_equal(node_ids, np.array([[7, 11], [8, 12]]))


def test_face_type_list():
    """Test using list as face arg."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, [0, 1], (4, 5))
    assert_array_equal(node_ids, np.array([[7, 11], [8, 12]]))


def test_face_type_scalar():
    """Test using scalar as face arg."""
    rmg = RasterModelGrid((4, 5))
    node_ids = rfuncs.neighbor_node_at_cell(rmg, 2, (4,))
    assert_array_equal(node_ids, np.array([17]))

    node_ids = rfuncs.neighbor_node_at_cell(rmg, 2, (4, 5))
    assert_array_equal(node_ids, np.array([17, 18]))

    node_ids = rfuncs.neighbor_node_at_cell(rmg, 1)
    assert_array_equal(node_ids, np.array([5, 6, 7, 10, 11, 12]))
