import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_equal, raises
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab.grid import raster_funcs as rfuncs


def setup_grid():
    """Set up test grid."""
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5),
        'values_at_nodes':  np.arange(20.),
    })


@with_setup(setup_grid)
def test_corner_with_scalar_cell_id():
    """Test using a scalar arg for cell."""
    node_ids = rfuncs.corner_node_at_cell(rmg, np.array([0]), 0)
    assert_array_equal(node_ids, np.array([[2]]))


@with_setup(setup_grid)
def test_face_with_scalar_cell_id():
    """Test using a scalar arg for cell."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0]), 0)
    assert_array_equal(node_ids, np.array([[1]]))


@with_setup(setup_grid)
def test_corner_with_iterable_cell_id():
    """Test using iterable arg for cell."""
    node_ids = rfuncs.corner_node_at_cell(rmg, np.array([0]), (5, ))
    assert_array_equal(node_ids, np.array([[9]]))


@with_setup(setup_grid)
def test_face_with_iterable_cell_id():
    """Test using iterable arg for cell."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0]), (5, ))
    assert_array_equal(node_ids, np.array([[8]]))


@with_setup(setup_grid)
def test_corner_with_no_cell_id():
    """Test without an arg for cells to mean all cells."""
    node_ids = rfuncs.corner_node_at_cell(rmg, np.array([0]))
    assert_array_equal(node_ids, np.array([[2, 3, 4, 7, 8, 9]]).T)


@with_setup(setup_grid)
def test_face_with_no_cell_id():
    """Test without an arg for cells to mean all cells."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0]))
    assert_array_equal(node_ids, np.array([[1, 2, 3, 6, 7, 8]]).T)


@with_setup(setup_grid)
def test_corner_multiple_corners():
    """Test getting nodes for more than one corner."""
    node_ids = rfuncs.corner_node_at_cell(rmg, np.array([0, 1]), (4, ))
    assert_array_equal(node_ids, np.array([[8, 6]]))


@with_setup(setup_grid)
def test_face_multiple_faces():
    """Test getting nodes for more than one corner."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0, 1]), (4, ))
    assert_array_equal(node_ids, np.array([[7, 11]]))


@with_setup(setup_grid)
def test_corner_multiple_corners_and_cells():
    """Test getting nodes for more than one corner and more than one cell."""
    node_ids = rfuncs.corner_node_at_cell(rmg, np.array([0, 1]), (4, 5))
    assert_array_equal(node_ids, np.array([[8, 6], [9, 7]]))


@with_setup(setup_grid)
def test_face_multiple_corners_and_cells():
    """Test getting nodes for more than one corner and more than one cell."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, np.array([0, 1]), (4, 5))
    assert_array_equal(node_ids, np.array([[7, 11], [8, 12]]))


@with_setup(setup_grid)
def test_corner_type_tuple():
    """Test using tuple as corner arg."""
    node_ids = rfuncs.corner_node_at_cell(rmg, (0, 1), (4, 5))
    assert_array_equal(node_ids, np.array([[8, 6], [9, 7]]))


@with_setup(setup_grid)
def test_face_type_tuple():
    """Test using tuple as face arg."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, (0, 1), (4, 5))
    assert_array_equal(node_ids, np.array([[7, 11], [8, 12]]))


@with_setup(setup_grid)
def test_corner_type_list():
    """Test using list as corner arg."""
    node_ids = rfuncs.corner_node_at_cell(rmg, [0, 1], (4, 5))
    assert_array_equal(node_ids, np.array([[8, 6], [9, 7]]))


@with_setup(setup_grid)
def test_face_type_list():
    """Test using list as face arg."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, [0, 1], (4, 5))
    assert_array_equal(node_ids, np.array([[7, 11], [8, 12]]))


@with_setup(setup_grid)
def test_corner_type_scalar():
    """Test using scalar as corner arg."""
    node_ids = rfuncs.corner_node_at_cell(rmg, 2, (4, ))
    assert_array_equal(node_ids, np.array([16]))

    node_ids = rfuncs.corner_node_at_cell(rmg, 2, (4, 5))
    assert_array_equal(node_ids, np.array([16, 17]))

    node_ids = rfuncs.corner_node_at_cell(rmg, 1)
    assert_array_equal(node_ids, np.array([0, 1, 2, 5, 6, 7]))


@with_setup(setup_grid)
def test_face_type_scalar():
    """Test using scalar as face arg."""
    node_ids = rfuncs.neighbor_node_at_cell(rmg, 2, (4, ))
    assert_array_equal(node_ids, np.array([17]))

    node_ids = rfuncs.neighbor_node_at_cell(rmg, 2, (4, 5))
    assert_array_equal(node_ids, np.array([17, 18]))

    node_ids = rfuncs.neighbor_node_at_cell(rmg, 1)
    assert_array_equal(node_ids, np.array([5, 6, 7, 10, 11, 12]))
