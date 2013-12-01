import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup, raises

from landlab import RasterModelGrid


def setup():
    globals()['rmg'] = RasterModelGrid(4, 4)


@with_setup(setup)
def test_unit_grid_all_cells():
    assert_array_equal(rmg.cell_areas, np.ones(4))


@with_setup(setup)
def test_unit_grid_one_cell():
    assert_array_equal(rmg.cell_areas[0], 1.)


@with_setup(setup)
def test_unit_grid_last_cell():
    assert_array_equal(rmg.cell_areas[-1], 1.)


@with_setup(setup)
@raises(IndexError)
def test_out_of_range():
    rmg.cell_areas[5]


@with_setup(setup)
@raises(ValueError)
def test_is_immutable():
    rmg.cell_areas[0] = 0.


def test_all_cells_with_spacing():
    rmg = RasterModelGrid(4, 4, 10.)
    assert_array_equal(rmg.cell_areas, 100 * np.ones(4))
