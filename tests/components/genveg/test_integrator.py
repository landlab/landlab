import numpy as np

from numpy.testing import assert_allclose
from landlab.components import GenVeg


def test_check_for_dispersal_success(
    example_plant_array, one_cell_grid, example_input_params
):
    # Modify integrator function to accept existing plant arrays or set default to example plant array
    genveg_obj = GenVeg(
        one_cell_grid,
        1,
        180,
        example_input_params,
        plant_array=example_plant_array,
    )
    dispersed_plant_array = genveg_obj.test_for_dispersal_success(example_plant_array)
    # filter = np.nonzero(not np.isnan(example_plant_array["pup_x_loc"]))
    dispersed_size = np.size(dispersed_plant_array, axis=1)
    assert dispersed_size == 13
