import numpy as np

from landlab.components import GenVeg


def test_check_for_dispersal_success(
    example_plant_array, one_cell_grid, example_input_params
):
    # Modify integrator function to accept existing plant arrays or set default to example plant array
    genveg_obj = GenVeg(
        one_cell_grid,
        1,
        np.datetime64("2010-06-28"),
        example_input_params,
        plant_array=example_plant_array,
    )
    dispersed_plant_array = genveg_obj.check_for_dispersal_success(example_plant_array)
    # filter = np.nonzero(not np.isnan(example_plant_array["pup_x_loc"]))
    (dispersed_size,) = dispersed_plant_array.shape
    return example_plant_array.shape
