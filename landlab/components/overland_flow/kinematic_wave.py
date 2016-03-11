#!/usr/bin/env python

import numpy as np

from landlab import Component
from ...utils.decorators import use_file_name_or_kwds


class KinematicWave(Component):

    """DESCRIPTION
    This code is based on an overland flow model by Francis Rengers and
    colleagues, , after Julien et al., 1995.
    It was implemented in Landlab by DEJH, March '16. Please cite
    Rengers et al., in review, Model Predictions of Water Runoff in Steep
    Catchments after Wildfire, WRR.

    Construction::

        KinematicWave(grid, )

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.

    Examples
    --------
    >>> from landlab import RasterModelGrid

    """

    _name = 'KinematicWave'

    _input_var_names = (
    )

    _output_var_names = (
    )

    _var_units = {
    }

    _var_mapping = {
    }

    _var_type = {
    }

    _var_doc = {
    }

    @use_file_name_or_kwds
    def __init__(self, grid, mannings_n=0.03,  **kwds):
        """Initialize the kinematic wave approximation overland flow component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.

        """

        self._grid = grid



    def update_one_timestep(self):
        """Update fields with current hydrologic conditions.

        Parameters
        ----------

        """


