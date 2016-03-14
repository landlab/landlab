#!/usr/bin/env python

import numpy as np

from landlab import Component, RasterModelGrid
from ...utils.decorators import use_file_name_or_kwds


class KinematicWave(Component):

    """
    This code is based on an overland flow model by Francis Rengers and
    colleagues, after Julien et al., 1995.
    It was implemented in Landlab by DEJH, March '16. Please cite
    Rengers et al., in review, Model Predictions of Water Runoff in Steep
    Catchments after Wildfire, WRR.

    Note this module assumes that the topography DOES NOT change during the
    run. If it does, call XXXXXXXX to update the component to the new topo.

    Construction::

        KinematicWave(grid, mannings_n=0.03)

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    mannings_n : float
        A value to use for Manning's n in the Manning discharge equation.

    Examples
    --------
    >>> from landlab import RasterModelGrid

    """

    _name = 'KinematicWave'

    _input_var_names = (
        'topographic__elevation',
        'surface_water__depth',
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
        """

        assert isinstance(grid, RasterModelGrid), 'grid must be regular'
        self._grid = grid
        active = self.grid.status_at_node != CLOSED_BOUNDARY
        self._active_depths = mg.at_node['surface_water__depth'][active]
        all_grads = self.grid.calculate_gradients_at_links(
            'topographic__elevation')


    def update_one_timestep(self, update_topography=False):
        """Update fields with current hydrologic conditions.

        Parameters
        ----------

        """


