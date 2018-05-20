#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 08:01:52 2018

@author: barnhark
"""

#!/usr/env/python

"""
RichDemDepressionFiller.py: Component to accumulate flow and calculate drainage area.

Provides the FlowAccumulator component which accumulates flow and calculates
drainage area. FlowAccumulator supports multiple methods for calculating flow
direction. Optionally a depression finding component can be specified and flow
directing, depression finding, and flow routing can all be accomplished
together.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid # for type tests
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message

from landlab import BAD_INDEX_VALUE
import six
import numpy as np

import richdem as rd


class RichDemDepressionFiller(Component):

    """

    Parameters
    ----------
    grid : ModelGrid
        A grid of type Raster.
    surface : field name at node or array of length node
        The surface to direct flow across.
        Expectation is that this has been pit filled.
    route_method : string
        A string of method. Options are:
        D4
        D8
        Rho8
        Quinn
        Freeman (also needs an exponent)
        Dinf

    Examples
    --------


    """
    _name = 'RichDemDepressionFiller'
    _input_var_names = ()
    _output_var_names = ()
    _var_units = {}
    _var_mapping = {}
    _var_doc = {}

    def __init__(self,
                 grid,
                 surface = 'topographic__elevation',
                 filled_surface = 'topographic__elevation',
                 fill_method='D8',
                 no_data = -9999,
                 **kwargs):
        """
        Initialize the RichDemDepressionFiller component.

        Saves the grid, tests grid type, tests imput types and compatability
        for the flow_director and depression_finder keyword arguments, tests
        the argument of runoff_rate, and initializes new fields.

        modify water__unit_flux_in in order to specify spatially or temporarlly variable
        water unit flux in.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RichDemDepressionFiller
        >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> _ = mg.add_field('topographic__elevation',
        ...                  mg.node_x + mg.node_y,
        ...                  at = 'node')
        """

        super(RichDemDepressionFiller, self).__init__(grid)
        # Keep a local reference to the grid
        self._grid = grid

        # Grid type testing
        if not isinstance(self._grid, RasterModelGrid):
            raise ValueError('RichDemFlowAccumulator only works with RasterModelGrids.')

        self.kwargs = kwargs

        # filled surface
        self.surface = surface
        self.surface_values = return_array_at_node(grid, surface)

        # cell area
        self.node_cell_area = self._grid.cell_area_at_node.copy()
        self.node_cell_area[self._grid.closed_boundary_nodes] = 0.

        try:
            self.water__unit_flux_in = grid.add_ones('water__unit_flux_in', at='node',
                                                      dtype=float)
        except FieldError:
            self.water__unit_flux_in = grid.at_node['water__unit_flux_in']

        try:
            self.drainage_area = grid.add_zeros('drainage_area', at='node',
                                                dtype=float)
        except FieldError:
            self.drainage_area = grid.at_node['drainage_area']

        try:
            self.discharges = grid.add_zeros('surface_water__discharge',
                                             at='node', dtype=float)
        except FieldError:
            self.discharges = grid.at_node['surface_water__discharge']

        try:
            self.upstream_ordered_nodes = grid.add_field('flow__upstream_node_order',
                                                         BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                         at='node', dtype=int)
        except FieldError:
            self.upstream_ordered_nodes = grid.at_node[
                'flow__upstream_node_order']

        # need to specify no data correctly.
        self._rich_dem_array = rd.rdarray(self.surface_values.reshape(grid.shape), no_data=no_data)

        # check that the method works.
        self.method = route_method
        self.route_kwargs = {}

        if route_method == 'Freeman':
            self.route_kwargs['exponent'] = self.kwargs.get('exponent', 4.0)

        # depression/breaching syntax.
        #beau_filled    = rd.FillDepressions(beau, in_place=False)
        #beau_epsilon   = rd.FillDepressions(beau, epsilon=True, in_place=False)
        #beau_breached    = rd.BreachDepressions(beau, in_place=False)

    @property
    def filled_surface(self):
        """Return the drainage area."""
        return self._grid['node']['drainage_area']


    def accumulate_flow(self):
        """
        """
        


    def run_one_step(self):
        """
        """
        pass
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
