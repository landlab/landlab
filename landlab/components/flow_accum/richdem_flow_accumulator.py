#!/usr/env/python

"""
Component to accumulate flow and calculate drainage area using RichDem XXX Link here.

Provides the RichDemFlowAccumulator component which accumulates flow and
calculates drainage area. This component only works with RasterModelGrids and
supports those options available in the RichDem package XXX link here.
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

_VALID_METHODS = ['D4',
                  'D8',
                  'Rho8',
                  'Quinn',
                  'Freeman',
                  'Dinf']


class RichDemFlowAccumulator(Component):

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
    kwargs : additional keyword arguments
        Some RichDem methods take additional arguments. They are as follows:
        The 'Freeman' method takes a 'exponent' argument (default value 4.0).

    Examples
    --------


    """
    _name = 'RichDemFlowAccumulator'
    _input_var_names = ()
    _output_var_names = ()
    _var_units = {}
    _var_mapping = {}
    _var_doc = {}

    def __init__(self,
                 grid,
                 surface = 'topographic__elevation',
                 route_method='D8',
                 no_data = -9999,
                 **kwargs):
        """
        Initialize the RichDemFlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and compatability
        for the flow_director and depression_finder keyword arguments, tests
        the argument of runoff_rate, and initializes new fields.

        If spatially variable water unit flux in is desired, then the user
        must modify the at-node model grid field ``water__unit_flux_in``. This
        is multiplied by the cell area and accumulated.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RichDemFlowAccumulator
        >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> _ = mg.add_field('topographic__elevation',
        ...                  mg.node_x + mg.node_y,
        ...                  at = 'node')
        """
        super(RichDemFlowAccumulator, self).__init__(grid)

        # Keep a local reference to the grid
        self._grid = grid

        # check that the method works.
        if route_method in _VALID_METHODS:
            self.method = route_method
        else:
            raise ValueError('RichDemFlowAccumulator: route_method invalid.')

        # Grid type testing
        if not isinstance(self._grid, RasterModelGrid):
            raise ValueError('RichDemFlowAccumulator only works with RasterModelGrids.')

        # get kwargs and no data.
        self._kwargs = kwargs
        self._no_data = no_data

        # get the surface to route flow on.
        self.surface = surface
        self.surface_values = return_array_at_node(grid, surface)

        # get and save status at node
        self._status = grid.status_at_node

        # cell area
        self._node_cell_area = self._grid.cell_area_at_node.copy()
        self._node_cell_area[self._grid.closed_boundary_nodes] = 0.

        # setup required fields.
        self._setup_accumulator_fields()

        # get the routing kwargs
        self._route_kwargs = {}

        # set routing kwargs for Freeman.
        if route_method == 'Freeman':
            self._route_kwargs['exponent'] = self._kwargs.get('exponent', 4.0)

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid['node']['drainage_area']

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid['node']['surface_water__discharge']

    @property
    def node_order_upstream(self):
        """Return the upstream node order (drainage stack)."""
        return self._grid['node']['flow__upstream_node_order']

    def _setup_accumulator_fields(self):
        if 'water__unit_flux_in' in grid.at_node:
            self.water__unit_flux_in = grid.at_node['water__unit_flux_in']
        else:
            self.water__unit_flux_in = grid.add_ones('water__unit_flux_in', at='node',
                                                      dtype=float)

        if 'drainage_area' in grid.at_node:
            self.drainage_area = grid.at_node['drainage_area']
        else:
            self.drainage_area = grid.add_zeros('drainage_area', at='node',
                                                dtype=float)

        if 'surface_water__discharge' in grid.at_node:
            self.discharges = grid.at_node['surface_water__discharge']
        else:
            self.discharges = grid.add_zeros('surface_water__discharge',
                                             at='node', dtype=float)

        if 'flow__upstream_node_order' in grid.at_node:
            self.upstream_ordered_nodes = grid.at_node['flow__upstream_node_order']
        else:
            self.upstream_ordered_nodes = grid.add_field('flow__upstream_node_order',
                                                         BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                         at='node', dtype=int)

    def accumulate_flow(self):
        """
        """
        # need to specify no data correctly.
        vals = self.surface_values.copy()
        vals[self._status == 4] = self._no_data
        self._rich_dem_array = rd.rdarray(vals.reshape(self._grid.shape), no_data=self.no_data)

        # some things to figure out:
        # 1) How to deal with boundary conditions in RichDem (e.g. closed or
        # fixed boundary nodes). We don't want to set those to no-data and I
        # don't think that setting them to have zero area will be sufficient.

        # run twice if water__unit_flux_in is not  1.0 everywhere.
        a = rd.FlowAccumulation(self._rich_dem_array,
                                method=self.method,
                                weights = self._node_cell_area.reshape(self._grid.shape),
                                **self._route_kwargs).flatten()
        a[self._status == 4] = 0.0
        if np.any(self.water__unit_flux_in != 1.0) :
            q = rd.FlowAccumulation(self._rich_dem_array,
                                    method=self.method,
                                    weights = (self.water__unit_flux_in.reshape(self._grid.shape) *
                                               self._node_cell_area.reshape(self._grid.shape)),
                                    **self.route_kwargs).flatten()
            q[self._status == 4] = 0
        else:
            q = a

        #self._grid['node']['flow__upstream_node_order'][:] = s
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['surface_water__discharge'][:] = q

        return (a, q)

    def run_one_step(self):
        """
        """
        self.accumulate_flow()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
