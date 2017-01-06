#!/usr/env/python

from __future__ import print_function

import warnings

import landlab
from landlab import FieldError, Component
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.decorators import use_field_name_or_array
from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY, \
    BAD_INDEX_VALUE
import numpy as np


@use_field_name_or_array('node')
def return_surface(grid, surface):
    return(surface)


class FlowAccumulator(Component):

    """
    Base class for accumulating flow.

    This component is not meant to be used directly in modeling efforts.
    Instead it has the functionality that all flow accumulators need
    to initialize, check boundary conditions, and create all necesary fields.

    Stores as ModelGrid fields:

        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure:
            *flow__data_structure_delta*. The first element is always zero.
        -  Link array of the D data structure: *flow__data_structure_D

    The primary method of this class, :func:`run_one_step` is not implemented.


    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.
    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_accum.flow_accumulator import FlowAccumulator
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y, at='node')
    >>> fa = FlowAccumulator(mg, 'topographic__elevation')
    >>> fa.elevs
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> sorted(list(mg.at_node.keys()))
    ['drainage_area', 'flow__data_structure_delta',
       'flow__upstream_node_order', 'surface_water__discharge',
       'topographic__elevation', 'water__unit_flux_in']


    """
    _name = 'FlowAccumulator'

    _input_var_names = ('topographic__elevation',
                        'water__unit_flux_in'
                        )

    _output_var_names = ('drainage_area',
                         'surface_water__discharge',
                         'flow__upstream_node_order',
                         'flow__nodes_not_in_stack',
                         'flow__data_structure_delta',
                         'flow__data_structure_D'
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_node': 'm',
                  'water__unit_flux_in': 'm/s',
                  'drainage_area': 'm**2',
                  'surface_water__discharge': 'm**3/s',
                  'flow__upstream_node_order': '-',
                  'flow__data_structure_delta': '-',
                  'flow__data_structure_D': '-',
                  'flow__nodes_not_in_stack': '-'
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'water__unit_flux_in': 'node',
                    'drainage_area': 'node',
                    'surface_water__discharge': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__nodes_not_in_stack': 'grid',
                    'flow__data_structure_delta': 'node',
                    'flow__data_structure_D': 'link',
                    }
    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'surface_water__discharge': 'Discharge of water through each node',
        'water__unit_flux_in':
            'External volume water per area per time input to each node '
            '(e.g., rainfall rate)',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'flow__data_structure_delta':
            'Node array containing the elements delta[1:] of the data '
            'structure "delta" used for construction of the downstream-to-'
            'upstream node array',
        'flow__data_structure_D':
            'Link array containing the data structure D used for construction'
            'of the downstream-to-upstream node array',
        'flow__nodes_not_in_stack':
            'Boolean value indicating if there are any nodes that have not yet'
            'been added to the stack stored in flow__upstream_node_order.'
            }

    def __init__(self, grid, surface, runoff_rate=None):
        # We keep a local reference to the grid
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if hasattr(self, 'method') == False:
            self.method = 'base'
            
        
        self.updated_boundary_conditions()

        # START: Testing of input values, supplied either in function call or
        # as part of the grid.

        # testing input for runoff rate, can be None, a string associated with
        # a field at node, a single float or int, or an array of size number of
        # nodes.
        if runoff_rate is not None:
            if type(runoff_rate) is str:
                runoff_rate = grid.at_node[runoff_rate]
            elif type(runoff_rate) in (float, int):
                pass
            else:
                assert runoff_rate.size == grid.number_of_nodes

        # test for water__unit_flux_in
        try:
            grid.at_node['water__unit_flux_in']
        except FieldError:
            if runoff_rate is None:
                # assume that if runoff rate is not supplied, that the value
                # should be set to one everywhere.
                grid.add_ones('node', 'water__unit_flux_in', dtype=float)
            else:
                if type(runoff_rate) in (float, int):
                    grid.add_empty('node', 'water__unit_flux_in', dtype=float)
                    grid.at_node['water__unit_flux_in'].fill(runoff_rate)
                else:
                    grid.at_node['water__unit_flux_in'] = runoff_rate
        else:
            if runoff_rate is not None:
                print ("FlowAccumulator found both the field " +
                       "'water__unit_flux_in' and a provided float or " +
                       "array for the runoff_rate argument. THE FIELD IS " +
                       "BEING OVERWRITTEN WITH THE SUPPLIED RUNOFF_RATE!")
                if type(runoff_rate) in (float, int):
                    grid.at_node['water__unit_flux_in'].fill(runoff_rate)
                else:
                    grid.at_node['water__unit_flux_in'] = runoff_rate

        # perform a test (for politeness!) that the old name for the water_in
        # field is not present:
        try:
            grid.at_node['water__discharge_in']
        except FieldError:
            pass
        else:
            warnings.warn("This component formerly took 'water__discharge" +
                          "_in' as an input field. However, this field is " +
                          "now named 'water__unit_flux_in'. You are still " +
                          "using a field with the old name. Please update " +
                          "your code if you intended the FlowRouter to use " +
                          "that field.", DeprecationWarning)

        # save elevations and node_cell_area to class properites.
        self.surface = surface
        surf = return_surface(grid, surface)

        # add elevations as a local variable.
        self.elevs = surf

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.

        self.node_cell_area = node_cell_area

        # This component will track of the following variables.
        # Attempt to create each, if they already exist, assign the existing
        # version to the local copy.

        #   - drainage area at each node
        #   - receiver of each node
        #   - D array
        #   - delta array
        #   - missing nodes in stack.
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
            self.upstream_ordered_nodes = grid.add_field(
                'flow__upstream_node_order',
                BAD_INDEX_VALUE*grid.ones(at='node'),
                at='node', dtype=int)

        except FieldError:
            self.upstream_ordered_nodes = grid.at_node[
                'flow__upstream_node_order']

        try:
            self.delta_structure = grid.add_field(
                'flow__data_structure_delta',
                BAD_INDEX_VALUE*grid.ones(at='node'),
                at='node', dtype=int)
        except FieldError:
            self.delta_structure = grid.at_node['flow__data_structure_delta']

        try:
            # needs to be BAD_INDEX_VALUE
            self.D_structure = grid.add_field(
                'flow__data_structure_D',
                BAD_INDEX_VALUE*grid.ones(at='link'),
                at='link', dtype=int)

        except FieldError:
            self.D_structure = grid.at_link['flow__data_structure_D']

        self.nodes_not_in_stack = True

    def updated_boundary_conditions(self):
        """
        Call this if boundary conditions on the grid are updated after the
        component is instantiated.
        """
        # We'll also keep track of the active links; if raster, then these are
        # the "D8" links; otherwise, it's just activelinks
        if self.method == 'D8':
            dal, d8t, d8h = self.grid._d8_active_links()
            self._active_links = dal
            self._activelink_tail = d8t
            self._activelink_head = d8h
            # needs modifying in the loop if D4 (now done)
        else:
            self._active_links = self.grid.active_links
            self._activelink_tail = self.grid.node_at_link_tail[
                self.grid.active_links]
            self._activelink_head = self.grid.node_at_link_head[
                self.grid.active_links]

    def run_one_step(self):
        raise NotImplementedError('run_one_step()')

    @property
    def node_drainage_area(self):
        return self._grid['node']['drainage_area']

    @property
    def node_water_discharge(self):
        return self._grid['node']['surface_water__discharge']

    @property
    def node_order_upstream(self):
        return self._grid['node']['flow__upstream_node_order']

    @property
    def node_D_structure(self):
        return self._grid['node']['flow__data_structure_D']

    @property
    def node_delta_structure(self):
        return self._grid['node']['flow__data_structure_delta']

if __name__ == '__main__':
    import doctest
    doctest.testmod()
