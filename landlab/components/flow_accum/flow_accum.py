#!/usr/env/python

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid  # for type tests
from landlab.utils.decorators import use_file_name_or_kwds


class FlowAccumulator(Component):

    """
    FlowAccumulator Base class Description text here
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
                  'drainage_area'   : 'm**2',
                  'surface_water__discharge': 'm**3/s',
                  'flow__upstream_node_order':'-',
                  'flow__data_structure_delta':'-',
                  'flow__data_structure_D':'-',
                  'flow__nodes_not_in_stack': '-'
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'water__unit_flux_in': 'node',
                    'drainage_area': 'node',
                    'surface_water__discharge': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__nodes_not_in_stack':'grid',
                    'flow__data_structure_delta':'node',
                    'flow__data_structure_D':'link',
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
            'Node array containing the elements delta[1:] of the data structure'
            'delta used for construction of the downstream-to-upstream node'
            'array',
        'flow__data_structure_D':
            'Link array containing the data structure D used for construction'
            'of the downstream-to-upstream node array',
        'flow__nodes_not_in_stack':
            'Boolean value indicating if there are any nodes that have not yet'
            'been added to the stack stored in flow__upstream_node_order.'
            }

   
    @use_file_name_or_kwds
    def __init__(self, grid, surface='topographic__elevation', runoff_rate=None, **kwds):
        # We keep a local reference to the grid
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self.method = None

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
        
        
        # make sure the surface that will be routed/accumulated is present
        # surface must either be a field at node or be array like of size 
        # number_of nodes
        
        if surface!='topographic__elevation' and grid.at_node['topographic__elevation']:
            print ("FlowAccumulator found both the field " +
                   "topographic__elevation' and a provided string or " +
                   "array for the surface argument. THE FIELD "+
                   "topographic__elevation WILL NOT BE USED FOR FLOW " +
                   "DIRECTION OR ACCUMULATION! THE FIELD CODE OR ARRAY" +
                   "GIVEN IN surface WILL BE USED INSTEAD.")
                   
        if surface is str:
            # note that this will test grid.at_node['topographic__elevation']
            # if user doesn't supply a value for surface. 
            grid.at_node[surface]
        else:
            # or be a number_of_nodes sized array
            assert surface.size == grid.number_of_nodes
            
            
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
            self.discharges = grid.add_zeros('surface_water__discharge', at='node',
                                             dtype=float)
        except FieldError:
            self.discharges = grid.at_node['surface_water__discharge']
        
        
        try:
            self.upstream_ordered_nodes = grid.add_zeros(
                'flow__upstream_node_order', at='node', dtype=int)
        except FieldError:
            self.upstream_ordered_nodes = grid.at_node[
                'flow__upstream_node_order']

        try:
            self.delta_structure = grid.add_zeros('flow__data_structure_delta', at='node',
                                                dtype=float)
        except FieldError:
            self.delta_structure = grid.at_node['flow__data_structure_delta']
    
            
        try:
            self.D_structure = grid.add_zeros('flow__data_structure_D', at='link',
                                                dtype=int)
        except FieldError:
            self.D_structure = grid.at_node['flow__data_structure_D']
        
        
        self.nodes_not_in_stack = True 
        
            
    def updated_boundary_conditions(self):
        """
        Call this if boundary conditions on the grid are updated after the
        component is instantiated.
        """
        # We'll also keep track of the active links; if raster, then these are
        # the "D8" links; otherwise, it's just activelinks
        if self._is_raster:
            dal, d8t, d8h = self.grid._d8_active_links()
            self._active_links = dal
            self._activelink_tail = d8t
            self._activelink_head = d8h
            # needs modifying in the loop if D4 (now done)
        else:
            self._active_links = self.grid.active_links
            self._activelink_tail = self.grid.node_at_link_tail[self.grid.active_links]
            self._activelink_head = self.grid.node_at_link_head[self.grid.active_links]

    def run_one_step():
        print ('You have called run_one_step() on the base FlowAccumulator class. ' +
               'This component does not actually accumulate flow, but instead ' +
               'sets up all of the core functionallity of the FlowAccumulators. ' +
               'You probably want to initiate a flow accumulator with a name like:\n'+
               'FlowAccumulator_method\n such as FlowAccumulator_D4 or FlowAccumulator_D8')
              

    if __name__ == '__main__':
        import doctest
        doctest.testmod()


class FlowAccumulator_D4(FlowAccumulator):
    """ 
    Info here
    """
    
    _name = 'FlowAccumulator_D4'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, depressionFinder=None, **kwds):
        super(FlowAccumulator_D4, self).__init__(**kwds)

        # save method as attribute
        self.method = 'D4'
        
        # load correct flow d
        from landlab.flow_direction import FlowDirection_D4 as fd
        from landlab.flow_accum import flow_accum_bw as flow_accum        
        
    def run_one_step():
        
        # step 1. Find flow directions by specified method
        
        # step 2. Stack, D, delta construction
        
        # step 3. Depression finder/router if called for
        
        # step 4. Accumulate (to one or to N depending on direction method. )
    