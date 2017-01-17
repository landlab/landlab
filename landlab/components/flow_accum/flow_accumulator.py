#!/usr/env/python

from __future__ import print_function

import warnings

import landlab
from landlab import FieldError, Component
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.decorators import use_field_name_or_array

from landlab.components.flow_accum import flow_accum_bw
from landlab.components.flow_accum import flow_accum_to_n

from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY, \
    BAD_INDEX_VALUE
import numpy as np
import six

@use_field_name_or_array('node')
def return_surface(grid, surface):
    return(surface)

class FlowAccumulator(Component):

    """
    Component to accumulate flow and calculate drainage area. 
 
    This is accomplished by first finding flow directions by a user-specified
    method and then calculating the drainage area and discharge. 
    
    Optionally, spatially variable runoff can be set either by the model grid 
    field 'water__unit_flux_in' or the input variable *runoff_rate**.
    
    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished 
    together. 
    
    
    NOTE: The perimeter nodes  NEVER contribute to the accumulating flux, even 
    if the  gradients from them point inwards to the main body of the grid. 
    This is because under Landlab definitions, perimeter nodes lack cells, so 
    cannot accumulate any discharge.


    FlowAccumulator stores as ModelGrid fields:
               
        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure: 
            *flow__data_structure_delta*. The first element is always zero.
        -  Link array of the D data structure: *flow__data_structure_D*
        
    The FlowDirector component will add additional ModelGrid fields. 
    DirectToOne methods(Steepest/D4 and D8) and DirectToMany(NAMES HERE) use
    different model grid fields. 
    
    DirectToOne Methods (Steeptest/D4 and D8) store the following as ModelGrid 
    fields:
        
        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Node array containing ID of link that leads from each node to its
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    DirectToMany Methods (NAMES HERE) store the following as ModelGrid 
    fields:
        
        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_nodes'*
        -  Node array of receiver proportions: *'flow__receiver_proportions'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*
        
    The primary method of this class is :func:`run_one_step`

    Parameters
    ----------
    grid : ModelGrid
        A grid of type Voroni.
    surface : field name at node or array of length node
        The surface to direct flow across.  
    flow_director : string, class, instance of class. 
        A string of method or class name (e.g. 'D8' or 'FlowDirectorD8'), an 
        uninstantiated FlowDirector class, or an instance of a FlowDirector 
        class. This sets the method used to calculate flow directions. 
        Default is 'FlowDirectorSteepest'
    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.  
    depression_finder : string, class, instance of class, optional
         A string of class name (e.g., 'DepressionFinderAndRouter'), an 
         uninstantiated DepressionFinder class, or an instance of a 
         DepressionFinder class. 
         This sets the method for depression finding. 
        

    Examples
    --------

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

    def __init__(self, grid, surface = 'topographic__elevation', flow_director = 'FlowDirectorD4', runoff_rate=None, depression_finder = None):
        # Keep a local reference to the grid
        self._grid = grid

        # Grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)

        # STEP 1: Testing of input values, supplied either in function call or
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
            
        # STEP 2: 
        # identify Flow Director method, save name, import and initialize the correct
        # flow director component if necessary        
        PERMITTED_DIRECTORS = ['FlowDirectorSteepest',
                               'FlowDirectorD8']
        
        PERMITTED_DEPRESSION_FINDERS = ['DepressionFinderAndRouter']

        # flow director is provided as a string.            
        if isinstance(flow_director, six.string_types):
            if flow_director[:12] == 'FlowDirector':
                flow_director = flow_director[12:]
            
            from landlab.components.flow_director import FlowDirectorSteepest, FlowDirectorD8
            DIRECTOR_METHODS = {'D4': FlowDirectorSteepest,
                            'Steepest': FlowDirectorSteepest,
                            'D8': FlowDirectorD8,
                            }
                
            try:
                FlowDirector = DIRECTOR_METHODS[flow_director]
            except KeyError:
                raise ValueError('String provided in flow_director is not a valid method or component name. '
                                 'The following components are valid imputs:\n'+str(PERMITTED_DIRECTORS))
                
            self.fd = FlowDirector(self._grid, self.elevs)
        # flow director is provided as an instantiated flow director   
        elif isinstance(flow_director, Component):
             if flow_director._name in PERMITTED_DIRECTORS:
                 self.fd = flow_director
             else:
                 raise ValueError('Component provided in flow_director is not a valid component. '
                                 'The following components are valid imputs:\n'+str(PERMITTED_DIRECTORS))
        # flow director is provided as an uninstantiated flow director 
        else:
            
            if flow_director._name in PERMITTED_DIRECTORS:
                FlowDirector = flow_director
                self.fd = FlowDirector(self._grid, self.elevs)
            else:
                raise ValueError('Component provided in flow_director is not a valid component. '
                                 'The following components are valid imputs:\n'+str(PERMITTED_DIRECTORS))
                           
        # save method as attribute    
        self.method = self.fd.method
                
        # now do a similar thing for the depression finder. 
        self.depression_finder = depression_finder
        if self.depression_finder:
            # depression finder is provided as a string.            
            if isinstance(self.depression_finder, six.string_types):
                
                from landlab.components import DepressionFinderAndRouter
                DEPRESSION_METHODS = {'DepressionFinderAndRouter': DepressionFinderAndRouter
                                    }
                    
                try:
                    DepressionFinder = DEPRESSION_METHODS[self.depression_finder]
                except KeyError:
                    raise ValueError('Component provided in depression_finder is not a valid component. '
                                     'The following components are valid imputs:\n'+str(PERMITTED_DEPRESSION_FINDERS))
                    
                self.df = DepressionFinder(self._grid)
            # flow director is provided as an instantiated depression finder   
            elif isinstance(self.depression_finder, Component):  
                
                if self.depression_finder._name in PERMITTED_DEPRESSION_FINDERS:
                    self.df = self.depression_finder
                else:
                    raise ValueError('Component provided in depression_finder is not a valid component. '
                                     'The following components are valid imputs:\n'+str(PERMITTED_DEPRESSION_FINDERS))
            # depression_fiuner is provided as an uninstantiated depression finder
            else:
                                       
                if self.depression_finder._name in PERMITTED_DEPRESSION_FINDERS:
                    DepressionFinder = self.depression_finder
                    self.df = DepressionFinder(self._grid)
                else:
                    raise ValueError('Component provided in depression_finder is not a valid component. '
                    'The following components are valid imputs:\n'+str(PERMITTED_DEPRESSION_FINDERS))
 
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
            self.upstream_ordered_nodes = grid.add_field('flow__upstream_node_order',
                                                         BAD_INDEX_VALUE*grid.ones(at='node'),
                                                         at='node', dtype=int)
            
        except FieldError:
            self.upstream_ordered_nodes = grid.at_node[
                'flow__upstream_node_order']

        try:
            self.delta_structure = grid.add_field('flow__data_structure_delta',
                                                  BAD_INDEX_VALUE*grid.ones(at='node'),
                                                  at='node', dtype=int)
        except FieldError:
            self.delta_structure = grid.at_node['flow__data_structure_delta']

        try:
            # needs to be BAD_INDEX_VALUE
            self.D_structure = grid.add_field('flow__data_structure_D',
                                              BAD_INDEX_VALUE*grid.ones(at='link'),
                                              at='link', dtype=int)
        except FieldError:
            self.D_structure = grid.at_link['flow__data_structure_D']

        self.nodes_not_in_stack = True


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


    def run_one_step(self):
        
        # step 1. Find flow directions by specified method
        self.fd.run_one_step()
        
        # further steps vary depending on how many recievers are present
        # one set of steps is for route to one (D8, Steepest/D4)
        if self.fd.to_n_receivers == 'one':
        
            # step 2. Get r (and potentially p) array(s)        
            r = self._grid['node']['flow__receiver_node']
            
            # step 2. Stack, D, delta construction
            nd = flow_accum_bw._make_number_of_donors_array(r)
            delta = flow_accum_bw._make_delta_array(nd)
            D = flow_accum_bw._make_array_of_donors(r, delta)
            s = flow_accum_bw.make_ordered_node_array(r, self.fd.sink)
            
            # put theese in grid so that depression finder can use it.         
            # store the generated data in the grid
            self._grid['node']['flow__data_structure_delta'][:] = delta[1:]
            self._grid['link']['flow__data_structure_D'][:len(D)] = D
            self._grid['node']['flow__upstream_node_order'][:] = s
            
            # step 3. Run depression finder if passed 
            # at present this must go at the end. 
            
            # step 4. Accumulate (to one or to N depending on direction method. )
            a, q = flow_accum_bw.find_drainage_area_and_discharge(s, 
                                                                  r, 
                                                                  self.node_cell_area,
                                                                  self._grid.at_node['water__unit_flux_in'])

            self._grid['node']['drainage_area'][:] = a
            self._grid['node']['surface_water__discharge'][:] = q

        # at the moment, this is where the depression finder needs to live. 
        if self.depression_finder:
            self.df.map_depressions()

            
if __name__ == '__main__':
    import doctest
    doctest.testmod()    