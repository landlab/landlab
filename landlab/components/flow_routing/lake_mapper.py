# -*- coding: utf-8 -*-
"""
Created on Sat May 30 14:01:10 2015

@author: gtucker
"""
from __future__ import print_function

from landlab import ModelParameterDictionary, Component, FieldError, \
                    FIXED_VALUE_BOUNDARY
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.grid.base import BAD_INDEX_VALUE
import numpy  # Used for: count_nonzero
import landlab


# Codes for depression status
_UNFLOODED = 0
_PIT = 1
_CURRENT_LAKE = 2
_FLOODED = 3


class DepressionFinderAndRouter(Component):
    """
    This component identifies depressions in a topographic surface, finds an
    outlet for each depression, and [when GT or someone else finishes it]
    will modify the drainage directions accordingly.
    """
    _name = 'DepressionFinderAndRouter'
    
    _input_var_names = set(['topographic__elevation',
                            ])
    
    _output_var_names = set(['depression__depth',  # depth below spill point
                             'depression__outlet_node_ID',
                             ])
                             
    _var_units = {'depression__depth' : 'm',
                  'depression__outlet_node_ID' : '-'
                  }
    
    _var_mapping = {'depression__depth' : 'node',
                    'depression__outlet_node_ID' : 'node'
                    }
    
    _var_defs = {'topographic__elevation' : 'Surface topographic elevation',
                 'depression__depth' : 'Depth of depression below its spillway point'
                  }
    
    def __init__(self, grid, input_stream=None, current_time=0.):
        """
        Constructor assigns a copy of the grid, sets the current time, and
        calls the initialize method.
        """
        self._grid = grid
        self.current_time = current_time
        self.initialize(input_stream)


    def initialize(self, input_stream=None):
        """
        The BMI-style initialize method takes an optional input_stream
        parameter, which may be either a ModelParameterDictionary object or
        an input stream from which a ModelParameterDictionary can read values.
        """
        # Create a ModelParameterDictionary for the inputs
        if input_stream is None:
            inputs = None
        elif type(input_stream) == ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)
        
        # Make sure the grid includes elevation data. This means either:
        #  1. The grid has a node field called 'topographic__elevation', or
        #  2. The input file has an item called 'ELEVATION_FIELD_NAME' *and*
        #     a field by this name exists in the grid.
        try:
            self._elev = self._grid.at_node['topographic__elevation']
        except FieldError:
            try:
                topo_field_name = inputs.read_string('ELEVATION_FIELD_NAME')
            except AttributeError:
                print('Error: Because your grid does not have a node field')
                print('called "topographic__elevation", you need to pass the')
                print('name of a text input file or ModelParameterDictionary,')
                print(' and this file or dictionary needs to include the name')
                print(' of another field in your grid that contains your')
                print('elevation data.')
                raise AttributeError
            except MissingKeyError:
                print('Error: Because your grid does not have a node field')
                print('called "topographic__elevation", your input file (or')
                print('ModelParameterDictionary) must include an entry with')
                print('the key "ELEVATION_FIELD_NAME", which gives the name')
                print('of a field in your grid that contains your elevation')
                print('data.')
                raise MissingKeyError('ELEVATION_FIELD_NAME')
            try:
                self._elev = self._grid.at_node[topo_field_name]
            except AttributeError:
                print('Your grid does not seem to have a node field called', \
                      topo_field_name)
                
        # Create output variables.
        #
        # Note that we initialize depression depth to -1 (negative values make
        # no sense, so this is a clue to non-flooded nodes), and depression
        # outlet ID to BAD_INDEX_VALUE (which is a major clue!)
        self.depression_depth = self._grid.add_zeros('node', \
                                                     'depression__depth') - 1.0             
        self.depression_outlet = self._grid.add_zeros('node', \
                                                'depression__outlet_node_id', \
                                                dtype=int) + BAD_INDEX_VALUE               
                
        # Later on, we'll need a number that's guaranteed to be larger than the
        # highest elevation in the grid.
        self._BIG_ELEV = numpy.amax(self._elev) + 1
        
        # We'll also need a handy copy of the node neighbor lists
        # TODO: presently, this grid method seems to only exist for Raster grids.
        # We need it for *all* grids!
        self._node_nbrs = self._grid.get_neighbor_list()
        if type(self._grid) is landlab.grid.raster.RasterModelGrid:
            diag_nbrs = self._grid.get_diagonal_list()
            self._node_nbrs = numpy.concatenate((self._node_nbrs, diag_nbrs), 1)
            #print('NN',self._node_nbrs)
                
                
    def find_pits(self):
        """
        Locate local depressions ("pits") in a gridded elevation field.
        
        Uses
        ----
        self._elev, self._grid
        
        Creates
        -------
        self.is_pit : node array of booleans
            Flag indicating whether the node is a pit
        self.number_of_pits : int
            Number of pits found
        self.pit_node_ids : node array of ints
            IDs of the nodes that are pits
        
        Notes
        -----
        A node is defined as being a pit if and only if: 
            1. All neighboring core nodes have equal or greater elevation, and
            2. Any neighboring open boundary nodes have a greater elevation.
        
        The algorithm starts off assuming that all core nodes are pits. We then
        loop through all active links. For each link, if one node is higher
        than the other, the higher one cannot be a pit, so we flag it False.
        We also look at cases in which an active link's nodes have equal 
        elevations. If one is an open boundary, then the other must be a core
        node, and we declare the latter not to be a pit (via rule 2 above).
        """
        # Create the is_pit array, with all core nodes initialized to True and
        # all boundary nodes initialized to False.
        self.is_pit = self._grid.add_ones('node', 'is_pit', dtype=bool)
        self.is_pit[self._grid.boundary_nodes] = False
        
#        # Get a list of active links; in a raster, this (TODO: optionally)
#        # includes diagonals.
#        if type(self._grid) is landlab.grid.raster.RasterModelGrid:
#            (active_links, tails, heads) = self._grid.d8_active_links()
#        else:
#            active_links = self._grid.active_links
#            tails = self._grid.node_index_at_link_tail
#            heads = self._grid.node_index_at_link_head
#        print 'AL',active_links
#        print 'T',tails
#        print 'H',heads
        
        # Loop over all active links: if one of a link's two nodes is higher
        # than the other, the higher one is not a pit. Also, if they have
        # equal elevations and one is an open boundary, the other is not a pit.
        for link in self._grid.active_links:
            h = self._grid.node_at_link_head[link]
            t = self._grid.node_at_link_tail[link]
            if self._elev[h] > self._elev[t]:
                self.is_pit[h] = False
            elif self._elev[t] > self._elev[h]:
                self.is_pit[t] = False
            elif self._elev[h] == self._elev[t]:
                if self._grid.node_boundary_status[h] == FIXED_VALUE_BOUNDARY:
                    self.is_pit[t] = False
                elif self._grid.node_boundary_status[t] == FIXED_VALUE_BOUNDARY:
                    self.is_pit[h] = False
                    
        # If we have a raster grid, handle the diagonal active links too
        # (At the moment, their data structure is a bit different)
        # TODO: update the diagonal link data structures
        if type(self._grid) is landlab.grid.raster.RasterModelGrid:
            if not self._grid._diagonal_links_created:
                self._grid._setup_diagonal_links()
            for i in range(len(self._grid._diag_active_links)):
                h = self._grid._diag_activelink_tonode[i]
                t = self._grid._diag_activelink_fromnode[i]
                if self._elev[h] > self._elev[t]:
                    self.is_pit[h] = False
                elif self._elev[t] > self._elev[h]:
                    self.is_pit[t] = False
                elif self._elev[h] == self._elev[t]:
                    if self._grid.node_boundary_status[h] == FIXED_VALUE_BOUNDARY:
                        self.is_pit[t] = False
                    elif self._grid.node_boundary_status[t] == FIXED_VALUE_BOUNDARY:
                        self.is_pit[h] = False

        # Record the number of pits and the IDs of pit nodes.
        self.number_of_pits = numpy.count_nonzero(self.is_pit)
        (self.pit_node_ids, ) = numpy.where(self.is_pit)
        
        
    def find_lowest_node_on_lake_perimeter(self, nodes_this_depression):
        """
        Locate the lowest node on the margin of the 'lake'.
        """
        # Start with the first node on the list, and an arbitrarily large elev
        lowest_node = nodes_this_depression[0]
        lowest_elev = self._BIG_ELEV
        
        for n in nodes_this_depression:
            
            #print '  working on node',n,'neighbors are:'
            for nbr in self._node_nbrs[n]:
                #print nbr
                if nbr!=BAD_INDEX_VALUE:
                    if self.flood_status[nbr] == _UNFLOODED:
                        if self._elev[nbr] < lowest_elev:
                            lowest_node = nbr
                            lowest_elev = self._elev[nbr]
                    elif self.flood_status[nbr]==_PIT or \
                         self.flood_status[nbr] == _FLOODED:
                        nodes_this_depression.append(nbr)
                        self.flood_status[nbr] = _CURRENT_LAKE
                        
        # print '  lowest node on perim is',lowest_node,'with elev',lowest_elev
        return lowest_node
        
        
    def node_can_drain(self, the_node, nodes_this_depression):
        """
        Determine whether the given node has drainage away from the current
        lake/depression.
        """
        #print '    checking drainage for node',the_node
        for nbr in self._node_nbrs[the_node]:
            #print '    checking nbr',nbr
            if nbr != BAD_INDEX_VALUE:
                #print '    this nbr is valid'
                #print '    its flood status is',self.flood_status[nbr]
                if self._elev[nbr] < self._elev[the_node] and \
                   self.flood_status[nbr] != _CURRENT_LAKE and \
                   self.flood_status[nbr] != _FLOODED:  # caveat about outlet elevation...
                    #print '    Eureka!'
                    return True
        
        #print '    no, this node cannot drain'
        return False
            
        
    def is_valid_outlet(self, the_node, nodes_this_depression):
        """
        Determine whether the given node is a valid outlet for the depression.
        """
        if self._grid.node_boundary_status[the_node]==FIXED_VALUE_BOUNDARY:
            #print '   this node is an open boundary, so of course it can drain'
            return True
            
        if self.node_can_drain(the_node, nodes_this_depression):
            #print '   this node can drain!'
            return True
        #else:
            #print '   it cannot drain'
            
            
    def record_depression_depth_and_outlet(self, nodes_this_depression, outlet_id):
        """
        Record information about this depression/lake in the flood_status,
        depression_depth, and depression_outlet arrays.
        """
        for n in nodes_this_depression:
            self.flood_status[n] = _FLOODED
            self.depression_depth[n] = self._elev[outlet_id] - self._elev[n]
            self.depression_outlet[n] = outlet_id
            

    def find_depression_from_pit(self, pit_node):
        """
        Identify extent of depression/lake whose lowest point is the node
        pit_node (which is a itself a pit, a.k.a., closed depression).
        """
        #print ' processing pit at',pit_node
        
        # Place pit_node at top of depression list
        nodes_this_depression = []
        nodes_this_depression.insert(0, pit_node)
        
        # Flag the pit as being _CURRENT_LAKE (it's the first node in the
        # current lake)
        self.flood_status[pit_node] = _CURRENT_LAKE
        
        # This flag keeps track of when we're done with this depression
        found_outlet = False
        
        # Safety check
        count = 0
        max_count = self._grid.number_of_nodes + 1
        
        while not found_outlet:
            
            lowest_node_on_perimeter = self.find_lowest_node_on_lake_perimeter(nodes_this_depression)
            
            found_outlet = self.is_valid_outlet(lowest_node_on_perimeter, nodes_this_depression)
            
            if not found_outlet:
                
                # Add lowest_node to the lake list
                nodes_this_depression.append(lowest_node_on_perimeter)
                
                # Flag it as being part of the current lake/depression
                self.flood_status[lowest_node_on_perimeter] = _CURRENT_LAKE
                
                #print ' lake list:',nodes_this_depression
                
            #else:
                #print ' node',lowest_node_on_perimeter,'has drainage and'
                #print ' therefore it is the outlet point'                
            
            # Safety check, in case a bug (ha!) puts us in an infinite loop
            assert (count<max_count), 'too many iterations in lake filler!'
        
        # Now that we've mapped this depression, record it in the arrays
        # depression_depth, depression_outlet, and flood_status
        self.record_depression_depth_and_outlet(nodes_this_depression, lowest_node_on_perimeter)
        
        # TODO: ideally we need a way to keep track of the number, area extent,
        # and average depth of depressions. Tricky thing is that one might be
        # devoured by another, so would need to be removed from the list.
        
        
    def identify_depressions_and_outlets(self):
        """
        Find and map the depressions/lakes in a topographic surface,
        given a previously identified list of pits (if any) in the surface.
        """
        for pit_node in self.pit_node_ids:
            #print 'mapping lake starting from pit',pit_node
            self.find_depression_from_pit(pit_node)
      
                
    def map_depressions(self):
        """
        Map depressions/lakes in a topographic surface.        
        
        Examples
        --------
        Test #1: 5x5 raster grid with a diagonal lake.
        
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter

        >>> rg = RasterModelGrid(5, 5)
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = np.array([100.,100.,95.,100.,100.,100.,101.,92.,1.,100.,100.,101.,2.,101.,100.,100.,3.,101.,101.,100.,90.,95.,100.,100.,100])
        >>> df = DepressionFinderAndRouter(rg)
        >>> df.map_depressions()
        >>> df.display_depression_map()
        . . . . . 
        . . . ~ . 
        . . ~ . . 
        . ~ . . . 
        o . . . . 
        """
        # Locate nodes with pits
        self.find_pits()
        
        # Set up "lake code" array
        self.flood_status = self._grid.add_zeros('node', 'flood_status_code') \
                            + _UNFLOODED
        self.flood_status[self.pit_node_ids] = _PIT
        
        self.identify_depressions_and_outlets()
        
        
    def display_depression_map(self):
        """
        Display on screen a simple character-based map of depressions/lakes.
        """
        # Find the outlet nodes (just for display purposes)
        is_outlet = numpy.zeros(self._grid.number_of_nodes, dtype=bool)
        for i in self._grid.core_nodes:
            if self.flood_status[i]==_FLOODED:
                is_outlet[self.depression_outlet[i]] = True
        
        n=0
        for r in range(self._grid.number_of_node_rows):
            for c in range(self._grid.number_of_node_columns):
                if is_outlet[n]:
                    print('o', end=' ')
                elif self.flood_status[n]==_UNFLOODED:
                    print('.', end=' ')
                else:
                    print('~', end=' ')
                n+=1
            print()
        
        
        
    

def main():
    """
    temporary: test.
    """
    print('howdy')
    from landlab import RasterModelGrid
    from numpy.random import rand
    grid = RasterModelGrid(4, 5, 1.0)
    z = grid.add_zeros('node', 'topographic__elevation')
    z[:] = rand(len(z))*100
    print(z)
    dep_finder = DepressionFinderAndRouter(grid, '/Users/gtucker/Dev/Landlab/gt_tests/test_inputs_for_depression_mapper.txt')
    #dep_finder.initialize()
    dep_finder.map_depressions()
    
    n=0
    for r in range(grid.number_of_node_rows-1, -1, -1):
        for c in range(grid.number_of_node_columns):
            print(int(z[n]),'(',dep_finder.is_pit[n],')', end=' ')
            n+=1
        print()
        
    n=0
    for r in range(grid.number_of_node_rows):
        for c in range(grid.number_of_node_columns):
            if dep_finder.depression_outlet[n]==BAD_INDEX_VALUE:
                print(dep_finder.depression_depth[n],'( X )', end=' ')
            else:
                print(dep_finder.depression_depth[n],'(',dep_finder.depression_outlet[n],')', end=' ')
            n+=1
        print()

    dep_finder.display_depression_map()
    
    # Now, route flow through.
    # First, find flow dirs without lakes. Then, adjust.
    from landlab.components.flow_routing.flow_direction_DN import grid_flow_directions
    (rcvr, ss) = grid_flow_directions(grid, z)
    print(z)
    print(rcvr)
    print(ss)
                    
    
    
    
if __name__=='__main__':
    import doctest
    doctest.testmod()
    main()

    
