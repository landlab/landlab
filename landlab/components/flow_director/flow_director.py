#! /usr/env/python
from __future__ import print_function


#from landlab.components.flow_routing import flow_direction_DN
from landlab import FieldError, Component
from landlab import RasterModelGrid
from landlab.utils.decorators import use_file_name_or_kwds

class FlowDirector(Component):

    """
    """

    _name = 'FlowDirector'

    @use_file_name_or_kwds
    def __init__(self, grid, surface='topographic_elevation'):
        # We keep a local reference to the grid
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        
        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self.method = None

        self.updated_boundary_conditions()

        
        # test input variables are present:
        grid.at_node[surface]
        
        # Keep track of the following variables:
        #   - flow__sink_flag
        #  Note: all other output fields are method specific.  
        try:
            self.sink_flag = grid.add_zeros('flow__sink_flag', at='node',
                                                dtype=float)
        except FieldError:
            self.sink_flag = grid.at_node['flow__sink_flag']   

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
    
    def run_one_step(self):
        print ('You have called run_one_step() on the base FlowDirection class. ' +
               'This component does not actually direct flow, but instead ' +
               'sets up all of the core functionallity of the FlowDirection components. ' +
               'You probably want to initiate a flow direction finder with a name like:\n'+
               'FlowDirection*method*\n such as FlowDirection_D4 or ' +
               'FlowDirection_D8')

        
if __name__ == '__main__':
    import doctest
    doctest.testmod()