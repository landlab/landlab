#! /usr/env/python
from __future__ import print_function


#from landlab.components.flow_routing import flow_direction_DN
import landlab
import warnings

from landlab import FieldError, Component
from landlab import ModelParameterDictionary
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.decorators import use_field_name_or_array
import numpy


class FlowDirector(Component):

    """
    Base class for calculating flow directions. 
    
    This component is not meant to be used directly in modeling efforts. 
    Instead it has the functionality that all flow direction calculators need
    
    

    The primary method of this class, :func:`run_one_step` is not implemented.


    Parameters
    ----------
    grid : ModelGrid
        A grid.
    method : {'D8', 'D4'}, optional
        Routing method ('D8' is the default). This keyword has no effect for a
        Voronoi-based grid.
    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.
    """
    """

    _name = 'FlowDirector'

    @use_field_name_or_array('node')
    def __init__(self, grid, surface):
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
        
        # add elevations as a local variable.
        self.elevs = surface        

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
        raise NotImplementedError('run_one_step()')

        
if __name__ == '__main__':
    import doctest
    doctest.testmod()