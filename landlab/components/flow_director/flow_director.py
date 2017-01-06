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


@use_field_name_or_array('node')
def return_surface(grid, surface):
    return(surface)

class FlowDirector(Component):

    """
    Base class for calculating flow directions. 
    
    This component is not meant to be used directly in modeling efforts. 
    Instead it has the functionality that all flow direction calculators need
    to initialize and check boundary conditions.
    
    The primary method of this class, :func:`run_one_step` is not implemented.


    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.   
        
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_director.flow_director import FlowDirector
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
    >>> fd=FlowDirector(mg, 'topographic__elevation')
    >>> fd.elevs
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> list(mg.at_node.keys())
    ['topographic__elevation']

    """

    _name = 'FlowDirector'
    
    def __init__(self, grid, surface):
        # We keep a local reference to the grid
        if hasattr(self, 'method') == False:
            self.method = 'base'
            
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        
        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self.method = None


        # test input variables are present:
        self.surface=surface
        surf=return_surface(grid, surface)
        
        # add elevations as a local variable.
        self.elevs = surf        
    
    def run_one_step(self):
        raise NotImplementedError('run_one_step()')

        
if __name__ == '__main__':
    import doctest
    doctest.testmod()