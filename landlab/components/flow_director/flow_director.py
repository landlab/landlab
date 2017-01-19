#! /usr/env/python

"""
Short title.

Description text.
"""

from __future__ import print_function
from landlab import Component
from landlab import RasterModelGrid  # for type tests
from landlab.utils.decorators import use_field_name_or_array

@use_field_name_or_array('node')
def return_surface(grid, surface):
    """
    """
    return surface


class _FlowDirector(Component):
    """
    Private class for creating components to calculate flow directions.

    This class is not meant to be used directly in modeling efforts.
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
    >>> from landlab.components.flow_director.flow_director import(
    ... _FlowDirector)
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at = 'node')
    >>> fd = _FlowDirector(mg, 'topographic__elevation')
    >>> fd.elevs
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> list(mg.at_node.keys())
    ['topographic__elevation']
    """

    _name = '_FlowDirector'

    def __init__(self, grid, surface):
        # We keep a local reference to the grid
        super(_FlowDirector, self).__init__(grid)

        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self.method = None

        # test input variables are present:
        self.surface = surface
        surf = return_surface(grid, surface)

        # add elevations as a local variable.
        self.elevs = surf

    def run_one_step(self):

        """
        Short Description.

        Long description.
        """
        raise NotImplementedError('run_one_step()')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
