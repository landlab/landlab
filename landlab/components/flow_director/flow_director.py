#! /usr/env/python

"""
flow_director.py provides a private class to help create FlowDirectors.

Provides the _FlowDirector component which does grid type testing, adds the
surface over which flow will be routed to the component, and sets up part of
the boundary condition testing.
"""

from __future__ import print_function
from landlab import Component
from landlab import RasterModelGrid  # for type tests
from landlab.utils.return_array import return_array_at_node
import six


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
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> list(mg.at_node.keys())
    ['topographic__elevation']

    _FlowDirector also works if you pass it an array instead of a field name.
    >>> import numpy as np
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> z = np.array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> fd = _FlowDirector(mg, z)
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    """

    _name = '_FlowDirector'

    def __init__(self, grid, surface):
        """Initialize the _FlowDirector class."""
        # We keep a local reference to the grid
        super(_FlowDirector, self).__init__(grid)

        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self.method = None

        # save elevations as class properites.
        self.surface = surface
        self.surface_values = return_array_at_node(grid, surface)

    def _changed_surface(self):
        """Check if the surface values have changed.

        If the surface values are stored as a field, it is important to check
        if they have changed since the component was instantiated.
        """
        if isinstance(self.surface, six.string_types):
            self.surface_values = return_array_at_node(self._grid, self.surface)

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError('run_one_step()')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
