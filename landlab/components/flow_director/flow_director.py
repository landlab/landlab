#! /usr/env/python

"""flow_director.py provides a private class to help create FlowDirectors.

Provides the _FlowDirector component which does grid type testing, adds
the surface over which flow will be routed to the component, and sets up
part of the boundary condition testing.
"""

from landlab import Component
from landlab import RasterModelGrid  # for type tests
from landlab.utils.return_array import return_array_at_node


class _FlowDirector(Component):
    """Private class for creating components to calculate flow directions.

    This class is not meant to be used directly in modeling efforts.
    Instead it has the functionality that all flow direction calculators need
    to initialize and check boundary conditions.

    It also creates the following field used by all FlowDirectors.

    -  Link array identifing if flow goes with (1) or against (-1) the link
       direction: *'flow__link_direction'*

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
    >>> from landlab.components.flow_director.flow_director import _FlowDirector
    >>> mg = RasterModelGrid((3, 3), xy_spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field(
    ...     "topographic__elevation",
    ...     mg.node_x + mg.node_y,
    ...     at="node",
    ... )
    >>> fd = _FlowDirector(mg, "topographic__elevation")
    >>> fd.surface_values
    array([0., 1., 2., 1., 2., 3., 2., 3., 4.])
    >>> "topographic__elevation" in mg.at_node.keys()
    True
    >>> "flow__sink_flag" in mg.at_node.keys()
    True

    _FlowDirector also works if you pass it an array instead of a field name.

    >>> import numpy as np
    >>> mg = RasterModelGrid((3, 3), xy_spacing=(1, 1))
    >>> z = np.array([0.0, 1.0, 2.0, 1.0, 2.0, 3.0, 2.0, 3.0, 4.0])
    >>> fd = _FlowDirector(mg, z)
    >>> fd.surface_values
    array([0., 1., 2., 1., 2., 3., 2., 3., 4.])
    """

    _name = "_FlowDirector"

    _unit_agnostic = True

    _info = {
        "flow__sink_flag": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Boolean array, True at local lows",
        }
    }

    def __init__(self, grid, surface):
        """Initialize the _FlowDirector class."""
        # We keep a local reference to the grid
        super().__init__(grid)

        self._bc_set_code = self._grid.bc_set_code

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self._method = None

        # save elevations as class properites.
        self._surface = surface
        self._surface_values = return_array_at_node(grid, surface)

        grid.add_zeros("flow__sink_flag", at="node", dtype=bool, clobber=True)

    @property
    def surface_values(self):
        """Values of the surface over which flow is directed."""
        return self._surface_values

    def _changed_surface(self):
        """Check if the surface values have changed.

        If the surface values are stored as a field, it is important to
        check if they have changed since the component was instantiated.
        """
        if isinstance(self._surface, str):
            self._surface_values = return_array_at_node(self._grid, self._surface)

    def _check_updated_bc(self):
        # step 0. Check and update BCs
        if self._bc_set_code != self._grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self._grid.bc_set_code

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError("run_one_step()")

    @property
    def sink_flag(self):
        """Return the array with sink flags."""
        return self._grid["node"]["flow__sink_flag"]

    @property
    def node_steepest_slope(self):
        """Return the steepest link slope at a node."""
        return self._grid["node"]["topographic__steepest_slope"]

    @property
    def link_to_flow_receiving_node(self):
        """Return the link id along the link transporting flow."""
        return self._grid["node"]["flow__link_to_receiver_node"]

    @property
    def node_receiving_flow(self):
        """Return the node ids of the nodes receiving flow."""
        return self._grid["node"]["flow__receiver_node"]


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
