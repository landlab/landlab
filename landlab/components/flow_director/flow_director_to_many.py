#! /usr/env/python

"""flow_director_to_many.py provides a private class to help create
FlowDirectors.

Provides the _FlowDirectorToMany component which makes sure all model
grid fields are set up correctly.
"""
from landlab.components.flow_director.flow_director import _FlowDirector


class _FlowDirectorToMany(_FlowDirector):

    """Private class for creating components to calculate flow directions.

    This class is not meant to be used directly in modeling efforts. It
    inherits from the _FlowDirector class and builds on it to provide the
    functionality that all flow direction calculators need if they direct flow
    only to multiple nodes, as in  D infinity or MFD direction finding. It
    exists in contrast to the other intermediate flow director class
    _FlowDirectorToOne which provides equivalent functionality for flow
    direction algorithms such as D8 or steepest descent which directs flow only
    to one other node. As the primary difference between these two methods is
    the names of the fields they create and use, the primary function of this
    class is to create model grid fields.

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
    >>> from landlab.components.flow_director.flow_director_to_many import(
    ... _FlowDirectorToMany)
    >>> mg = RasterModelGrid((3,3), xy_spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field(
    ...     'topographic__elevation',
    ...     mg.node_x + mg.node_y,
    ...     at = 'node'
    ... )
    >>> fd = _FlowDirectorToMany(mg, 'topographic__elevation')
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> sorted(list(mg.at_node.keys()))
    ['flow__sink_flag', 'topographic__elevation']
    """

    _name = "FlowDirectorToMany"

    _input_var_names = ("topographic__elevation",)

    _output_var_names = ("flow__sink_flag",)

    _var_units = {"topographic__elevation": "m", "flow__sink_flag": "-"}

    _var_mapping = {"topographic__elevation": "node", "flow__sink_flag": "node"}

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "flow__sink_flag": "Boolean array, True at local lows",
    }

    def __init__(self, grid, surface):
        """Initialize the _FlowDirectorTo_One class."""
        # run init for the inherited class
        super(_FlowDirectorToMany, self).__init__(grid, surface)
        self.to_n_receivers = "many"
        # initialize new fields

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError("run_one_step()")

    @property
    def proportions_of_flow(self):
        """Return the proportions of flow going to recievers."""
        return self._grid["node"]["flow__receiver_proportions"]


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
