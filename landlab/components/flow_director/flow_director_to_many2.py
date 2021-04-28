#! /usr/env/python

"""flow_director_to_many2.py provides a private class to help create
FlowDirectors.

Provides the _FlowDirectorToMany2 component which makes sure all model
grid fields are set up correctly. Similar to _FlowDirectorToMany but uses an
alternative field naming convention that can preserve the fields used by the
single-direction components AS single (steepest) direction indicators, with
the same data structure size, while adding extra data structures to represent
multiple receivers.
"""
import numpy as np

from landlab.components.flow_director.flow_director import _FlowDirector


_OUTPUT_FIELDS_WITH_1D_ARRAY = [
    "flow__link_to_receiver_node",
    "flow__receiver_node",
    "topographic__steepest_slope",
]

_OUTPUT_FIELDS_WITH_2D_ARRAY = [
    "flow__links_to_receiver_node",
    "flow__receiver_nodes",
    "flow__receiver_proportions",
    "topographic__downhill_slopes",
]


class _FlowDirectorToMany2(_FlowDirector):

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
    >>> from landlab.components.flow_director.flow_director_to_many2 import(
    ... _FlowDirectorToMany2)
    >>> mg = RasterModelGrid((3,3), xy_spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field(
    ...     "topographic__elevation",
    ...     mg.node_x + mg.node_y,
    ...     at="node",
    ... )
    >>> fd = _FlowDirectorToMany2(mg, 'topographic__elevation')
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> sorted(list(mg.at_node.keys()))
    ['flow__link_to_receiver_node', 'flow__links_to_receiver_node',
     'flow__receiver_node', 'flow__receiver_nodes', 'flow__receiver_proportions',
     'flow__sink_flag', 'topographic__downhill_slopes', 'topographic__elevation',
     'topographic__steepest_slope']
    """

    _name = "FlowDirectorToMany2"

    _unit_agnostic = True

    _info = {
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__links_to_receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "IDs of all link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of one receiver (node that receives flow from current node)",
        },
        "flow__receiver_nodes": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of all receivers (node that receives flow from current node)",
        },
        "flow__receiver_proportions": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "flow__sink_flag": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Boolean array, True at local lows",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "topographic__downhill_slopes": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The *downhill* slopes from each node",
        },
    }

    _max_receivers = 2

    def __init__(self, grid, surface):
        """Initialize the _FlowDirectorToMany2 class."""
        # run init for the inherited class
        super().__init__(grid, surface)
        self._to_n_receivers = "many"

        # Create output fields with the right size. Those in the list of fields
        # with 2D arrays have size (number of nodes) x (max # receivers)
        for name in _OUTPUT_FIELDS_WITH_1D_ARRAY:
            if name not in grid["node"]:
                self.initialize_one_output_field(name, at="node")
        for name in _OUTPUT_FIELDS_WITH_2D_ARRAY:
            if name not in grid["node"]:
                self.initialize_one_output_field(
                    name, at="node", values_per_element=self._max_receivers
                )

        self._receivers = grid.at_node["flow__receiver_nodes"]
        if np.all(self._receivers == 0):
            self._receivers.fill(self._grid.BAD_INDEX)

        self._receiver_links = grid.at_node["flow__links_to_receiver_node"]
        if np.all(self._receiver_links == 0):
            self._receiver_links.fill(self._grid.BAD_INDEX)

        self._proportions = grid.at_node["flow__receiver_proportions"]
        self._steepest_slope = grid.at_node["topographic__downhill_slopes"]

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError("run_one_step()")

    def map_steepest_directions(self):
        """Identify and store the receivers, links, and slopes in the steepest
        direction.
        
        A many-direction component stores multiple receivers, receiver links,
        and downhill slopes for each node. This method finds, for each node,
        the receiver, link, and downhill gradient (slope) in the steepest among
        the multiple directions. The multiple receivers, links, and slopes are
        in the fields flow__receiver_nodes, flow__links_to_receiver_node, and
        topographic__downhill_slopes, respectively (note the plural). The
        *steepest* of these, for each node, is recorded in the fields
        flow__receiver_node, flow__link_to_receiver_node, and 
        topographic__steepest_slope, respectively (note the singular).
        """
        index_array = np.argmax(self._receivers, axis=-1)
        steepest_rcvr = self._grid.at_node["flow__receiver_node"]
        steepest_rcvr[:] = np.take_along_axis(
            self._receivers, np.expand_dims(index_array, axis=-1), axis=-1
        ).squeeze(axis=-1)
        link_to_steep_rcvr = self._grid.at_node["flow__link_to_receiver_node"]
        link_to_steep_rcvr[:] = np.take_along_axis(
            self._receiver_links, np.expand_dims(index_array, axis=-1), axis=-1
        ).squeeze(axis=-1)
        steepest_slope = self._grid.at_node["topographic__steepest_slope"]
        steepest_slope[:] = np.take_along_axis(
            self._grid.at_node["topographic__downhill_slopes"], np.expand_dims(index_array, axis=-1), axis=-1
        ).squeeze(axis=-1)

    @property
    def proportions_of_flow(self):
        """Return the proportions of flow going to recievers."""
        return self._grid["node"]["flow__receiver_proportions"]


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
