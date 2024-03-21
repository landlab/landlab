#! /usr/env/python
"""Simple raster Landlab cellular automaton.

Simple raster Landlab cellular automaton, with
cell-pair transitions that depend on orientation (vertical or horizontal)

This file defines the OrientedRasterCTS class, which is a sub-class of
CellLabCTSModel that implements a simple, oriented, raster-grid
CA. Like its parent class, OrientedRasterCTS implements a continuous-time,
stochastic, pair-based CA.

Created GT Sep 2014
"""


import numpy as np

from ..grid import RasterModelGrid
from .celllab_cts import CellLabCTSModel


class OrientedRasterCTS(CellLabCTSModel):
    """Oriented raster CellLab-CTS model.

    RasterCTS constructor: sets number of orientations to 2 and calls
    base-class constructor.

    Parameters
    ----------
    model_grid : Landlab ModelGrid object
        Reference to the model's grid
    node_state_dict : dict
        Keys are node-state codes, values are the names associated with
        these codes
    transition_list : list of Transition objects
        List of all possible transitions in the model
    initial_node_states : array of ints (x number of nodes in grid)
        Starting values for node-state grid
    prop_data : array (x number of nodes in grid) (optional)
        Array of properties associated with each node/cell
    prop_reset_value : number or object, optional
        Default or initial value for a node/cell property (e.g., 0.0).
        Must be same type as *prop_data*.
    seed : int (default 0)
        Seed for random number generator

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.ca.celllab_cts import Transition
    >>> from landlab.ca.oriented_raster_cts import OrientedRasterCTS

    >>> mg = RasterModelGrid((3, 4))
    >>> nsd = {0: "yes", 1: "no"}
    >>> xnlist = []
    >>> xnlist.append(Transition((0, 1, 0), (1, 1, 0), 1.0, "frogging"))
    >>> nsg = mg.add_zeros("node", "node_state_grid")
    >>> orcts = OrientedRasterCTS(mg, nsd, xnlist, nsg)
    """

    def __init__(
        self,
        model_grid,
        node_state_dict,
        transition_list,
        initial_node_states,
        prop_data=None,
        prop_reset_value=None,
        seed=0,
    ):
        """RasterCTS constructor: sets number of orientations to 2 and calls
        base-class constructor.

        Parameters
        ----------
        model_grid : Landlab ModelGrid object
            Reference to the model's grid
        node_state_dict : dict
            Keys are node-state codes, values are the names associated with
            these codes
        transition_list : list of Transition objects
            List of all possible transitions in the model
        initial_node_states : array of ints (x number of nodes in grid)
            Starting values for node-state grid
        prop_data : array (x number of nodes in grid) (optional)
            Array of properties associated with each node/cell
        prop_reset_value : number or object, optional
            Default or initial value for a node/cell property (e.g., 0.0).
            Must be same type as *prop_data*.
        """

        # Make sure caller has sent the right grid type
        if not isinstance(model_grid, RasterModelGrid):
            raise TypeError("model_grid must be a Landlab RasterModelGrid")

        # Define the number of distinct cell-pair orientations: here just 1,
        # because RasterLCA represents a non-oriented CA model.
        self.number_of_orientations = 2

        # Call the LandlabCellularAutomaton constructor to do the rest of
        # the initialization
        super().__init__(
            model_grid,
            node_state_dict,
            transition_list,
            initial_node_states,
            prop_data,
            prop_reset_value,
            seed,
        )

    def setup_array_of_orientation_codes(self):
        """Creates and configures an array that contain the orientation code
        for each active link (and corresponding cell pair).

        Notes
        -----
        **Creates**:

        * ``self.active_link_orientation``: 1D numpy array of ints
          Array of orientation codes for each cell pair (link)

        This overrides the method of the same name in landlab_ca.py.
        """
        # Create array for the orientation of each active link
        self.link_orientation = np.zeros(self.grid.number_of_links, dtype=np.int8)

        # Set its value according to the different in y coordinate between each
        # link's TO and FROM nodes (the numpy "astype" method turns the
        # resulting array into integer format)
        dy = (
            self.grid.node_y[self.grid.node_at_link_head]
            - self.grid.node_y[self.grid.node_at_link_tail]
        )
        self.link_orientation = dy.astype(np.int8)
