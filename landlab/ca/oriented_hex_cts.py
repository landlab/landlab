#! /usr/env/python
"""Simple hexagonal Landlab cellular automaton

This file defines the OrientedHexCTS class, which is a sub-class of
CellLabCTSModel that implements a simple, non-oriented, hex-grid
CA. Like its parent class, OrientedHexCTS implements a continuous-time,
stochastic, pair-based CA. The hex grid has 3 principal directions, rather
than 2 for a raster. Hex grids are often used in CA models because of their
symmetry.

Created GT Sep 2014
"""
import numpy as np

from ..grid import HexModelGrid
from .celllab_cts import CellLabCTSModel


class OrientedHexCTS(CellLabCTSModel):

    """Oriented hex-grid CellLab-CTS model.

    OrientedHexCTS constructor: sets number of orientations to 3 and calls
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
    prop_data : array (x number of nodes in grid), optional
        Array of properties associated with each node/cell
    prop_reset_value : number or object, optional
        Default or initial value for a node/cell property (e.g., 0.0).
        Must be same type as *prop_data*.
    seed : int (default 0)
        Seed for random number generator

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> from landlab.ca.oriented_hex_cts import OrientedHexCTS
    >>> from landlab.ca.celllab_cts import Transition

    >>> mg = HexModelGrid(4, 3, 1.0)
    >>> nsd = {0 : 'yes', 1 : 'no'}
    >>> xnlist = []
    >>> xnlist.append(Transition((0,1,0), (1,1,0), 1.0, 'frogging'))
    >>> nsg = mg.add_zeros('node', 'node_state_grid')
    >>> ohcts = OrientedHexCTS(mg, nsd, xnlist, nsg)
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
        """Initialize a OrientedHexCTS.

        OrientedHexCTS constructor: sets number of orientations to 3 and calls
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
        prop_data : array (x number of nodes in grid), optional
            Array of properties associated with each node/cell
        prop_reset_value : number or object, optional
            Default or initial value for a node/cell property (e.g., 0.0).
            Must be same type as *prop_data*.
        """

        # Make sure caller has sent the right grid type
        if not isinstance(model_grid, HexModelGrid):
            raise TypeError("model_grid must be a Landlab HexModelGrid")

        # Define the number of distinct cell-pair orientations: here 3,
        # representing
        self.number_of_orientations = 3

        # Call the LandlabCellularAutomaton.__init__() method to do the rest of
        # the initialization
        super(OrientedHexCTS, self).__init__(
            model_grid,
            node_state_dict,
            transition_list,
            initial_node_states,
            prop_data,
            prop_reset_value,
            seed,
        )

    def setup_array_of_orientation_codes(self):
        """
        Creates and configures an array that contain the orientation code for
        each active link (and corresponding cell pair).

        Notes
        -----
        **Creates**:

        * ``self.active_link_orientation``: 1D numpy array

        This overrides the method of the same name in celllab_cts.py. If the
        hex grid is oriented such that one of the 3 axes is vertical (a
        'vertical' grid), then the three orientations are:

        * 0 = vertical (0 degrees clockwise from vertical)
        * 1 = right and up (60 degrees clockwise from vertical)
        * 2 = right and down (120 degrees clockwise from vertical)

        If the grid is oriented with one principal axis horizontal
        ('horizontal' grid), then the orientations are:

        * 0 = up and left (30 degrees counter-clockwise from vertical)
        * 1 = up and right (30 degrees clockwise from vertical)
        * 2 = horizontal (90 degrees clockwise from vertical)
        """
        self.link_orientation = np.zeros(self.grid.number_of_links, dtype=np.int8)
        for i in range(self.grid.number_of_links):
            dy = (
                self.grid.node_y[self.grid.node_at_link_head[i]]
                - self.grid.node_y[self.grid.node_at_link_tail[i]]
            )
            dx = (
                self.grid.node_x[self.grid.node_at_link_head[i]]
                - self.grid.node_x[self.grid.node_at_link_tail[i]]
            )
            if dx <= 0.0:
                self.link_orientation[i] = 0
            elif dy <= 0.0:
                self.link_orientation[i] = 2
            elif dx > 0.0 and dy > 0.0:
                self.link_orientation[i] = 1
