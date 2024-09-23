#! /usr/env/python
"""Simple hexagonal Landlab cellular automaton.

This file defines the HexCTS class, which is a sub-class of
CellLabCTSModel that implements a simple, non-oriented, hex-grid
CA. Like its parent class, HexCTS implements a continuous-time, stochastic,
pair-based CA. The hex grid has 3 principal directions, rather than 2 for a
raster. Hex grids are often used in CA models because of their symmetry.

Created GT Sep 2014
"""

from ..grid import HexModelGrid
from .celllab_cts import CellLabCTSModel


class HexCTS(CellLabCTSModel):
    """Class HexCTS implements a non-oriented hex-grid CellLab-CTS model.

    HexCTS constructor: sets number of orientations to 1 and calls
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

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> from landlab.ca.celllab_cts import Transition
    >>> from landlab.ca.hex_cts import HexCTS

    >>> mg = HexModelGrid((4, 3), spacing=1.0)
    >>> nsd = {0: "yes", 1: "no"}
    >>> xnlist = []
    >>> xnlist.append(Transition((0, 1, 0), (1, 1, 0), 1.0, "frogging"))
    >>> nsg = mg.add_zeros("node", "node_state_grid")
    >>> hcts = HexCTS(mg, nsd, xnlist, nsg)
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
        """HexCTS constructor: sets number of orientations to 1 and calls base-
        class constructor.

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
        """

        # Make sure caller has sent the right grid type
        if not isinstance(model_grid, HexModelGrid):
            raise TypeError("model_grid must be a Landlab HexModelGrid")

        # Define the number of distinct cell-pair orientations: here just 1,
        # because HexCTS represents a non-oriented CA model.
        self.number_of_orientations = 1

        # Call the LandlabCellularAutomaton.__init__() method to do the rest of
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
