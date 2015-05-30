#! /usr/env/python
"""
raster_lca.py: simple raster Landlab cellular automaton

This file defines the RasterLCA class, which is a sub-class of 
LandlabCellularAutomaton that implements a simple, non-oriented, raster-grid
CA. Like its parent class, RasterLCA implements a continuous-time, stochastic,
pair-based CA.

Created GT Sep 2014, starting from link_ca.py.
"""

from .landlab_ca import LandlabCellularAutomaton, Transition
import landlab


class RasterLCA(LandlabCellularAutomaton):
    """
    Class RasterLCA implements a non-oriented raster CellLab-CTS model.
    """
    def __init__(self, model_grid, node_state_dict, transition_list,
                 initial_node_states, prop_data=None, prop_reset_value=None):
        """
        RasterLCA constructor: sets number of orientations to 1 and calls
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
        prop_reset_value : (scalar; same type as entries in prop_data) (optional)
            Default or initial value for a node/cell property (e.g., 0.0)
        """
        print 'WARNING: use of RasterLCA is deprecated.'
        print 'Use RasterCTS instead.'
        
        # Make sure caller has sent the right grid type        
        assert (type(model_grid) is landlab.grid.raster.RasterModelGrid), \
               'model_grid must be a Landlab RasterModelGrid'
               
        # Define the number of distinct cell-pair orientations: here just 1,
        # because RasterLCA represents a non-oriented CA model.
        self.number_of_orientations = 1
        
        # Call the LandlabCellularAutomaton.__init__() method to do the rest of
        # the initialization
        super(RasterLCA, self).__init__(model_grid, node_state_dict, 
            transition_list, initial_node_states, prop_data, prop_reset_value)
        

if __name__=='__main__':
    import doctest
    doctest.testmod()
