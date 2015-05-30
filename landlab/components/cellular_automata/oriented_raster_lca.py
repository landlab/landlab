#! /usr/env/python
"""
oriented_raster_lca.py: simple raster Landlab cellular automaton, with 
cell-pair transitions that depend on orientation (vertical or horizontal)

This file defines the OrientedRasterLCA class, which is a sub-class of 
LandlabCellularAutomaton that implements a simple, oriented, raster-grid
CA. Like its parent class, OrientedRasterLCA implements a continuous-time, 
stochastic, pair-based CA.

Created GT Sep 2014
"""
import warnings

from numpy import zeros
import six

from .landlab_ca import LandlabCellularAutomaton, Transition
import landlab

_DEBUG = False

class OrientedRasterLCA(LandlabCellularAutomaton):
    """
    Class OrientedRasterLCA implements an oriented raster CellLab-CTS model.
    """
    def __init__(self, model_grid, node_state_dict, transition_list,
                 initial_node_states, prop_data=None, prop_reset_value=None):
        """
        RasterLCA constructor: sets number of orientations to 2 and calls
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
        warnings.warn('use of OrientedRasterLCA is deprecated. '
                      'Use OrientedRasterCTS instead.')
                
        if _DEBUG:
            six.print_('OrientedRasterLCA.__init__ here')

        # Make sure caller has sent the right grid type        
        assert (type(model_grid) is landlab.grid.raster.RasterModelGrid), \
               'model_grid must be a Landlab RasterModelGrid'
               
        # Define the number of distinct cell-pair orientations: here just 1,
        # because RasterLCA represents a non-oriented CA model.
        self.number_of_orientations = 2
        
        # Call the LandlabCellularAutomaton constructor to do the rest of
        # the initialization
        super(OrientedRasterLCA, self).__init__(model_grid, node_state_dict, 
            transition_list, initial_node_states, prop_data, prop_reset_value)
            
        if _DEBUG:
            six.print_('ORLCA:')
            six.print_(self.n_xn)
            six.print_(self.xn_to)
            six.print_(self.xn_rate)
        

    def setup_array_of_orientation_codes(self):
        """
        Creates and configures an array that contain the orientation code for 
        each active link (and corresponding cell pair).
        
        Parameters
        ----------
        (none)
        
        Returns
        -------
        (none)
        
        Creates
        -------
        self.active_link_orientation : 1D numpy array of ints
            Array of orientation codes for each cell pair (link)
        
        Notes
        -----
        This overrides the method of the same name in landlab_ca.py.
        """
        # Create array for the orientation of each active link
        self.active_link_orientation = zeros(self.grid.number_of_active_links, dtype=int)
    
        # Set its value according to the different in y coordinate between each
        # link's TO and FROM nodes (the numpy "astype" method turns the
        # resulting array into integer format)
        dy = self.grid.node_y[self.grid.link_tonode[self.grid.active_links]] \
             - self.grid.node_y[self.grid.link_fromnode[self.grid.active_links]]
        self.active_link_orientation = dy.astype(int)
        
        if _DEBUG:
            six.print_(self.active_link_orientation)
            
            
if __name__=='__main__':
    import doctest
    doctest.testmod()
