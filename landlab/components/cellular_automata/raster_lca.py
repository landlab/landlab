#! /usr/env/python
"""
raster_lca.py: simple raster Landlab cellular automaton

This file defines the RasterLCA class, which is a sub-class of 
LandlabCellularAutomaton that implements a simple, non-oriented, raster-grid
CA. Like its parent class, RasterLCA implements a continuous-time, stochastic,
pair-based CA.

Created GT Sep 2014, starting from link_ca.py.
"""

from landlab_ca import LandlabCellularAutomaton, Transition
import landlab


class RasterLCA(LandlabCellularAutomaton):
    
    def __init__(self, model_grid, node_state_dict, transition_list,
                 initial_node_states):
        
        print 'RasterLCA.__init__ here'

        # Make sure caller has sent the right grid type        
        assert (type(model_grid) is landlab.grid.raster.RasterModelGrid), \
               'model_grid must be a Landlab RasterModelGrid'
               
        # Define the number of distinct cell-pair orientations: here just 1,
        # because RasterLCA represents a non-oriented CA model.
        self.number_of_orientations = 1
        
        # Call the LandlabCellularAutomaton.__init__() method to do the rest of
        # the initialization
        super(RasterLCA, self).__init__(model_grid, node_state_dict, 
            transition_list, initial_node_states)
        

if __name__=='__main__':
    print issubclass(RasterLCA, object)
    mg = landlab.RasterModelGrid(3, 4, 1.0)
    nsd = {0 : 'yes', 1 : 'no'}
    xnlist = []
    xnlist.append( Transition( (0,1,0), (1,1,0), 1.0, 'frogging' ) )
    nsg = mg.add_zeros('node', 'node_state_grid')
    rlca = RasterLCA(mg, nsd, xnlist, nsg)
    print rlca.__dict__
