#! /usr/env/python
"""
    A python flow accumulation module. It is designed to be general, and to operate across multiple grids and multiple flow direction patterns. However, at the moment, only a steepest descent (single path) routing scheme is implemented.
    
    There remain some outstanding issues with the handling of boundary cells, which this component has inherited from flow_routing_D8.
    
    Created DEJH, 8/2013
"""

from scipy import weave
import numpy as np
#import flow_routing_D8

class AccumFlow(object):
    """
    This class allows the routing of flow around a landscape according to a previously calculated flow direction vector. It is not sensitive to grid type. It will eventually be able to work with discharges which are split across more than one node, but at the moment, assumes a single line of descent for a given node.
    """
    def __init__(self, grid, data):
        self.initialize(grid, data)

    def initialize(self, grid, data):
        flow_accum_by_area = grid.create_active_cell_dvector()

    def calc_flowacc(self, grid, data):
        active_cell_ids = grid.get_active_cell_node_ids()
        try:
            height_order = np.argsort(data.elev[active_cell_ids])[::-1] #descending order
        except:
            print 'Cells could not be sorted by elevation. Does the data object contain the elevation vector?'

        try:
            sorted_flowdirs = data.flowdirs[height_order]
        except:
            print 'Flow directions could not be sorted by elevation. Does the data object contain the flow direction vector?'

        flow_accum_by_area = grid[active_cell_ids].cell_areas

        cpp_code_fragment = """
        for(int i=0; i<(Nsorted_flowdirs[0]-1); i++)
        {
            flow_accum_by_area[sorted_flowdirs[i]] += flow_accum_by_area[height_order[i]];
        }
        """
        #Note the -1 in the loop. Don't want to move flow out of the lowest cell. (Issues with BCells resurface here...?)
        #Shouldn't need to return_val, as we're handling mutable objects not ints

        weave.inline(cpp_code_fragment, ['flow_accum_by_area', 'height_order', 'sorted_flowdirs'])

        return flow_accum_by_area

#Treatment of active nodes very quickly becomes an issue! Will this work? flowdirs should be len(n_active_nodes), but will it be in the right order?

