#! /usr/env/python
"""
    A python flow accumulation module. It is designed to be general, and to operate across multiple grids and multiple flow direction patterns. However, at the moment, only a steepest descent (single path) routing scheme is implemented.
    
    There remain some outstanding issues with the handling of boundary cells, which this component has inherited from flow_routing_D8.
    
    Created DEJH, 8/2013
"""
try:
    from scipy import weave
except ImportError:
    import warnings
    warnings.warn('unable to import scipy.weave')

import numpy as np
#weave.test()
#import flow_routing_D8

class AccumFlow(object):
    """
    This class allows the routing of flow around a landscape according to a previously calculated flow direction vector. It is not sensitive to grid type. It will eventually be able to work with discharges which are split across more than one node, but at the moment, assumes a single line of descent for a given node.
    """
    def __init__(self, grid, data):
        self.initialize(grid, data)

    def initialize(self, grid, data):
        self.flow_accum_by_area = np.zeros(grid.number_of_nodes+1) #prefilled with zeros, size of WHOLE grid+1, to allow -1 ids
        #Test if this module is to produce the flowacc data. It should be!
        try:
            data.flowacc
        except:
            data.flowacc = self.flow_accum_by_area[:-1]
        else:
            print "Another module has created the flow accumulation grid. Undesirable conflicts may occur!"
            assert(len(data.flowacc) == len(self.flow_accum_by_area[:-1]))

    def calc_flowacc(self, grid, data):
        active_cell_ids = grid.get_active_cell_node_ids()
        #Perform test to see if the flowdir data is a single vector, or multidimensional, here. Several ways possible: 1. Is the vector multidimensional?, e.g., try: data.flowdirs.shape[1] 2. set a flag in flowdir.
        
        try:
            height_order_active_cells = np.argsort(data.elev[active_cell_ids])[::-1] #descending order
        except:
            print 'Cells could not be sorted by elevation. Does the data object contain the elevation vector?'

        try:
            sorted_flowdirs = (data.flowdirs[active_cell_ids])[height_order_active_cells]
        except:
            print 'Flow directions could not be sorted by elevation. Does the data object contain the flow direction vector?'
        #print grid.cell_areas
        self.flow_accum_by_area[active_cell_ids] = grid.cell_areas #This is only the active nodes == cells by definition

        #print len(height_order_active_cells), len(sorted_flowdirs), len(self.flow_accum_by_area)
        #print height_order_active_cells
        #print sorted_flowdirs
        #print data.flowdirs
        #print self.flow_accum_by_area.reshape(5,5)

#---
# Two ways of routing flow are provided. All route flow in descending height order.
#The first, using weave, is not working due to an installation-dependent with the C++ compiler weave uses. However, it will be a massive improvement over other methods
#The second is an inefficient but functional looped method.

        cpp_code_fragment = """
printf ('Test');
"""
        #Shouldn't need to return_val, as we're handling mutable objects not ints

#        flow_accum_by_area = self.flow_accum_by_area
#        a=1.
#        weave.inline(cpp_code_fragment, ['a'], compiler='gcc') #['flow_accum_by_area', 'height_order_active_cells', 'sorted_flowdirs', 'active_cell_ids']) #,verbose=2, compiler='gcc')

#---
        ##inefficient Python code to mimic the above weave:
        for i in xrange(len(sorted_flowdirs)):
            iter_height_order = height_order_active_cells[i]
            iter_sorted_fldirs = sorted_flowdirs[i]
            self.flow_accum_by_area[iter_sorted_fldirs] += (self.flow_accum_by_area[active_cell_ids])[iter_height_order]                     
                                                            
        return self.flow_accum_by_area[:-1]

        #int downhill_node;
        #int active_node;
        #PyObject *active_node_array[(Nactive_cell_ids[0])];
        #for (int i=0; i<(Nactive_cell_ids[0]); i++)
        #    {
        #        active_node_array[i] = &(flow_accum_by_area[(active_cell_ids[i])]);
        #    }
        #
        #for (int i=0; i<(Nsorted_flowdirs[0]); i++)
        #{
        #    downhill_node = sorted_flowdirs[i];
        #    active_node = height_order_active_cells[i];
        #    if downhill_node != active_cell_ids[active_node]
        #    {
        #        flow_accum_by_area[downhill_node] += *(active_node_array[active_node]);
        #    }
        #}
