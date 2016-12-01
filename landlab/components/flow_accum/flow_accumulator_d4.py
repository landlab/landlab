from landlab.components.flow_accum.flow_accumulator import FlowAccumulator 
from landlab.components.flow_director import FlowDirectorD4 as FlowDirector
from landlab.components.flow_accum import flow_accum_bw 
from landlab import VoronoiDelaunayGrid

class FlowAccumulatorD4(FlowAccumulator):
    """Single-path (steepest direction) flow routing by the D4 method.

    This class implements single-path (steepest direction) flow routing, and
    calculates flow directions, drainage area, and discharge.

    Note that for Voronoi-based grids there is no difference between
    D4 and D8 methods. For this reason the D4 method is not implemented for
    Voroni grids. Use FlowAccumulatorD8  instead. 
    
    The perimeter nodes  NEVER contribute to the accumulating flux, even if the 
    gradients from them point inwards to the main body of the grid. This is 
    because under Landlab definitions, perimeter nodes lack cells, so cannot 
    accumulate any discharge.

    This is accomplished by first finding D4 flow directions using the 
    component FlowDirectorD4, and then calculating the accumulation area and
    discharge. 
    
    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished 
    together. 
    
    Stores as ModelGrid fields:
               
        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure: 
            *flow__data_structure_delta*. The first element is always zero.
        -  Link array of the D data structure: *flow__data_structure_D*
        
    The FlowDirectorD4 component adds the additional ModelGrid fields:
        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Node array containing ID of link that leads from each node to its
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.   
    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.  
    depression_finder : Component, optional
        A depression finding component
        
    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulatorD4
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
    >>> fa=FlowAccumulator(mg, 'topographic__elevation')
    >>> fa.elevs
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> mg_2 = RasterModelGrid((5, 4), spacing=(1, 1))
    >>> elev = np.array([0.,  0.,  0., 0.,
    ...                  0., 21., 10., 0.,
    ...                  0., 31., 20., 0.,
    ...                  0., 32., 30., 0.,
    ...                  0.,  0.,  0., 0.])
    >>> _ = mg_2.add_field('node','topographic__elevation', elev)
    >>> mg_2.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa_2 = FlowAccumulatorD4(mg_2)
    >>> fa_2.run_one_step()
    >>> mg_2.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,  3,  
            4,  1,  2,  7,  
            8, 10,  6, 11,
           12, 14, 10, 15, 
           16, 17, 18, 19])

    """
    
    _name = 'FlowAccumulatorD4'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, surface='topographic__elevation', depression_finder=None):
        super(FlowAccumulatorD4, self).__init__(grid, surface)

        # save method as attribute
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        self.method = 'D4'
        if self._is_Voroni:
            raise NotImplementedError('FlowAccumulatorD4 not implemented for irregular grids, use FlowAccumulatorD8')
        # save 
        self.fd = FlowDirector(self._grid, self.elevs)
        self.df = depression_finder
        
    def run_one_step(self):
        
        # step 0. Check and update BCs
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code        
        
        # step 1. Find flow directions by specified method
        self.fd.run_one_step()
        
        # step 2. Get r (and potentially p) array(s)        
        r = self._grid['node']['flow__receiver_node']
        s = self.fd.baselevel_nodes       
        
        # step 2. Stack, D, delta construction
        nd = flow_accum_bw._make_number_of_donors_array(r)
        delta = flow_accum_bw._make_delta_array(nd)
        D = flow_accum_bw._make_array_of_donors(r, delta)
        s = flow_accum_bw.make_ordered_node_array(r, s)
        
        #put theese in grid so that depression finder can use it.         
        # store the generated data in the grid
        self._grid['node']['flow__data_structure_delta'][:] = delta[1:]
        self._grid['link']['flow__data_structure_D'][:len(D)] = D
        self._grid['node']['flow__upstream_node_order'][:] = s
        
        # step 3. Initialize and Run depression finder if passed 
        if self.df:
            df=self.df(self.grid)
            df.run_one_step()
        
        # step 4. Accumulate (to one or to N depending on direction method. )
        a, q = flow_accum_bw.find_drainage_area_and_discharge(self._grid['node']['flow__upstream_node_order'], 
                                                              self._grid['node']['flow__receiver_node'], 
                                                              self.node_cell_area,
                                                              self._grid.at_node['water__unit_flux_in'])
        
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['surface_water__discharge'][:] = q
        

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()    