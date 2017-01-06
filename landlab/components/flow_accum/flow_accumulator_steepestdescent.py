from landlab.components.flow_accum.flow_accumulator import FlowAccumulator 
from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest as FlowDirector
from landlab.components.flow_accum import flow_accum_bw 
from landlab import VoronoiDelaunayGrid

class FlowAccumulatorSteepestDescent(FlowAccumulator):
    """
    Single-path (steepest direction) flow routing by the steepest descent 
    method for irregular grids
    
    For Raster Grids use FlowAccumulatorD4 or FlowAccumulatorD8. 

    This class implements single-path (steepest direction) flow routing, and
    calculates flow directions, drainage area, and discharge.
    
    The perimeter nodes  NEVER contribute to the accumulating flux, even if the 
    gradients from them point inwards to the main body of the grid. This is 
    because under Landlab definitions, perimeter nodes lack cells, so cannot 
    accumulate any discharge.

    This is accomplished by first finding steepest descent flow directions 
    using the component FlowDirectorSteepestDescent, and then calculating the 
    accumulation area and discharge. 
    
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
        A grid of type Voroni.
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
    >>> from landlab import HexModelGrid
    >>> from landlab.components import FlowAccumulatorSteepestDescent
    
    For the example we will use a Hexagonal Model Grid, a special type of 
    Voroni Grid that has regularly spaced hexagonal cells. 
    
    >>> mg = HexModelGrid(5,3)
    >>> _ = mg.add_field('topographic__elevation', mg.node_x + np.round(mg.node_y), at = 'node')
    >>> fa=FlowAccumulatorSteepestDescent(mg, 'topographic__elevation')
    >>> fa.elevs
    array([ 0. ,  1. ,  2. ,  
            0.5,  1.5,  2.5,  3.5,  
            1. ,  2. ,  3. ,  4. ,  5. ,  
            2.5,  3.5,  4.5,  5.5,  
            3. ,  4. ,  5. ])
    
    For the second exampl, we will also set the dx spacing such that each cell 
    has an area of one. 
    
    >>> dx=(2./(3.**0.5))**0.5
    >>> mg_2 = HexModelGrid(5,3, dx)
    >>> _ = mg_2.add_field('topographic__elevation', mg.node_x + np.round(mg.node_y), at = 'node')
    >>> fa_2 = FlowAccumulatorSteepestDescent(mg_2)
    >>> fa_2.run_one_step()
    >>> mg_2.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,  
            3,  0,  1,  6,  
            7,  3,  4,  5, 11, 
           12,  8,  9, 15, 
           16, 17, 18])
    >>> mg_2.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 3.,  2.,  0., 
            2.,  3.,  2.,  0.,  
            0.,  2.,  2.,  1.,  0.,  
            0., 1.,  1.,  0.,  
            0.,  0.,  0.])

    Now let's change the cell area (100.) and the runoff rates:

    >>> mg_3 = HexModelGrid(5,3, dx*100.)

    Put the data back into the new grid.

    >>> _ = mg_3.add_field('topographic__elevation', mg_3.node_x + np.round(mg_3.node_y), at = 'node')
    >>> fa_3 = FlowAccumulatorSteepestDescent(mg_3)
    >>> runoff_rate = np.arange(mg_3.number_of_nodes)
    >>> _ = mg_3.add_field('node', 'water__unit_flux_in', runoff_rate,
    ...                  noclobber=False)
    >>> fa_3.run_one_step()
    >>> mg_3.at_node['surface_water__discharge']
    array([ 270000.,  150000.,       0.,  
            210000.,  270000.,  150000.,  0.,
                 0.,  210000.,  230000.,  100000.,       0.,
                 0.,  130000.,  140000.,       0.,       
                 0.,       0.,       0.])
    
    Finally, see what happens when there is a depression. 

    >>> mg_4 = HexModelGrid(9,5, dx)
    >>> _ = mg_4.add_field('topographic__elevation', mg_4.node_x + np.round(mg_4.node_y), at = 'node')
    >>> depression_ids=[21,22,29,30,31, 38, 39]
    >>> mg_4.at_node['topographic__elevation'][depression_ids] *= 0.1

    This model grid has a depression in the center. 
    
    >>> mg_4.at_node['topographic__elevation']
    array([  0.        ,   1.07456993,   2.14913986,   3.2237098 ,   4.29827973,   
             0.46271503,   1.53728497,   2.6118549 ,   3.68642483,   4.76099476,   5.83556469,   
             0.92543007,   2.        ,   3.07456993,   4.14913986,   5.2237098 ,   6.29827973,   7.37284966,   
             1.3881451 ,   2.46271503,   3.53728497,   0.46118549,   0.56864248,   6.76099476,   7.83556469,   8.91013463,   
             1.85086014,   2.92543007,   4.        ,   0.50745699,   0.61491399,   0.72237098,   8.29827973,   9.37284966,  10.44741959,   
             3.3881451 ,   4.46271503,   5.53728497,   0.66118549,   0.76864248,   8.76099476,   9.83556469,  10.91013463,   
             4.92543007,   6.        ,   7.07456993,   8.14913986,   9.2237098 ,  10.29827973,  11.37284966,   
             6.46271503,   7.53728497,   8.6118549 ,   9.68642483,  10.76099476,  11.83556469,
             7.        ,   8.07456993,   9.14913986,  10.2237098 ,  11.29827973])
    >>> fa_4 = FlowAccumulatorSteepestDescent(mg_4) 
    >>> fa_4.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg_4.at_node['flow__receiver_node']
    array([ 0,   1,  2,  3,  4,  
            5,   0,  1,  2,  3, 10, 
            11,  5, 21, 21, 22,  9, 17, 
            18, 11, 21, 21, 21, 22, 16, 25, 
            26, 18, 29, 21, 21, 22, 31, 24, 34, 
            35, 27, 29, 29, 30, 31, 32, 42, 
            43, 36, 38, 38, 39, 40, 49, 
            50, 44, 45, 46, 47, 55, 
            56, 57, 58, 59, 60])
    >>> mg_4.at_node['drainage_area']
    array([  1.,   1.,   1.,   4.,   0.,  
             1.,   1.,   1.,   1.,   4.,   0.,
             1.,   1.,   1.,   1.,   1.,   3.,   0.,   
             4.,   1.,   1.,  24.,   8.,   1.,   2.,   0.,   
             0.,   4.,   1.,   8.,   4.,   5.,   2.,   1.,   0.,  
             0.,   3.,   1.,   5.,   3.,   2.,   1.,   0.,   
             0.,   2.,   2.,   2.,   2.,   1.,   0.,   
             0.,   1.,   1.,   1.,   1.,   0.,   
             0.,   0.,   0.,   0.,   0.])
    
    Because of the depression, the flow 'got stuck' in the hole in the center
    of the grid. We can fix this by using a depression finder, such as 
    DepressionFinderAndRouter.
    
    >>> from landlab.components import DepressionFinderAndRouter
    
    We can either run the depression finder separately from the flow 
    accumulator or we can specify the depression finder and router when we
    instantiate the accumulator and it will run automatically. 
    
    First let's try running them separately. 
    
    >>> df_4 = DepressionFinderAndRouter(mg_4)
    >>> df_4.map_depressions()
    >>> mg_4.at_node['flow__receiver_node']
    array([ 0,  1,  2,  3,  4,  
            5,  0,  1,  2,  3, 10,
           11,  5,  6, 21, 22,  9, 17, 
           18, 11, 21, 13, 21, 22, 16, 25, 
           26, 18, 29, 21, 21, 22, 31, 24, 34, 
           35, 27, 29, 30, 30, 31, 32, 42, 
           43, 36, 38, 38, 39, 40, 49, 
           50, 44, 45, 46, 47, 55, 
           56, 57, 58, 59, 60])
    >>> mg_4.at_node['drainage_area']
    array([ 25.,   1.,   1.,   4.,   0.,   
             1.,  25.,   1.,   1.,   4.,   0.,
             1.,   1.,  24.,   1.,   1.,   3.,   0.,   
             4.,   1.,   1.,  23.,   8.,   1.,   2.,   0.,   
             0.,   4.,   1.,   3.,   9.,   5.,   2.,   1.,   0.,  
             0.,   3.,   1.,   5.,   3.,   2.,   1.,   0.,   
             0.,   2.,   2.,   2.,   2.,   1.,   0.,   
             0.,   1.,   1.,   1.,   1.,   0.,   
             0.,   0.,   0.,   0.,   0.])
    
    Now the flow is routed correctly. The depression finder has properties that 
    including whether there is a lake at the node, which lake is at each node,
    the outlet node of each lake, and the area of each lake. 
    
    >>> df_4.lake_at_node
    array([False, False, False, False, False, 
           False, False, False, False, False, False, 
           False, False, False, False, False, False, False,
           False, False, False,  True,  True, False, False, False, 
           False, False, False,  True,  True,  True, False, False, False, 
           False, False, False,  True,  True, False, False, False, 
           False, False, False, False, False, False, False, 
           False, False, False, False, False, False, 
           False, False, False, False, False], dtype=bool)
    >>> df_4.lake_map
    array([-1, -1, -1, -1, -1,
           -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, 21, 21, -1, -1, -1, 
           -1, -1, -1, 21, 21, 21, -1, -1, -1, 
           -1, -1, -1, 21, 21, -1, -1, -1, 
           -1, -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, -1, -1])
    >>> df_4.lake_codes  # a unique code for each lake present on the grid
    array([21])
    >>> df_4.lake_outlets  # the outlet node of each lake in lake_codes
    array([13])
    >>> df_4.lake_areas  # the area of each lake in lake_codes
    array([ 7.])
    
    Alternatively, we can initialize a flow accumulator with a depression 
    finder specified. Calling run_one_step() will run both the accumulator
    and the depression finder with one call. 
    
    >>> mg_5 = HexModelGrid(9,5, dx)
    >>> _ = mg_5.add_field('topographic__elevation', mg_5.node_x + np.round(mg_5.node_y), at = 'node')
    >>> mg_5.at_node['topographic__elevation'][depression_ids] *= 0.1
    >>> fa_5 = FlowAccumulatorSteepestDescent(mg_5, depression_finder=DepressionFinderAndRouter)
    >>> fa_5.run_one_step()
    
    This has the same effect of first calling the accumulator and then calling
    the depression finder. 
    
    >>> mg_5.at_node['flow__receiver_node']
    array([ 0,  1,  2,  3,  4,  
            5,  0,  1,  2,  3, 10,
           11,  5,  6, 21, 22,  9, 17, 
           18, 11, 21, 13, 21, 22, 16, 25, 
           26, 18, 29, 21, 21, 22, 31, 24, 34, 
           35, 27, 29, 30, 30, 31, 32, 42, 
           43, 36, 38, 38, 39, 40, 49, 
           50, 44, 45, 46, 47, 55, 
           56, 57, 58, 59, 60])
    >>> mg_5.at_node['drainage_area']
    array([ 25.,   1.,   1.,   4.,   0.,   
             1.,  25.,   1.,   1.,   4.,   0.,
             1.,   1.,  24.,   1.,   1.,   3.,   0.,   
             4.,   1.,   1.,  23.,   8.,   1.,   2.,   0.,   
             0.,   4.,   1.,   3.,   9.,   5.,   2.,   1.,   0.,  
             0.,   3.,   1.,   5.,   3.,   2.,   1.,   0.,   
             0.,   2.,   2.,   2.,   2.,   1.,   0.,   
             0.,   1.,   1.,   1.,   1.,   0.,   
             0.,   0.,   0.,   0.,   0.])
    
    The depression finder is stored as part of the flow accumulator, so its 
    properties can be accessed through the depression finder. 
    
    >>> fa_5.df.lake_at_node
    array([False, False, False, False, False, 
           False, False, False, False, False, False, 
           False, False, False, False, False, False, False,
           False, False, False,  True,  True, False, False, False, 
           False, False, False,  True,  True,  True, False, False, False, 
           False, False, False,  True,  True, False, False, False, 
           False, False, False, False, False, False, False, 
           False, False, False, False, False, False, 
           False, False, False, False, False], dtype=bool)
    >>> fa_5.df.lake_map
    array([-1, -1, -1, -1, -1,
           -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, 21, 21, -1, -1, -1, 
           -1, -1, -1, 21, 21, 21, -1, -1, -1, 
           -1, -1, -1, 21, 21, -1, -1, -1, 
           -1, -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, -1, -1, -1, 
           -1, -1, -1, -1, -1])
    >>> fa_5.df.lake_codes  # a unique code for each lake present on the grid
    array([21])
    >>> fa_5.df.lake_outlets  # the outlet node of each lake in lake_codes
    array([13])
    >>> fa_5.df.lake_areas  # the area of each lake in lake_codes
    array([ 7.])
    
    """
    
    _name = 'FlowAccumulatorSteepestDescent'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, surface='topographic__elevation', depression_finder=None):
        self.method = 'SteepestDescent'
        super(FlowAccumulatorSteepestDescent, self).__init__(grid, surface)

        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        if not self._is_Voroni:
            raise NotImplementedError('FlowAccumulatorSteepestDescent not implemented for regular grids, use FlowAccumulatorD4 or FlowAccumulatorD8 instead')
        # save method as attribute
        
        self.fd = FlowDirector(self._grid, self.elevs)
        self.df_component = depression_finder
        
        if self.df_component:
            self.df=self.df_component(self.grid)
        
    def run_one_step(self):
        
        # step 0. Check and update BCs
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code        
        
        # step 1. Find flow directions by specified method
        self.fd.run_one_step()
        
        # step 2. Get r (and potentially p) array(s)        
        r = self._grid['node']['flow__receiver_node']
        
        # step 2. Stack, D, delta construction
        nd = flow_accum_bw._make_number_of_donors_array(r)
        delta = flow_accum_bw._make_delta_array(nd)
        D = flow_accum_bw._make_array_of_donors(r, delta)
        s = flow_accum_bw.make_ordered_node_array(r, self.fd.sink)
        
        #put theese in grid so that depression finder can use it.         
        # store the generated data in the grid
        self._grid['node']['flow__data_structure_delta'][:] = delta[1:]
        self._grid['link']['flow__data_structure_D'][:len(D)] = D
        self._grid['node']['flow__upstream_node_order'][:] = s
        
        # step 3. Initialize and Run depression finder if passed 
        # at present this must go at the end. 

       
        
        # step 4. Accumulate (to one or to N depending on direction method. )
        a, q = flow_accum_bw.find_drainage_area_and_discharge(s, 
                                                              r, 
                                                              self.node_cell_area,
                                                              self._grid.at_node['water__unit_flux_in'])
      
        
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['surface_water__discharge'][:] = q
        
        if self.df_component:
            self.df.map_depressions()


    
if __name__ == '__main__':
    import doctest
    doctest.testmod()    