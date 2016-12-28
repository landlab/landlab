from landlab.components.flow_accum.flow_accumulator import FlowAccumulator 
from landlab.components.flow_director import FlowDirectorD4 as FlowDirector
from landlab.components.flow_accum import flow_accum_bw 
from landlab import VoronoiDelaunayGrid

class FlowAccumulatorD4(FlowAccumulator):
    """
    Single-path (steepest direction) flow routing by the D4 method on raster
    grids. This method considers flow on the four links that connect a given 
    node across faces (no flow on diagonal links). The method that considers 
    diagonal links for raster grids is FlowAccumulatorD8.
    
    For Voroni Grids use FlowAccumulatorSteepestDecent. 

    This class implements single-path (steepest direction) flow routing, and
    calculates flow directions, drainage area, and discharge.
    
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
    >>> mg_2.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  1.,  5.,  0.,
            0.,  1.,  5.,  0.,  
            0.,  1.,  4.,  0.,  
            0.,  1.,  2.,  0.,  
            0.,  0.,  0.,  0.])

    Now let's change the cell area (100.) and the runoff rates:

    >>> mg_3 = RasterModelGrid((5, 4), spacing=(10., 10))

    Put the data back into the new grid.

    >>> _ = mg_3.add_field('node','topographic__elevation', elev)
    >>> mg_3.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa_3 = FlowAccumulatorD4(mg_3)
    >>> runoff_rate = np.arange(mg_3.number_of_nodes)
    >>> _ = mg_3.add_field('node', 'water__unit_flux_in', runoff_rate,
    ...                  noclobber=False)
    >>> fa_3.run_one_step()
    >>> mg_3.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
    array([    0.,   500.,  5200.,     0.,
               0.,   500.,  5200.,     0.,
               0.,   900.,  4600.,     0.,     
               0.,  1300.,  2700.,     0.,
               0.,     0.,     0.,     0.])
    
    Finally, see what happens when there is a depression. 
    
    >>> mg_4 = RasterModelGrid((7, 7), 0.5)
    >>> z = mg_4.add_field('node', 'topographic__elevation', mg_4.node_x.copy())
    >>> z += 0.01 * mg_4.node_y
    >>> mg_4.at_node['topographic__elevation'].reshape(mg_4.shape)[2:5, 2:5] *= 0.1
    >>> mg_4.set_closed_boundaries_at_grid_edges(True, True, False, True)

    This model grid has a depression in the center. 
    
    >>> mg_4.at_node['topographic__elevation'].reshape(mg_4.shape)
    array([[ 0.    ,  0.5   ,  1.    ,  1.5   ,  2.    ,  2.5   ,  3.    ],
           [ 0.005 ,  0.505 ,  1.005 ,  1.505 ,  2.005 ,  2.505 ,  3.005 ],
           [ 0.01  ,  0.51  ,  0.101 ,  0.151 ,  0.201 ,  2.51  ,  3.01  ],
           [ 0.015 ,  0.515 ,  0.1015,  0.1515,  0.2015,  2.515 ,  3.015 ],
           [ 0.02  ,  0.52  ,  0.102 ,  0.152 ,  0.202 ,  2.52  ,  3.02  ],
           [ 0.025 ,  0.525 ,  1.025 ,  1.525 ,  2.025 ,  2.525 ,  3.025 ],
           [ 0.03  ,  0.53  ,  1.03  ,  1.53  ,  2.03  ,  2.53  ,  3.03  ]])
    >>> fa_4 = FlowAccumulatorD4(mg_4) 
    >>> fa_4.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg_4.at_node['flow__receiver_node'].reshape(mg_4.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14, 16, 16, 17, 18, 20],
           [21, 21, 16, 23, 24, 25, 27],
           [28, 28, 23, 30, 31, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg_4.at_node['drainage_area'].reshape(mg_4.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  5.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  3.  ,  0.75,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  2.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    
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
    >>> mg_4.at_node['flow__receiver_node'].reshape(mg_4.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg_4.at_node['drainage_area'].reshape(mg_4.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 5.25,  5.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  5.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.75,  2.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.5 ,  0.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    
    Now the flow is routed correctly. The depression finder has properties that 
    including whether there is a lake at the node, which lake is at each node,
    the outlet node of each lake, and the area of each lake. 
    
    >>> df_4.lake_at_node.reshape(mg_4.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]], dtype=bool)
    >>> df_4.lake_map.reshape(mg_4.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1]])
    >>> df_4.lake_codes  # a unique code for each lake present on the grid
    array([16])
    >>> df_4.lake_outlets  # the outlet node of each lake in lake_codes
    array([8])
    >>> df_4.lake_areas  # the area of each lake in lake_codes
    array([ 2.25])
    
    Alternatively, we can initialize a flow accumulator with a depression 
    finder specified. Calling run_one_step() will run both the accumulator
    and the depression finder with one call. 
    
    >>> mg_5 = RasterModelGrid((7, 7), 0.5)
    >>> z = mg_5.add_field('node', 'topographic__elevation', mg_5.node_x.copy())
    >>> z += 0.01 * mg_5.node_y
    >>> mg_5.at_node['topographic__elevation'].reshape(mg_5.shape)[2:5, 2:5] *= 0.1
    >>> fa_5 = FlowAccumulatorD4(mg_5, depression_finder=DepressionFinderAndRouter)
    >>> fa_5.run_one_step()
    
    This has the same effect of first calling the accumulator and then calling
    the depression finder. 
    
    >>> mg_5.at_node['flow__receiver_node'].reshape(mg_5.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg_5.at_node['drainage_area'].reshape(mg_5.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 5.25,  5.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  5.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.75,  2.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.5 ,  0.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    
    The depression finder is stored as part of the flow accumulator, so its 
    properties can be accessed through the depression finder. 
    
    >>> fa_5.df.lake_at_node.reshape(mg_5.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]], dtype=bool)
    >>> fa_5.df.lake_map.reshape(mg_5.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1]])
    >>> fa_5.df.lake_codes  # a unique code for each lake present on the grid
    array([16])
    >>> fa_5.df.lake_outlets  # the outlet node of each lake in lake_codes
    array([8])
    >>> fa_5.df.lake_areas  # the area of each lake in lake_codes
    array([ 2.25])
    """
    
    _name = 'FlowAccumulatorD4'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, surface='topographic__elevation', depression_finder=None):
        super(FlowAccumulatorD4, self).__init__(grid, surface)

        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        if self._is_Voroni:
            raise NotImplementedError('FlowAccumulatorD4 not implemented for irregular grids, use FlowAccumulatorSteepestDecent')
        # save method as attribute
        self.method = 'D4'
        self.fd = FlowDirector(self._grid, self.elevs)
        self.df_component = depression_finder
        
    def run_one_step(self):
        
        # step 0. Check and update BCs
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code        
        
        # step 1. Find flow directions by specified method
        self.fd.run_one_step()
        self.baselevel_nodes = self.fd.baselevel_nodes
        self.sink = self.fd.sink
        
        # step 2. Get r (and potentially p) array(s)        
        r = self._grid['node']['flow__receiver_node']
        
        # step 2. Stack, D, delta construction
        nd = flow_accum_bw._make_number_of_donors_array(r)
        delta = flow_accum_bw._make_delta_array(nd)
        D = flow_accum_bw._make_array_of_donors(r, delta)
        s = flow_accum_bw.make_ordered_node_array(r, self.sink)
        
        #put theese in grid so that depression finder can use it.         
        # store the generated data in the grid
        self._grid['node']['flow__data_structure_delta'][:] = delta[1:]
        self._grid['link']['flow__data_structure_D'][:len(D)] = D
        self._grid['node']['flow__upstream_node_order'][:] = s
        
        # step 3. Initialize and Run depression finder if passed 
        # at present this might need to go at the very end... also need to 
        #make sure that the df. properties work. 
  
        
        # step 4. Accumulate (to one or to N depending on direction method. )
        a, q = flow_accum_bw.find_drainage_area_and_discharge(s, 
                                                              r, 
                                                              self.node_cell_area,
                                                              self._grid.at_node['water__unit_flux_in'])
  
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['surface_water__discharge'][:] = q
        
        if self.df_component:
            self.df=self.df_component(self.grid)
            self.df.map_depressions()


    
if __name__ == '__main__':
    import doctest
    doctest.testmod()    