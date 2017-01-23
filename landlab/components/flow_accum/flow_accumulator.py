#!/usr/env/python

"""
flow_accumulator.py: Component to accumulate flow and calculate drainage area.

Provides the FlowAccumulator component which accumulates flow and calculates
drainage area. FlowAccumulator supports multiple methods for calculating flow
direction. Optionally a depression finding component can be specified and flow
directing, depression finding, and flow routing can all be accomplished
together.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.decorators import use_field_name_or_array

from landlab.components.flow_accum import flow_accum_bw
from landlab.components.flow_accum import flow_accum_to_n

from landlab import BAD_INDEX_VALUE
import six


@use_field_name_or_array('node')
def _return_surface(grid, surface):

    """
    Private function to return the surface to direct flow over.

    This function exists to take advantange of the 'use_field_name_or_array
    decorator which permits providing the surface as a field name or array.
    """
    return(surface)


class FlowAccumulator(Component):

    """
    Component to accumulate flow and calculate drainage area.

    This is accomplished by first finding flow directions by a user-specified
    method and then calculating the drainage area and discharge.

    Optionally, spatially variable runoff can be set either by the model grid
    field 'water__unit_flux_in' or the input variable *runoff_rate**.

    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished
    together.


    NOTE: The perimeter nodes  NEVER contribute to the accumulating flux, even
    if the  gradients from them point inwards to the main body of the grid.
    This is because under Landlab definitions, perimeter nodes lack cells, so
    cannot accumulate any discharge.

    FlowAccumulator stores as ModelGrid fields:

        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure:
            *flow__data_structure_delta*. The first element is always zero.
        -  Link array of the D data structure: *flow__data_structure_D*

    The FlowDirector component will add additional ModelGrid fields.
    DirectToOne methods(Steepest/D4 and D8) and DirectToMany(NAMES HERE) use
    different model grid fields.

    DirectToOne Methods (Steeptest/D4 and D8) store the following as ModelGrid
    fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Node array containing ID of link that leads from each node to its
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    DirectToMany Methods (NAMES HERE) store the following as ModelGrid
    fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_nodes'*
        -  Node array of receiver proportions: *'flow__receiver_proportions'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`

    `run_one_step` takes the optional argument update_flow_director (default is
    True) that determines if the flow_director is re-run before flow is
    accumulated.

    Parameters
    ----------
    grid : ModelGrid
        A grid of type Voroni.
    surface : field name at node or array of length node
        The surface to direct flow across.
    flow_director : string, class, instance of class.
        A string of method or class name (e.g. 'D8' or 'FlowDirectorD8'), an
        uninstantiated FlowDirector class, or an instance of a FlowDirector
        class. This sets the method used to calculate flow directions.
        Default is 'FlowDirectorSteepest'
    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.
    depression_finder : string, class, instance of class, optional
         A string of class name (e.g., 'DepressionFinderAndRouter'), an
         uninstantiated DepressionFinder class, or an instance of a
         DepressionFinder class.
         This sets the method for depression finding.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at = 'node')

    The FlowAccumulator component accumulates flow and calculates drainage using
    all of the different methods for directing flow in Landlab. These include
    steepest descent (also known as D4 for the case of a raster grid) and D8 (on
    raster grids only). The method for flow director can be specified as a
    string (e.g., 'D8' or 'FlowDirectorD8'), as an uninstantiated FlowDirector
    component or as an instantiated FlowDirector component.

    The default method is to use FlowDirectorSteepest.

    First let's look at the three ways to instantiate a FlowAccumulator. The
    following four methods are all equivalent. First, we can pass the entire
    name of a flow director as a string to the argument `flow_director`:

    >>> fa = FlowAccumulator(mg, 'topographic__elevation',
    ...                      flow_director='FlowDirectorSteepest')

    Second, we can pass just the method name as a string to the argument
    `flow_director`:

    >>> fa = FlowAccumulator(mg, 'topographic__elevation',
    ...                      flow_director='Steepest')

    Third, we can import a FlowDirector component from Landlab and pass it to
    `flow_director`:
    >>> from landlab.components import FlowDirectorSteepest
    >>> fa = FlowAccumulator(mg, 'topographic__elevation',
    ...                      flow_director=FlowDirectorSteepest)

    Finally, we can instantiate a FlowDirector component and pass this
    instantiated version to `flow_director`. You might want to do this if you
    used a FlowDirector in order to set up something before starting a
    time loop and then want to use the same flow director within the loop.

    >>> fd = FlowDirectorSteepest(mg, 'topographic__elevation')
    >>> fa = FlowAccumulator(mg, 'topographic__elevation',
    ...                      flow_director=FlowDirectorSteepest)

    Now let's look at what FlowAccumulator does. Even before we run
    FlowAccumulator it has the property `surface_values` that stores the values
    of the surface over which flow is directed and accumulated.

    >>> fa.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])

    Now let's make a more complicated elevation grid for the next examples.

    >>> mg_2 = RasterModelGrid((5, 4), spacing=(1, 1))
    >>> topographic__elevation = np.array([0.,  0.,  0., 0.,
    ...                                    0., 21., 10., 0.,
    ...                                    0., 31., 20., 0.,
    ...                                    0., 32., 30., 0.,
    ...                                    0.,  0.,  0., 0.])
    >>> _ = mg_2.add_field('node',
    ...                    'topographic__elevation',
    ...                    topographic__elevation)
    >>> mg_2.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa_2 = FlowAccumulator(mg_2, 'topographic__elevation',
    ...                        flow_director=FlowDirectorSteepest)
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

    >>> _ = mg_3.add_field('node',
    ...                    'topographic__elevation',
    ...                    topographic__elevation)
    >>> mg_3.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa_3 = FlowAccumulator(mg_3, 'topographic__elevation',
    ...                        flow_director=FlowDirectorSteepest)
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

    The FlowAccumulator component will work for both raster grids and irregular
    grids. For the example we will use a Hexagonal Model Grid, a special type of
    Voroni Grid that has regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid(5,3)
    >>> _ = hmg.add_field('topographic__elevation',
    ...                   hmg.node_x + np.round(hmg.node_y),
    ...                   at = 'node')
    >>> fa = FlowAccumulator(hmg, 'topographic__elevation',
    ...                      flow_director=FlowDirectorSteepest)
    >>> fa.surface_values
    array([ 0. ,  1. ,  2. ,
            0.5,  1.5,  2.5,  3.5,
            1. ,  2. ,  3. ,  4. ,  5. ,
            2.5,  3.5,  4.5,  5.5,
            3. ,  4. ,  5. ])

    Next, let's set the dx spacing such that each cell has an area of one.

    >>> dx=(2./(3.**0.5))**0.5
    >>> hmg_2 = HexModelGrid(5,3, dx)
    >>> _ = hmg_2.add_field('topographic__elevation',
    ...                     hmg_2.node_x**2 + np.round(hmg_2.node_y)**2,
    ...                     at = 'node')
    >>> fa_2 = FlowAccumulator(hmg_2, 'topographic__elevation',
    ...                        flow_director=FlowDirectorSteepest)
    >>> fa_2.run_one_step()
    >>> hmg_2.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,
            3,  0,  1,  6,
            7,  3,  4,  5, 11,
           12,  8,  9, 15,
           16, 17, 18])
    >>> hmg_2.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 3.,  2.,  0.,
            2.,  3.,  2.,  0.,
            0.,  2.,  2.,  1.,  0.,
            0., 1.,  1.,  0.,
            0.,  0.,  0.])

    Now let's change the cell area (100.) and the runoff rates:

    >>> hmg_3 = HexModelGrid(5,3, dx*10.)

    Put the data back into the new grid.

    >>> _ = hmg_3.add_field('topographic__elevation',
    ...                     hmg_3.node_x**2 + np.round(hmg_3.node_y)**2,
    ...                     at = 'node')
    >>> fa_3 = FlowAccumulator(hmg_3, 'topographic__elevation',
    ...                        flow_director=FlowDirectorSteepest)
    >>> fa_3.run_one_step()
    >>> hmg_3.at_node['surface_water__discharge']
    array([ 500.,    0.,    0.,
            200.,  500.,  200.,    0.,
              0.,  200.,  200.,  100.,    0.,
              0.,  100.,  100.,    0.,
              0.,    0.,    0.])

    Next, let's see what happens to a raster grid when there is a depression.

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
    >>> fa_4 = FlowAccumulator(mg_4, 'topographic__elevation',
    ...                        flow_director=FlowDirectorSteepest)
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
    instantiate the accumulator and it will run automatically. Similar to
    specifying the FlowDirector we can provide a depression finder in multiple
    three ways.

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
    and the depression finder with one call. For this example, we will pass the
    class DepressionFinderAndRouter to the parameter `depression_finder`.

    >>> mg_5 = RasterModelGrid((7, 7), 0.5)
    >>> z = mg_5.add_field('node', 'topographic__elevation', mg_5.node_x.copy())
    >>> z += 0.01 * mg_5.node_y
    >>> mg_5.at_node['topographic__elevation'].reshape(mg_5.shape)[2:5, 2:5] *= 0.1
    >>> fa_5 = FlowAccumulator(mg_5, 'topographic__elevation',
    ...                        flow_director=FlowDirectorSteepest,
    ...                        depression_finder=DepressionFinderAndRouter)
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

    >>> fa_5.depression_finder.lake_at_node.reshape(mg_5.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]], dtype=bool)
    >>> fa_5.depression_finder.lake_map.reshape(mg_5.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1]])
    >>> fa_5.depression_finder.lake_codes  # a unique code for each lake present on the grid
    array([16])
    >>> fa_5.depression_finder.lake_outlets  # the outlet node of each lake in lake_codes
    array([8])
    >>> fa_5.depression_finder.lake_areas  # the area of each lake in lake_codes
    array([ 2.25])
    """

    _name = 'FlowAccumulator'

    _input_var_names = ('topographic__elevation',
                        'water__unit_flux_in'
                        )

    _output_var_names = ('drainage_area',
                         'surface_water__discharge',
                         'flow__upstream_node_order',
                         'flow__nodes_not_in_stack',
                         'flow__data_structure_delta',
                         'flow__data_structure_D'
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_node': 'm',
                  'water__unit_flux_in': 'm/s',
                  'drainage_area': 'm**2',
                  'surface_water__discharge': 'm**3/s',
                  'flow__upstream_node_order': '-',
                  'flow__data_structure_delta': '-',
                  'flow__data_structure_D': '-',
                  'flow__nodes_not_in_stack': '-'
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'water__unit_flux_in': 'node',
                    'drainage_area': 'node',
                    'surface_water__discharge': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__nodes_not_in_stack': 'grid',
                    'flow__data_structure_delta': 'node',
                    'flow__data_structure_D': 'link',
                    }
    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'surface_water__discharge': 'Discharge of water through each node',
        'water__unit_flux_in':
            'External volume water per area per time input to each node '
            '(e.g., rainfall rate)',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'flow__data_structure_delta':
            'Node array containing the elements delta[1:] of the data '
            'structure "delta" used for construction of the downstream-to-'
            'upstream node array',
        'flow__data_structure_D':
            'Link array containing the data structure D used for construction'
            'of the downstream-to-upstream node array',
        'flow__nodes_not_in_stack':
            'Boolean value indicating if there are any nodes that have not yet'
            'been added to the stack stored in flow__upstream_node_order.'
            }

    def __init__(self,
                 grid,
                 surface='topographic__elevation',
                 flow_director='FlowDirectorD4',
                 runoff_rate=None,
                 depression_finder=None):
        """
        Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and compatability
        for the flow_director and depression_finder keyword arguments, tests
        the argument of runoff_rate, and initializes new fields.
        """
        super(FlowAccumulator, self).__init__(grid)
        # Keep a local reference to the grid
        self._grid = grid

        # Grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)

        # STEP 1: Testing of input values, supplied either in function call or
        # as part of the grid.
        self.test_water_inputs(grid, runoff_rate)

        # save elevations and node_cell_area to class properites.
        self.surface = surface
        self.surface_values = _return_surface(grid, surface)

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.

        self.node_cell_area = node_cell_area

        # STEP 2:
        # identify Flow Director method, save name, import and initialize the correct
        # flow director component if necessary
        self.add_director(flow_director)
        self.add_depression_finder(depression_finder)

        # This component will track of the following variables.
        # Attempt to create each, if they already exist, assign the existing
        # version to the local copy.

        #   - drainage area at each node
        #   - receiver of each node
        #   - D array
        #   - delta array
        #   - missing nodes in stack.
        try:
            self.drainage_area = grid.add_zeros('drainage_area', at='node',
                                                dtype=float)
        except FieldError:
            self.drainage_area = grid.at_node['drainage_area']

        try:
            self.discharges = grid.add_zeros('surface_water__discharge',
                                             at='node', dtype=float)
        except FieldError:
            self.discharges = grid.at_node['surface_water__discharge']

        try:
            self.upstream_ordered_nodes = grid.add_field('flow__upstream_node_order',
                                                         BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                         at='node', dtype=int)

        except FieldError:
            self.upstream_ordered_nodes = grid.at_node[
                'flow__upstream_node_order']

        try:
            self.delta_structure = grid.add_field('flow__data_structure_delta',
                                                  BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                  at='node', dtype=int)
        except FieldError:
            self.delta_structure = grid.at_node['flow__data_structure_delta']

        try:
            # needs to be BAD_INDEX_VALUE
            self.D_structure = grid.add_field('flow__data_structure_D',
                                              BAD_INDEX_VALUE*grid.ones(at='link'),
                                              at='link', dtype=int)
        except FieldError:
            self.D_structure = grid.at_link['flow__data_structure_D']

        self.nodes_not_in_stack = True

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid['node']['drainage_area']

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid['node']['surface_water__discharge']

    @property
    def node_order_upstream(self):
        """Return the upstream node order (drainage stack)."""
        return self._grid['node']['flow__upstream_node_order']

    def test_water_inputs(self, grid, runoff_rate):
        """Test inputs for runoff_rate and water__unit_flux_in."""
        # testing input for runoff rate, can be None, a string associated with
        # a field at node, a single float or int, or an array of size number of
        # nodes.
        if runoff_rate is not None:
            if type(runoff_rate) is str:
                runoff_rate = grid.at_node[runoff_rate]
            elif type(runoff_rate) in (float, int):
                pass
            else:
                assert runoff_rate.size == grid.number_of_nodes

        # test for water__unit_flux_in
        try:
            grid.at_node['water__unit_flux_in']
        except FieldError:
            if runoff_rate is None:
                # assume that if runoff rate is not supplied, that the value
                # should be set to one everywhere.
                grid.add_ones('node', 'water__unit_flux_in', dtype=float)
            else:
                if type(runoff_rate) in (float, int):
                    grid.add_empty('node', 'water__unit_flux_in', dtype=float)
                    grid.at_node['water__unit_flux_in'].fill(runoff_rate)
                else:
                    grid.at_node['water__unit_flux_in'] = runoff_rate
        else:
            if runoff_rate is not None:
                print ("FlowAccumulator found both the field " +
                       "'water__unit_flux_in' and a provided float or " +
                       "array for the runoff_rate argument. THE FIELD IS " +
                       "BEING OVERWRITTEN WITH THE SUPPLIED RUNOFF_RATE!")
                if type(runoff_rate) in (float, int):
                    grid.at_node['water__unit_flux_in'].fill(runoff_rate)
                else:
                    grid.at_node['water__unit_flux_in'] = runoff_rate

        # perform a test (for politeness!) that the old name for the water_in
        # field is not present:
        try:
            grid.at_node['water__discharge_in']
        except FieldError:
            pass
        else:
            warnings.warn("This component formerly took 'water__discharge" +
                          "_in' as an input field. However, this field is " +
                          "now named 'water__unit_flux_in'. You are still " +
                          "using a field with the old name. Please update " +
                          "your code if you intended the FlowRouter to use " +
                          "that field.", DeprecationWarning)

    def add_director(self, flow_director):
        """Test and add the flow director component."""
        PERMITTED_DIRECTORS = ['FlowDirectorSteepest',
                               'FlowDirectorD8']

        # flow director is provided as a string.
        if isinstance(flow_director, six.string_types):
            if flow_director[:12] == 'FlowDirector':
                flow_director = flow_director[12:]

            from landlab.components.flow_director import FlowDirectorSteepest, FlowDirectorD8
            DIRECTOR_METHODS = {'D4': FlowDirectorSteepest,
                                'Steepest': FlowDirectorSteepest,
                                'D8': FlowDirectorD8,
                                }

            try:
                FlowDirector = DIRECTOR_METHODS[flow_director]
            except KeyError:
                raise ValueError('String provided in flow_director is not a valid method or component name. '
                                 'The following components are valid imputs:\n'+str(PERMITTED_DIRECTORS))

            self.flow_director = FlowDirector(self._grid, self.surface_values)
        # flow director is provided as an instantiated flow director
        elif isinstance(flow_director, Component):
             if flow_director._name in PERMITTED_DIRECTORS:
                 self.flow_director = flow_director
             else:
                 raise ValueError('Component provided in flow_director is not a valid component. '
                                 'The following components are valid imputs:\n'+str(PERMITTED_DIRECTORS))
        # flow director is provided as an uninstantiated flow director
        else:

            if flow_director._name in PERMITTED_DIRECTORS:
                FlowDirector = flow_director
                self.flow_director = FlowDirector(self._grid, self.surface_values)
            else:
                raise ValueError('Component provided in flow_director is not a valid component. '
                                 'The following components are valid imputs:\n'+str(PERMITTED_DIRECTORS))

        # save method as attribute
        self.method = self.flow_director.method
    def add_depression_finder(self,depression_finder):
        """Test and add the depression finder component."""
        PERMITTED_DEPRESSION_FINDERS = ['DepressionFinderAndRouter']

        # now do a similar thing for the depression finder.
        self.depression_finder_provided = depression_finder
        if self.depression_finder_provided:
            # depression finder is provided as a string.
            if isinstance(self.depression_finder_provided, six.string_types):

                from landlab.components import DepressionFinderAndRouter
                DEPRESSION_METHODS = {'DepressionFinderAndRouter': DepressionFinderAndRouter
                                    }

                try:
                    DepressionFinder = DEPRESSION_METHODS[self.depression_finder_provided]
                except KeyError:
                    raise ValueError('Component provided in depression_finder is not a valid component. '
                                     'The following components are valid imputs:\n'+str(PERMITTED_DEPRESSION_FINDERS))

                self.depression_finder = DepressionFinder(self._grid)
            # flow director is provided as an instantiated depression finder
            elif isinstance(self.depression_finder_provided, Component):

                if self.depression_finder_provided._name in PERMITTED_DEPRESSION_FINDERS:
                    self.depression_finder = self.depression_finder_provided
                else:
                    raise ValueError('Component provided in depression_finder is not a valid component. '
                                     'The following components are valid imputs:\n'+str(PERMITTED_DEPRESSION_FINDERS))
            # depression_fiuner is provided as an uninstantiated depression finder
            else:

                if self.depression_finder_provided._name in PERMITTED_DEPRESSION_FINDERS:
                    DepressionFinder = self.depression_finder_provided
                    self.depression_finder = DepressionFinder(self._grid)
                else:
                    raise ValueError('Component provided in depression_finder is not a valid component. '
                    'The following components are valid imputs:\n'+str(PERMITTED_DEPRESSION_FINDERS))


    def run_one_step(self, update_flow_director=True):
        """
        Function to make FlowAccumulator calculate drainage area and discharge.

        Running run_one_step() results in the following to occur:
            1. Flow directions are updated (unless update_flow_director is set
            as False).
            2. Intermediate steps that analyse the drainage network topology
            and create datastructures for efficient drainage area and discharge
            calculations.
            3. Calculation of drainage area and discharge.
            4. Depression finding and mapping, which updates drainage area and
            discharge.
        """
        # step 1. Find flow directions by specified method
        if update_flow_director == True:
            self.flow_director.run_one_step()

        # further steps vary depending on how many recievers are present
        # one set of steps is for route to one (D8, Steepest/D4)
        if self.flow_director.to_n_receivers == 'one':

            # step 2. Get r (and potentially p) array(s)
            r = self._grid['node']['flow__receiver_node']

            # step 2. Stack, D, delta construction
            nd = flow_accum_bw._make_number_of_donors_array(r)
            delta = flow_accum_bw._make_delta_array(nd)
            D = flow_accum_bw._make_array_of_donors(r, delta)
            s = flow_accum_bw.make_ordered_node_array(r)

            # put theese in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid['node']['flow__data_structure_delta'][:] = delta[1:]
            self._grid['link']['flow__data_structure_D'][:len(D)] = D
            self._grid['node']['flow__upstream_node_order'][:] = s

            # step 3. Run depression finder if passed
            # at present this must go at the end.

            # step 4. Accumulate (to one or to N depending on direction method. )
            a, q = flow_accum_bw.find_drainage_area_and_discharge(s,
                                                                  r,
                                                                  self.node_cell_area,
                                                                  self._grid.at_node['water__unit_flux_in'])

            self._grid['node']['drainage_area'][:] = a
            self._grid['node']['surface_water__discharge'][:] = q

        # at the moment, this is where the depression finder needs to live.
        if self.depression_finder_provided:
            self.depression_finder.map_depressions()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
