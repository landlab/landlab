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

import numpy as np
import six

from landlab import (  # for type tests
    BAD_INDEX_VALUE,
    Component,
    FieldError,
    NetworkModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
)
from landlab.components.flow_accum import flow_accum_bw, flow_accum_to_n
from landlab.core.messages import warning_message
from landlab.core.utils import as_id_array
from landlab.utils.return_array import return_array_at_node


class FlowAccumulator(Component):

    """Component to accumulate flow and calculate drainage area.

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
        -  At Grid: D data structure: *flow__data_structure_D*

    The FlowDirector component will add additional ModelGrid fields.
    DirectToOne methods(Steepest/D4 and D8) and DirectToMany(DINF and MFD) use
    the same model grid field names. Some of these fields will be different
    shapes if a DirectToOne or a DirectToMany method is used.

    The FlowDirectors store the following as ModelGrid fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
            there is no receiver: *'flow__receiver_node'*. This array is 2D for
            RouteToMany methods and has the shape
            (n-nodes x max number of receivers).
        -  Node array of flow proportions: *'flow__receiver_proportions'*. This
            array is 2D, for RouteToMany methods and has the shape
            (n-nodes x max number of receivers).
        -  Node array of links carrying flow:  *'flow__link_to_receiver_node'*.
            This array is 2D for RouteToMany methods and has the shape
            (n-nodes x max number of receivers).
        -  Node array of downhill slopes from each receiver:
            *'topographic__steepest_slope'* This array is 2D for RouteToMany
            methods and has the shape (n-nodes x max number of receivers).
        -  Boolean node array of all local lows: *'flow__sink_flag'*
        -  Link array identifing if flow goes with (1) or against (-1) the link
           direction: *'flow_link_direction'*

    The primary method of this class is :func:`run_one_step`.

    `run_one_step` takes the optional argument update_flow_director (default is
    True) that determines if the flow_director is re-run before flow is
    accumulated.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid.
    surface : field name at node or array of length node
        The surface to direct flow across.
    flow_director : string, class, instance of class.
        A string of method or class name (e.g. 'D8' or 'FlowDirectorD8'), an
        uninstantiated FlowDirector class, or an instance of a FlowDirector
        class. This sets the method used to calculate flow directions.
        Default is 'FlowDirectorSteepest'
    runoff_rate : field name, array, or float, optional (m/time)
        If provided, sets the runoff rate and will be assigned to the grid field
        'water__unit_flux_in'. If a spatially and and temporally variable runoff
        rate is desired, pass this field name and update the field through model
        run time. If both the field and argument are present at the time of
        initialization, runoff_rate will *overwrite* the field. If neither are
        set, defaults to spatially constant unit input.
        Both a runoff_rate array and the 'water__unit_flux_in' field are
        permitted to contain negative values, in which case they mimic
        transmission losses rather than e.g. rain inputs.
    depression_finder : string, class, instance of class, optional
         A string of class name (e.g., 'DepressionFinderAndRouter'), an
         uninstantiated DepressionFinder class, or an instance of a
         DepressionFinder class.
         This sets the method for depression finding.
    **kwargs : any additional parameters to pass to a FlowDirector or
         DepressionFinderAndRouter instance (e.g., partion_method for
         FlowDirectorMFD). This will have no effect if an instantiated component
         is passed using the flow_director or depression_finder keywords.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> mg = RasterModelGrid((3,3))
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

    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director='FlowDirectorSteepest'
    ... )

    Second, we can pass just the method name as a string to the argument
    `flow_director`:

    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director='Steepest'
    ... )

    Third, we can import a FlowDirector component from Landlab and pass it to
    `flow_director`:
    >>> from landlab.components import FlowDirectorSteepest
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ... )

    Finally, we can instantiate a FlowDirector component and pass this
    instantiated version to `flow_director`. You might want to do this if you
    used a FlowDirector in order to set up something before starting a
    time loop and then want to use the same flow director within the loop.

    >>> fd = FlowDirectorSteepest(mg, 'topographic__elevation')
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ... )

    Now let's look at what FlowAccumulator does. Even before we run
    FlowAccumulator it has the property `surface_values` that stores the values
    of the surface over which flow is directed and accumulated.

    >>> fa.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])

    Now let's make a more complicated elevation grid for the next examples.

    >>> mg = RasterModelGrid((5, 4))
    >>> topographic__elevation = np.array([0.,  0.,  0., 0.,
    ...                                    0., 21., 10., 0.,
    ...                                    0., 31., 20., 0.,
    ...                                    0., 32., 30., 0.,
    ...                                    0.,  0.,  0., 0.])
    >>> _ = mg.add_field(
    ...     'node',
    ...     'topographic__elevation',
    ...     topographic__elevation
    ... )
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ...      )
    >>> fa.run_one_step()
    >>> mg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,  3,
            4,  1,  2,  7,
            8, 10,  6, 11,
           12, 14, 10, 15,
           16, 17, 18, 19])
    >>> mg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  1.,  5.,  0.,
            0.,  1.,  5.,  0.,
            0.,  1.,  4.,  0.,
            0.,  1.,  2.,  0.,
            0.,  0.,  0.,  0.])

    Now let's change the cell area (100.) and the runoff rates:

    >>> mg = RasterModelGrid((5, 4), xy_spacing=(10., 10))

    Put the data back into the new grid.

    >>> _ = mg.add_field(
    ...     'node',
    ...     'topographic__elevation',
    ...     topographic__elevation
    ... )
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ...      )
    >>> runoff_rate = np.arange(mg.number_of_nodes, dtype=float)
    >>> rnff = mg.add_field(
    ...        'node',
    ...        'water__unit_flux_in',
    ...        runoff_rate,
    ...        noclobber=False
    ... )
    >>> fa.run_one_step()
    >>> mg.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
    array([    0.,   500.,  5200.,     0.,
               0.,   500.,  5200.,     0.,
               0.,   900.,  4600.,     0.,
               0.,  1300.,  2700.,     0.,
               0.,     0.,     0.,     0.])

    The flow accumulator will happily work with a negative runoff rate, which
    could be used to allow, e.g., transmission losses:

    >>> runoff_rate.fill(1.)
    >>> fa.run_one_step()
    >>> mg.at_node['surface_water__discharge']
    array([   0.,  100.,  500.,    0.,
              0.,  100.,  500.,    0.,
              0.,  100.,  400.,    0.,
              0.,  100.,  200.,    0.,
              0.,    0.,    0.,    0.])
    >>> runoff_rate[:8] = -0.5
    >>> fa.run_one_step()
    >>> mg.at_node['surface_water__discharge']
    array([   0.,    0.,  350.,    0.,
              0.,    0.,  350.,    0.,
              0.,  100.,  400.,    0.,
              0.,  100.,  200.,    0.,
              0.,    0.,    0.,    0.])

    The drainage area array is unaffected, as you would expect:

    >>> mg.at_node['drainage_area']
    array([   0.,  100.,  500.,    0.,
              0.,  100.,  500.,    0.,
              0.,  100.,  400.,    0.,
              0.,  100.,  200.,    0.,
              0.,    0.,    0.,    0.])

    The FlowAccumulator component will work for both raster grids and irregular
    grids. For the example we will use a Hexagonal Model Grid, a special type
    of Voroni Grid that has regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid(5,3, xy_of_lower_left=(-1., 0.))
    >>> _ = hmg.add_field(
    ...     'topographic__elevation',
    ...     hmg.node_x + np.round(hmg.node_y),
    ...     at = 'node'
    ...     )
    >>> fa = FlowAccumulator(
    ...      hmg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ... )
    >>> fa.surface_values
    array([ 0. ,  1. ,  2. ,
            0.5,  1.5,  2.5,  3.5,
            1. ,  2. ,  3. ,  4. ,  5. ,
            2.5,  3.5,  4.5,  5.5,
            3. ,  4. ,  5. ])

    If the FlowDirector you want to use takes keyword arguments and you want
    to specify it using a string or uninstantiated FlowDirector class, include
    those keyword arguments when you create FlowAccumulator.

    For example, in the case of a raster grid, FlowDirectorMFD can use only
    orthogonal links, or it can use both orthogonal and diagonal links.

    >>> mg = RasterModelGrid((5, 5))
    >>> topographic__elevation = mg.node_y+mg.node_x
    >>> _ = mg.add_field(
    ...     'node',
    ...     'topographic__elevation',
    ...     topographic__elevation
    ... )
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director='MFD',
    ...      diagonals = True
    ... )
    >>> fa.run_one_step()
    >>> mg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([[ 0, -1, -1, -1, -1, -1, -1, -1],
           [ 1, -1, -1, -1, -1, -1, -1, -1],
           [ 2, -1, -1, -1, -1, -1, -1, -1],
           [ 3, -1, -1, -1, -1, -1, -1, -1],
           [ 4, -1, -1, -1, -1, -1, -1, -1],
           [ 5, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1,  5,  1, -1, -1,  0, -1],
           [-1, -1,  6,  2, -1, -1,  1, -1],
           [-1, -1,  7,  3, -1, -1,  2, -1],
           [ 9, -1, -1, -1, -1, -1, -1, -1],
           [10, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 10,  6, -1, -1,  5, -1],
           [-1, -1, 11,  7, -1, -1,  6, -1],
           [-1, -1, 12,  8, -1, -1,  7, -1],
           [14, -1, -1, -1, -1, -1, -1, -1],
           [15, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 15, 11, -1, -1, 10, -1],
           [-1, -1, 16, 12, -1, -1, 11, -1],
           [-1, -1, 17, 13, -1, -1, 12, -1],
           [19, -1, -1, -1, -1, -1, -1, -1],
           [20, -1, -1, -1, -1, -1, -1, -1],
           [21, -1, -1, -1, -1, -1, -1, -1],
           [22, -1, -1, -1, -1, -1, -1, -1],
           [23, -1, -1, -1, -1, -1, -1, -1],
           [24, -1, -1, -1, -1, -1, -1, -1]])
    >>> mg.at_node['drainage_area'].round(4) # doctest: +NORMALIZE_WHITESPACE
    array([ 1.4117,  2.065 ,  1.3254,  0.4038,  0.    ,
            2.065 ,  3.4081,  2.5754,  1.3787,  0.    ,
            1.3254,  2.5754,  2.1716,  1.2929,  0.    ,
            0.4038,  1.3787,  1.2929,  1.    ,  0.    ,
            0.    ,  0.    ,  0.    ,  0.    ,  0.    ])

    It may seem odd that there are no round numbers in the drainage area field.
    This is because flow is directed to all downhill boundary nodes and
    partitioned based on slope.

    To check that flow is conserved, sum along all boundary nodes.

    >>> round(sum(mg.at_node['drainage_area'][mg.boundary_nodes]), 4)
    9.0

    This should be the same as the number of core nodes --- as boundary nodes
    in landlab do not have area.

    >>> len(mg.core_nodes)
    9

    Next, let's set the dx spacing such that each cell has an area of one.

    >>> dx=(2./(3.**0.5))**0.5
    >>> hmg = HexModelGrid(5,3, dx, xy_of_lower_left=(-1.0745, 0.))
    >>> _ = hmg.add_field(
    ...     'topographic__elevation',
    ...     hmg.node_x**2 + np.round(hmg.node_y)**2,
    ...     at = 'node'
    ... )
    >>> fa = FlowAccumulator(
    ...      hmg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ... )
    >>> fa.run_one_step()
    >>> hmg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,
            3,  0,  1,  6,
            7,  3,  4,  5, 11,
           12,  8,  9, 15,
           16, 17, 18])
    >>> hmg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 3.,  2.,  0.,
            2.,  3.,  2.,  0.,
            0.,  2.,  2.,  1.,  0.,
            0., 1.,  1.,  0.,
            0.,  0.,  0.])

    Now let's change the cell area (100.) and the runoff rates:

    >>> hmg = HexModelGrid(5,3, dx*10., xy_of_lower_left=(-10.745, 0.))

    Put the data back into the new grid.

    >>> _ = hmg.add_field(
    ...     'topographic__elevation',
    ...     hmg.node_x**2 + np.round(hmg.node_y)**2,
    ...     at = 'node'
    ... )
    >>> fa = FlowAccumulator(
    ...      hmg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ...      )
    >>> fa.run_one_step()
    >>> hmg.at_node['surface_water__discharge']
    array([ 500.,    0.,    0.,
            200.,  500.,  200.,    0.,
              0.,  200.,  200.,  100.,    0.,
              0.,  100.,  100.,    0.,
              0.,    0.,    0.])

    Next, let's see what happens to a raster grid when there is a depression.

    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)

    This model grid has a depression in the center.

    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)
    array([[ 0.    ,  0.5   ,  1.    ,  1.5   ,  2.    ,  2.5   ,  3.    ],
           [ 0.005 ,  0.505 ,  1.005 ,  1.505 ,  2.005 ,  2.505 ,  3.005 ],
           [ 0.01  ,  0.51  ,  0.101 ,  0.151 ,  0.201 ,  2.51  ,  3.01  ],
           [ 0.015 ,  0.515 ,  0.1015,  0.1515,  0.2015,  2.515 ,  3.015 ],
           [ 0.02  ,  0.52  ,  0.102 ,  0.152 ,  0.202 ,  2.52  ,  3.02  ],
           [ 0.025 ,  0.525 ,  1.025 ,  1.525 ,  2.025 ,  2.525 ,  3.025 ],
           [ 0.03  ,  0.53  ,  1.03  ,  1.53  ,  2.03  ,  2.53  ,  3.03  ]])
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest
    ...      )
    >>> fa.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg.at_node['flow__receiver_node'].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14, 16, 16, 17, 18, 20],
           [21, 21, 16, 23, 24, 25, 27],
           [28, 28, 23, 30, 31, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
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

    >>> df_4 = DepressionFinderAndRouter(mg)
    >>> df_4.map_depressions()
    >>> mg.at_node['flow__receiver_node'].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
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

    >>> df_4.lake_at_node.reshape(mg.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]], dtype=bool)
    >>> df_4.lake_map.reshape(mg.shape)  # doctest: +NORMALIZE_WHITESPACE
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

    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director='FlowDirectorD8',
    ...      depression_finder=DepressionFinderAndRouter
    ... )
    >>> fa.run_one_step()

    This has the same effect of first calling the accumulator and then calling
    the depression finder.

    >>> mg.at_node['flow__receiver_node'].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 5.25,  5.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.25,  0.25,  5.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.75,  2.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.5 ,  0.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])

    The depression finder is stored as part of the flow accumulator, so its
    properties can be accessed through the depression finder.

    >>> fa.depression_finder.lake_at_node.reshape(mg.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]], dtype=bool)
    >>> fa.depression_finder.lake_map.reshape(mg.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1]])
    >>> fa.depression_finder.lake_codes  # a unique code for each lake present on the grid
    array([16])
    >>> fa.depression_finder.lake_outlets  # the outlet node of each lake in lake_codes
    array([8])
    >>> fa.depression_finder.lake_areas  # the area of each lake in lake_codes
    array([ 2.25])

    Finally, note that the DepressionFinderAndRouter takes a keyword argument
    *routing* ('D8', default; 'D4') that sets how connectivity is set between
    nodes. Similar to our ability to pass keyword arguments to the FlowDirector
    through FlowAccumulator, we can pass this keyword argument to the
    DepressionFinderAndRouter component.

    >>> fa = FlowAccumulator(
    ...      mg,
    ...      'topographic__elevation',
    ...      flow_director=FlowDirectorSteepest,
    ...      depression_finder=DepressionFinderAndRouter,
    ...      routing='D4'
    ... )

    FlowAccumulator was designed to work with all types of grids. However,
    NetworkModelGrid's have no cell area. Thus, in order for FlowAccumulator to
    this type of grid, an at-node array called ``cell_area_at_node`` must be
    present.

    >>> from landlab.grid.network import NetworkModelGrid
    >>> y_of_node = (0, 1, 2, 2)
    >>> x_of_node = (0, 0, -1, 1)
    >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
    >>> nmg = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> area = nmg.add_ones('node', 'cell_area_at_node')
    >>> z = nmg.add_field(
    ...     'topographic__elevation',
    ...     nmg.x_of_node + nmg.y_of_node,
    ...     at = 'node')
    >>> fa = FlowAccumulator(nmg)
    >>> fa.run_one_step()
    >>> nmg.at_node['flow__receiver_node']
    array([0, 0, 2, 1])
    """

    _name = "FlowAccumulator"

    _input_var_names = ("topographic__elevation", "water__unit_flux_in")

    _output_var_names = (
        "drainage_area",
        "surface_water__discharge",
        "flow__upstream_node_order",
        # "flow__nodes_not_in_stack",
        "flow__data_structure_delta",
        "flow__data_structure_D",
    )

    _var_units = {
        "topographic__elevation": "m",
        "flow__receiver_node": "m",
        "water__unit_flux_in": "m/s",
        "drainage_area": "m**2",
        "surface_water__discharge": "m**3/s",
        "flow__upstream_node_order": "-",
        "flow__data_structure_delta": "-",
        "flow__data_structure_D": "-",
        # "flow__nodes_not_in_stack": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "water__unit_flux_in": "node",
        "drainage_area": "node",
        "surface_water__discharge": "node",
        "flow__upstream_node_order": "node",
        # "flow__nodes_not_in_stack": "grid",
        "flow__data_structure_delta": "node",
        "flow__data_structure_D": "grid",
    }
    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
        "drainage_area": "Upstream accumulated surface area contributing to the node's "
        "discharge",
        "surface_water__discharge": "Discharge of water through each node",
        "water__unit_flux_in": "External volume water per area per time input to each node "
        "(e.g., rainfall rate)",
        "flow__upstream_node_order": "Node array containing downstream-to-upstream ordered list of "
        "node IDs",
        "flow__data_structure_delta": "Node array containing the elements delta[1:] of the data "
        'structure "delta" used for construction of the downstream-to-'
        "upstream node array",
        "flow__data_structure_D": "Array containing the data structure D used for construction"
        "of the downstream-to-upstream node array. Stored at Grid.",
        # "flow__nodes_not_in_stack": "Boolean value indicating if there are any nodes that have not yet"
        # "been added to the stack stored in flow__upstream_node_order.",
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        flow_director="FlowDirectorSteepest",
        runoff_rate=None,
        depression_finder=None,
        **kwargs
    ):
        """Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and
        compatability for the flow_director and depression_finder
        keyword arguments, tests the argument of runoff_rate, and
        initializes new fields.
        """
        super(FlowAccumulator, self).__init__(grid)
        # Keep a local reference to the grid
        self._grid = grid

        # Grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        self._is_Network = isinstance(self._grid, NetworkModelGrid)
        self.kwargs = kwargs
        # STEP 1: Testing of input values, supplied either in function call or
        # as part of the grid.
        self._test_water_inputs(grid, runoff_rate)

        # save elevations and node_cell_area to class properites.
        self.surface = surface
        self.surface_values = return_array_at_node(grid, surface)

        if self._is_Network:
            try:
                node_cell_area = self._grid.at_node["cell_area_at_node"]
            except FieldError:
                raise FieldError(
                    "In order for the FlowAccumulator to work, the "
                    "grid must have an at-node field called "
                    "cell_area_at_node."
                )
        else:
            node_cell_area = self._grid.cell_area_at_node.copy()
            node_cell_area[self._grid.closed_boundary_nodes] = 0.0

        self.node_cell_area = node_cell_area

        # STEP 2:
        # identify Flow Director method, save name, import and initialize the correct
        # flow director component if necessary
        self._add_director(flow_director)
        self._add_depression_finder(depression_finder)

        # This component will track of the following variables.
        # Attempt to create each, if they already exist, assign the existing
        # version to the local copy.

        #   - drainage area at each node
        #   - receiver of each node
        #   - D array
        #   - delta array
        #   - missing nodes in stack.
        if "drainage_area" not in grid.at_node:
            self.drainage_area = grid.add_zeros("drainage_area", at="node", dtype=float)
        else:
            self.drainage_area = grid.at_node["drainage_area"]

        if "surface_water__discharge" not in grid.at_node:
            self.discharges = grid.add_zeros(
                "surface_water__discharge", at="node", dtype=float
            )
        else:
            self.discharges = grid.at_node["surface_water__discharge"]

        if "flow__upstream_node_order" not in grid.at_node:
            self.upstream_ordered_nodes = grid.add_field(
                "flow__upstream_node_order",
                BAD_INDEX_VALUE * grid.ones(at="node", dtype=int),
                at="node",
                dtype=int,
            )
        else:
            self.upstream_ordered_nodes = grid.at_node["flow__upstream_node_order"]

        if "flow__data_structure_delta" not in grid.at_node:
            self.delta_structure = grid.add_field(
                "flow__data_structure_delta",
                BAD_INDEX_VALUE * grid.ones(at="node", dtype=int),
                at="node",
                dtype=int,
            )
        else:
            self.delta_structure = grid.at_node["flow__data_structure_delta"]

        try:
            D = BAD_INDEX_VALUE * grid.ones(at="link", dtype=int)
            D_structure = np.array([D], dtype=object)
            self.D_structure = grid.add_field(
                "flow__data_structure_D",
                D_structure,
                at="grid",
                dtype=object,
                noclobber=False,
            )

        except FieldError:
            self.D_structure = grid.at_grid["flow__data_structure_D"]

        self.nodes_not_in_stack = True

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid["node"]["drainage_area"]

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid["node"]["surface_water__discharge"]

    @property
    def node_order_upstream(self):
        """Return the upstream node order (drainage stack)."""
        return self._grid["node"]["flow__upstream_node_order"]

    def link_order_upstream(self):
        """Return the upstream order of active links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> mg = RasterModelGrid((5, 5))
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> _ = mg.add_field(
        ...     'topographic__elevation',
        ...     mg.node_x + mg.node_y,
        ...     at = 'node'
        ...     )
        >>> fa = FlowAccumulator(mg, 'topographic__elevation')
        >>> fa.run_one_step()
        >>> fa.link_order_upstream()
        array([ 5, 14, 23,  6, 15, 24,  7, 16, 25])

        This also works for route-to-many methods

        >>> mg = RasterModelGrid((5, 5))
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> _ = mg.add_field(
        ...     'topographic__elevation',
        ...     mg.node_x + mg.node_y,
        ...     at = 'node'
        ... )
        >>> fa = FlowAccumulator(mg,
        ...      'topographic__elevation',
        ...      flow_director='MFD')
        >>> fa.run_one_step()
        >>> fa.link_order_upstream()
        array([ 5, 14, 10,  6, 11,  7, 23, 19, 15, 20, 16, 28, 24, 29, 25])
        """
        downstream_links = self._grid["node"]["flow__link_to_receiver_node"][
            self.node_order_upstream
        ]
        out = downstream_links.flatten()
        return out[out != BAD_INDEX_VALUE]

    def headwater_nodes(self):
        """Return the headwater nodes.

        These are nodes that contribute flow and have no upstream nodes.

        Examples
        --------
        >>> from numpy.testing import assert_array_equal
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> mg = RasterModelGrid((5, 5))
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> _ = mg.add_field(
        ...     'topographic__elevation',
        ...     mg.node_x + mg.node_y,
        ...     at = 'node'
        ... )
        >>> fa = FlowAccumulator(mg, 'topographic__elevation')
        >>> fa.run_one_step()
        >>> assert_array_equal(fa.headwater_nodes(), np.array([16, 17, 18]))
        """
        delta = np.concatenate(([0], self.delta_structure))
        num_donors = np.diff(delta)
        # note closed nodes have a value of 1 here since they flow to
        # themselves
        source_nodes = np.where(num_donors == 0)[0]
        return source_nodes

    def _test_water_inputs(self, grid, runoff_rate):
        """Test inputs for runoff_rate and water__unit_flux_in."""
        if "water__unit_flux_in" not in grid.at_node:
            if runoff_rate is None:
                # assume that if runoff rate is not supplied, that the value
                # should be set to one everywhere.
                grid.add_ones("node", "water__unit_flux_in", dtype=float)
            else:
                runoff_rate = return_array_at_node(grid, runoff_rate)
                grid.at_node["water__unit_flux_in"] = runoff_rate
        else:
            if runoff_rate is not None:
                print(
                    "FlowAccumulator found both the field "
                    + "'water__unit_flux_in' and a provided float or "
                    + "array for the runoff_rate argument. THE FIELD IS "
                    + "BEING OVERWRITTEN WITH THE SUPPLIED RUNOFF_RATE!"
                )
                runoff_rate = return_array_at_node(grid, runoff_rate)
                grid.at_node["water__unit_flux_in"] = runoff_rate

        # perform a test (for politeness!) that the old name for the water_in
        # field is not present:
        if "water__discharge_in" in grid.at_node:
            warnings.warn(
                "This component formerly took 'water__discharge"
                + "_in' as an input field. However, this field is "
                + "now named 'water__unit_flux_in'. You are still "
                + "using a field with the old name. Please update "
                + "your code if you intended to use that field.",
                DeprecationWarning,
            )

    def _add_director(self, flow_director):
        """Test and add the flow director component."""
        PERMITTED_DIRECTORS = [
            "FlowDirectorSteepest",
            "FlowDirectorD8",
            "FlowDirectorMFD",
            "FlowDirectorDINF",
        ]

        potential_kwargs = ["partition_method", "diagonals"]
        kw = {}
        for p_k in potential_kwargs:
            if p_k in self.kwargs.keys():
                kw[p_k] = self.kwargs.pop(p_k)

        # flow director is provided as a string.
        if isinstance(flow_director, six.string_types):
            if flow_director[:12] == "FlowDirector":
                flow_director = flow_director[12:]

            from landlab.components.flow_director import (
                FlowDirectorSteepest,
                FlowDirectorD8,
                FlowDirectorMFD,
                FlowDirectorDINF,
            )

            DIRECTOR_METHODS = {
                "D4": FlowDirectorSteepest,
                "Steepest": FlowDirectorSteepest,
                "D8": FlowDirectorD8,
                "MFD": FlowDirectorMFD,
                "DINF": FlowDirectorDINF,
            }

            try:
                FlowDirector = DIRECTOR_METHODS[flow_director]
            except KeyError:
                raise ValueError(
                    "String provided in flow_director is not a "
                    "valid method or component name. The following"
                    "components are valid imputs:\n" + str(PERMITTED_DIRECTORS)
                )
            self.flow_director = FlowDirector(self._grid, self.surface, **kw)
        # flow director is provided as an instantiated flow director
        elif isinstance(flow_director, Component):
            if flow_director._name in PERMITTED_DIRECTORS:
                self.flow_director = flow_director
            else:
                raise ValueError(
                    "String provided in flow_director is not a "
                    "valid method or component name. The following"
                    "components are valid imputs:\n" + str(PERMITTED_DIRECTORS)
                )

            if len(kw) > 0:
                raise ValueError(
                    "flow_director provided as an instantiated ",
                    "component and keyword arguments provided. ",
                    "These kwargs would be ignored.",
                )

        # flow director is provided as an uninstantiated flow director
        else:
            if flow_director._name in PERMITTED_DIRECTORS:
                FlowDirector = flow_director
                self.flow_director = FlowDirector(self._grid, self.surface, **kw)
            else:
                raise ValueError(
                    "String provided in flow_director is not a "
                    "valid method or component name. The following"
                    "components are valid imputs:\n" + str(PERMITTED_DIRECTORS)
                )

        # save method as attribute
        self.method = self.flow_director.method

    def _add_depression_finder(self, depression_finder):
        """Test and add the depression finder component."""
        PERMITTED_DEPRESSION_FINDERS = ["DepressionFinderAndRouter"]

        # now do a similar thing for the depression finder.
        self.depression_finder_provided = depression_finder
        if self.depression_finder_provided is not None:

            # collect potential kwargs to pass to depression_finder
            # instantiation
            potential_kwargs = ["routing"]
            kw = {}
            for p_k in potential_kwargs:
                if p_k in self.kwargs.keys():
                    kw[p_k] = self.kwargs.pop(p_k)

            # NEED TO TEST WHICH FLOWDIRECTOR WAS PROVIDED.
            if self.flow_director._name in ("FlowDirectorMFD", "FlowDirectorDINF"):
                msg = (
                    "The depression finder only works with route "
                    "to one FlowDirectors such as "
                    "FlowDirectorSteepest and  FlowDirectorD8. "
                    "Provide a different FlowDirector."
                )
                raise NotImplementedError(msg)

            # if D4 is being used here and should be.
            if (
                (("routing" not in kw) or (kw["routing"] != "D4"))
                and isinstance(self._grid, RasterModelGrid)
                and (self.flow_director._name in ("FlowDirectorSteepest"))
            ):

                message = (
                    "You have specified \n"
                    "flow_director=FlowDirectorSteepest and\n"
                    "depression_finder=DepressionFinderAndRouter\n"
                    "in the instantiation of FlowAccumulator on a "
                    "RasterModelGrid. The default behavior of "
                    "DepressionFinderAndRouter is to use D8 connectivity "
                    "which is in conflict with D4 connectivity used by "
                    "FlowDirectorSteepest. \n"
                    "To fix this, provide the kwarg routing='D4', when "
                    "you instantiate FlowAccumulator."
                )

                raise ValueError(warning_message(message))

            # depression finder is provided as a string.
            if isinstance(self.depression_finder_provided, six.string_types):

                from landlab.components import DepressionFinderAndRouter

                DEPRESSION_METHODS = {
                    "DepressionFinderAndRouter": DepressionFinderAndRouter
                }

                try:
                    DepressionFinder = DEPRESSION_METHODS[
                        self.depression_finder_provided
                    ]
                except KeyError:
                    raise ValueError(
                        "Component provided in depression_finder "
                        "is not a valid component. The following "
                        "components are valid imputs:\n"
                        + str(PERMITTED_DEPRESSION_FINDERS)
                    )

                self.depression_finder = DepressionFinder(self._grid, **kw)
            # flow director is provided as an instantiated depression finder
            elif isinstance(self.depression_finder_provided, Component):

                if (
                    self.depression_finder_provided._name
                    in PERMITTED_DEPRESSION_FINDERS
                ):
                    self.depression_finder = self.depression_finder_provided
                else:
                    raise ValueError(
                        "Component provided in depression_finder "
                        "is not a valid component. The following "
                        "components are valid imputs:\n"
                        + str(PERMITTED_DEPRESSION_FINDERS)
                    )

                if len(kw) > 0:
                    raise ValueError(
                        "flow_director provided as an instantiated ",
                        "component and keyword arguments provided. ",
                        "These kwargs would be ignored.",
                    )

            # depression_fiuner is provided as an uninstantiated depression finder
            else:

                if (
                    self.depression_finder_provided._name
                    in PERMITTED_DEPRESSION_FINDERS
                ):
                    DepressionFinder = self.depression_finder_provided
                    self.depression_finder = DepressionFinder(self._grid, **kw)
                else:
                    raise ValueError(
                        "Component provided in depression_finder "
                        "is not a valid component. The following "
                        "components are valid imputs:\n"
                        + str(PERMITTED_DEPRESSION_FINDERS)
                    )
        else:
            self.depression_finder = None

    def accumulate_flow(self, update_flow_director=True):
        """Function to make FlowAccumulator calculate drainage area and
        discharge.

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
        # set a couple of aliases
        a = self._grid["node"]["drainage_area"]
        q = self._grid["node"]["surface_water__discharge"]

        # step 1. Find flow directions by specified method
        if update_flow_director:
            self.flow_director.run_one_step()

        # further steps vary depending on how many recievers are present
        # one set of steps is for route to one (D8, Steepest/D4)

        # step 2. Get r
        r = as_id_array(self._grid["node"]["flow__receiver_node"])

        if self.flow_director.to_n_receivers == "one":

            # step 2b. Run depression finder if passed
            # Depression finder reaccumulates flow at the end of its routine.
            # At the moment, no depression finders work with to-many, so it
            # lives here
            if self.depression_finder_provided is not None:
                self.depression_finder.map_depressions()

                # if FlowDirectorSteepest is used, update the link directions
                if self.flow_director._name == "FlowDirectorSteepest":
                    self.flow_director._determine_link_directions()

            # step 3. Stack, D, delta construction
            nd = as_id_array(flow_accum_bw._make_number_of_donors_array(r))
            delta = as_id_array(flow_accum_bw._make_delta_array(nd))
            D = as_id_array(flow_accum_bw._make_array_of_donors(r, delta))
            s = as_id_array(flow_accum_bw.make_ordered_node_array(r))

            # put these in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid["node"]["flow__data_structure_delta"][:] = delta[1:]
            self._grid["grid"]["flow__data_structure_D"] = np.array([D], dtype=object)
            self._grid["node"]["flow__upstream_node_order"][:] = s

            # step 4. Accumulate (to one or to N depending on direction method)
            a[:], q[:] = self._accumulate_A_Q_to_one(s, r)

        else:
            # Get p
            p = self._grid["node"]["flow__receiver_proportions"]

            # step 3. Stack, D, delta construction
            nd = as_id_array(flow_accum_to_n._make_number_of_donors_array_to_n(r, p))
            delta = as_id_array(flow_accum_to_n._make_delta_array_to_n(nd))
            D = as_id_array(flow_accum_to_n._make_array_of_donors_to_n(r, p, delta))
            s = as_id_array(flow_accum_to_n.make_ordered_node_array_to_n(r, p))

            # put theese in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid["node"]["flow__data_structure_delta"][:] = delta[1:]
            self._grid["grid"]["flow__data_structure_D"][0] = np.array(
                [D], dtype=object
            )
            self._grid["node"]["flow__upstream_node_order"][:] = s
            self._grid["node"]["flow__upstream_node_order"][:] = s

            # step 4. Accumulate (to one or to N depending on direction method)
            a[:], q[:] = self._accumulate_A_Q_to_n(s, r, p)

        return (a, q)

    def _accumulate_A_Q_to_one(self, s, r):
        """Accumulate area and discharge for a route-to-one scheme.

        Note this can be overridden in inherited components.
        """
        a, q = flow_accum_bw.find_drainage_area_and_discharge(
            s, r, self.node_cell_area, self._grid.at_node["water__unit_flux_in"]
        )
        return (a, q)

    def _accumulate_A_Q_to_n(self, s, r, p):
        """Accumulate area and discharge for a route-to-many scheme.

        Note this can be overridden in inherited components.
        """
        a, q = flow_accum_to_n.find_drainage_area_and_discharge_to_n(
            s, r, p, self.node_cell_area, self._grid.at_node["water__unit_flux_in"]
        )
        return (a, q)

    def run_one_step(self):
        """Accumulate flow and save to the model grid.

        run_one_step() checks for updated boundary conditions, calculates
        slopes on links, finds baselevel nodes based on the status at node,
        calculates flow directions, and accumulates flow and saves results to
        the grid.

        An alternative to run_one_step() is accumulate_flow() which does the
        same things but also returns the drainage area and discharge.
        """
        self.accumulate_flow()


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
