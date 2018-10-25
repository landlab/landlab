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
from landlab import RasterModelGrid, VoronoiDelaunayGrid, ModelGrid
# ^^for type tests
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message

from landlab.components.flow_accum import flow_accum_bw
from landlab.components.flow_accum import flow_accum_to_n
from landlab.components.flow_accum import FlowAccumulator

from landlab import BAD_INDEX_VALUE
import six
import numpy as np


class LossyFlowAccumulator(FlowAccumulator):

    """
    Component to calculate drainage area and accumulate flow, while permitting
    dynamic loss of flow downstream.

    This component is closely related to the FlowAccumulator, in that
    this is accomplished by first finding flow directions by a user-specified
    method and then calculating the drainage area and discharge. However,
    this component additionally requires the passing of a function that
    describes how discharge is lost downstream, f(Qw, nodeID, linkID).

    Optionally, spatially variable runoff can be set either by the model grid
    field 'water__unit_flux_in' or the input variable *runoff_rate**.

    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished
    together. Note that the DepressionFinderAndRouter is not particularly
    intelligent when running on lossy streams, and in particular, it will
    reroute flow around pits even when they are in fact not filled due to loss.

    NOTE: The perimeter nodes  NEVER contribute to the accumulating flux, even
    if the  gradients from them point inwards to the main body of the grid.
    This is because under Landlab definitions, perimeter nodes lack cells, so
    cannot accumulate any discharge.

    LossyFlowAccumulator stores as ModelGrid fields:

        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure:
            *flow__data_structure_delta*. The first element is always zero.
        -  Link array of the D data structure: *flow__data_structure_D*

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
    runoff_rate : field name, array, or float, optional (m/time)
        If provided, sets the runoff rate and will be assigned to the grid field
        'water__unit_flux_in'. If a spatially and and temporally variable runoff
        rate is desired, pass this field name and update the field through model
        run time. If both the field and argument are present at the time of
        initialization, runoff_rate will *overwrite* the field. If neither are
        set, defaults to spatially constant unit input.
    depression_finder : string, class, instance of class, optional
         A string of class name (e.g., 'DepressionFinderAndRouter'), an
         uninstantiated DepressionFinder class, or an instance of a
         DepressionFinder class.
         This sets the method for depression finding.
    loss_function : Python function, optional
        A function of the form f(Qw, [node_ID, [linkID, [grid]]]), where Qw is
        the discharge at a node, node_ID the ID of the node at which the loss
        is to be calculated, linkID is the ID of the link down which the
        outflow drains, and grid is a Landlab ModelGrid. Note that if a linkID
        is needed, a nodeID must also be specified, even if only as a dummy
        parameter; similarly, if a grid is to be passed, all of the preceding
        parameters must be specified. Both nodeID and linkID are required to
        permit spatially variable losses, and also losses dependent on flow
        path geometry (e.g., flow length). The grid is passed to allow fields
        or grid properties describing values across the grid to be accessed
        for the loss calculation (see examples).
        This function should take (float, [int, [int, [ModelGrid]]]), and
        return a single float.
    **kwargs : any additional parameters to pass to a FlowDirector or
         DepressionFinderAndRouter instance (e.g., partion_method for
         FlowDirectorMFD). This will have no effect if an instantiated component
         is passed using the flow_director or depression_finder keywords.

    Examples
    --------
    (Note that examples of the use of the loss function specifically are at the
    end of these examples.)

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import LossyFlowAccumulator
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at = 'node')

    The LossyFlowAccumulator component accumulates flow and calculates drainage
    using all of the different methods for directing flow in Landlab. These
    include steepest descent (also known as D4 for the case of a raster grid)
    and D8 (on raster grids only). The method for flow director can be
    specified as a string (e.g., 'D8' or 'FlowDirectorD8'), as an
    uninstantiated FlowDirector component or as an instantiated FlowDirector
    component.

    The default method is to use FlowDirectorSteepest.

    First let's look at the three ways to instantiate a LossyFlowAccumulator.
    The following four methods are all equivalent. First, we can pass the
    entire name of a flow director as a string to the argument `flow_director`:

    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director='FlowDirectorSteepest')

    Second, we can pass just the method name as a string to the argument
    `flow_director`:

    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director='Steepest')

    Third, we can import a FlowDirector component from Landlab and pass it to
    `flow_director`:
    >>> from landlab.components import FlowDirectorSteepest
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)

    Finally, we can instantiate a FlowDirector component and pass this
    instantiated version to `flow_director`. You might want to do this if you
    used a FlowDirector in order to set up something before starting a
    time loop and then want to use the same flow director within the loop.

    >>> fd = FlowDirectorSteepest(mg, 'topographic__elevation')
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)

    Now let's look at what LossyFlowAccumulator does. Even before we run
    LossyFlowAccumulator it has the property `surface_values` that stores the
    values of the surface over which flow is directed and accumulated.

    >>> fa.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])

    Now let's make a more complicated elevation grid for the next examples.

    >>> mg = RasterModelGrid((5, 4), spacing=(1, 1))
    >>> topographic__elevation = np.array([0.,  0.,  0., 0.,
    ...                                    0., 21., 10., 0.,
    ...                                    0., 31., 20., 0.,
    ...                                    0., 32., 30., 0.,
    ...                                    0.,  0.,  0., 0.])
    >>> _ = mg.add_field('node',
    ...                    'topographic__elevation',
    ...                    topographic__elevation)
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)
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

    >>> mg = RasterModelGrid((5, 4), spacing=(10., 10))

    Put the data back into the new grid.

    >>> _ = mg.add_field('node',
    ...                    'topographic__elevation',
    ...                    topographic__elevation)
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)
    >>> runoff_rate = np.arange(mg.number_of_nodes)
    >>> _ = mg.add_field('node', 'water__unit_flux_in', runoff_rate,
    ...                  noclobber=False)
    >>> fa.run_one_step()
    >>> mg.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
    array([    0.,   500.,  5200.,     0.,
               0.,   500.,  5200.,     0.,
               0.,   900.,  4600.,     0.,
               0.,  1300.,  2700.,     0.,
               0.,     0.,     0.,     0.])

    The LossyFlowAccumulator component will work for both raster grids and
    irregular grids. For the example we will use a Hexagonal Model Grid, a
    special type of Voroni Grid that has regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid(5,3)
    >>> _ = hmg.add_field('topographic__elevation',
    ...                   hmg.node_x + np.round(hmg.node_y),
    ...                   at = 'node')
    >>> fa = LossyFlowAccumulator(hmg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)
    >>> fa.surface_values
    array([ 0. ,  1. ,  2. ,
            0.5,  1.5,  2.5,  3.5,
            1. ,  2. ,  3. ,  4. ,  5. ,
            2.5,  3.5,  4.5,  5.5,
            3. ,  4. ,  5. ])

    If the FlowDirector you want to use takes keyword arguments and you want
    to specify it using a string or uninstantiated FlowDirector class, include
    those keyword arguments when you create LossyFlowAccumulator.

    For example, in the case of a raster grid, FlowDirectorMFD can use only
    orthogonal links, or it can use both orthogonal and diagonal links.

    >>> mg = RasterModelGrid((5, 5), spacing=(1, 1))
    >>> topographic__elevation = mg.node_y+mg.node_x
    >>> _ = mg.add_field('node',
    ...                  'topographic__elevation',
    ...                   topographic__elevation)
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director='MFD',
    ...                           diagonals = True)
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
    >>> mg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 1.41168825,  2.06497116,  1.3253788 ,  0.40380592,  0.        ,
            2.06497116,  3.40811691,  2.5753788 ,  1.37867966,  0.        ,
            1.3253788 ,  2.5753788 ,  2.17157288,  1.29289322,  0.        ,
            0.40380592,  1.37867966,  1.29289322,  1.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ])

    It may seem odd that there are no round numbers in the drainage area field.
    This is because flow is directed to all downhill boundary nodes and
    partitioned based on slope.

    To check that flow is conserved, sum along all boundary nodes.

    >>> sum(mg.at_node['drainage_area'][mg.boundary_nodes])
    9.0000000000000018

    This should be the same as the number of core nodes --- as boundary nodes
    in landlab do not have area.

    >>> len(mg.core_nodes)
    9

    Next, let's set the dx spacing such that each cell has an area of one.

    >>> dx=(2./(3.**0.5))**0.5
    >>> hmg = HexModelGrid(5,3, dx)
    >>> _ = hmg.add_field('topographic__elevation',
    ...                     hmg.node_x**2 + np.round(hmg.node_y)**2,
    ...                     at = 'node')
    >>> fa = LossyFlowAccumulator(hmg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)
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

    >>> hmg = HexModelGrid(5,3, dx*10.)

    Put the data back into the new grid.

    >>> _ = hmg.add_field('topographic__elevation',
    ...                     hmg.node_x**2 + np.round(hmg.node_y)**2,
    ...                     at = 'node')
    >>> fa = LossyFlowAccumulator(hmg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)
    >>> fa.run_one_step()
    >>> hmg.at_node['surface_water__discharge']
    array([ 500.,    0.,    0.,
            200.,  500.,  200.,    0.,
              0.,  200.,  200.,  100.,    0.,
              0.,  100.,  100.,    0.,
              0.,    0.,    0.])

    Next, let's see what happens to a raster grid when there is a depression.

    >>> mg = RasterModelGrid((7, 7), 0.5)
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
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest)
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

    >>> mg = RasterModelGrid((7, 7), 0.5)
    >>> z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director='FlowDirectorD8',
    ...                           depression_finder=DepressionFinderAndRouter)
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
    through LossyFlowAccumulator, we can pass this keyword argument to the
    DepressionFinderAndRouter component.

    Now, some examples of the loss happening. Here's an example of a 50% loss
    of flow every time flow moves along a node:

    >>> mg = RasterModelGrid((3, 5), (1, 2))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    >>> z = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at='node')

    >>> def mylossfunction(qw):
    ...     return 0.5 * qw

    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest,
    ...                           routing='D4', loss_function=mylossfunction)
    >>> fa.run_one_step()

    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 6.,  6.,  4.,  2.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    >>> mg.at_node['surface_water__discharge'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 1.75,  3.5 ,  3.  ,  2.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])

    >>> def mylossfunction2(Qw, nodeID):
    ...     return
    """

    _name = "LossyFlowAccumulator"

    _input_var_names = ("topographic__elevation", "water__unit_flux_in")

    _output_var_names = (
        "drainage_area",
        "surface_water__discharge",
        "flow__upstream_node_order",
        "flow__nodes_not_in_stack",
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
        "flow__nodes_not_in_stack": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "water__unit_flux_in": "node",
        "drainage_area": "node",
        "surface_water__discharge": "node",
        "flow__upstream_node_order": "node",
        "flow__nodes_not_in_stack": "grid",
        "flow__data_structure_delta": "node",
        "flow__data_structure_D": "link",
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
        "flow__data_structure_D": "Link array containing the data structure D used for construction"
        "of the downstream-to-upstream node array",
        "flow__nodes_not_in_stack": "Boolean value indicating if there are any nodes that have not yet"
        "been added to the stack stored in flow__upstream_node_order.",
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        flow_director="FlowDirectorSteepest",
        runoff_rate=None,
        depression_finder=None,
        loss_function=None,
        **kwargs
    ):
        """
        Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and compatability
        for the flow_director and depression_finder keyword arguments, tests
        the argument of runoff_rate, and initializes new fields.
        """
        super(LossyFlowAccumulator, self).__init__(
            grid, surface=surface, flow_director=flow_director,
            runoff_rate=runoff_rate, depression_finder=depression_finder,
            **kwargs)

        if loss_function is not None:
            # save the func for loss, and do a quick test on its inputs:
            if loss_function.func_code.co_argcount == 1:
                # check the func takes a single value and turns it into a new
                # single value:
                if not isinstance(loss_function(1.), float):
                    raise TypeError(
                        'The loss_function should take a float, and return a ' +
                        'float.')
                # now, for logical consistency in our calls to
                # find_drainage_area_and_discharge, wrap the func so it has two
                # arguments:
                def lossfunc(Qw, dummyn, dummyl, dummygrid):
                    return float(loss_function(Qw))
                self._lossfunc = lossfunc

            elif loss_function.func_code.co_argcount == 2:
                # check the func takes a single value and turns it into a new
                # single value:
                if not isinstance(loss_function(1., 0), float):
                    raise TypeError(
                        'The loss_function should take (float, int), and ' +
                        'return a float.')
                # now, for logical consistency in our calls to
                # find_drainage_area_and_discharge, wrap the func so it has two
                # arguments:
                def lossfunc(Qw, nodeID, dummyl, dummygrid):
                    return float(loss_function(Qw, nodeID))
                self._lossfunc = lossfunc

            elif loss_function.func_code.co_argcount == 3:
                # check the func takes (float, int) and turns it into a new
                # single value:
                if not isinstance(loss_function(1., 0, 0), float):
                    raise TypeError(
                        'The loss_function should take (float, int, int), ' +
                        'and return a float.')
                def lossfunc(Qw, nodeID, linkID, dummygrid):
                    return float(loss_function(Qw, nodeID, linkID))
                self._lossfunc = lossfunc

            elif loss_function.func_code.co_argcount == 4:
                # this time, the test is too hard to implement cleanly so just
                self._lossfunc = loss_function
            else:
                raise ValueError(
                    'The loss_function must have only a single argument, ' +
                    'which should be the discharge at a node; a pair of ' +
                    'arguments, which should be the discharge at a node and ' +
                    'the node ID; or three arguments, which should be the ' +
                    'discharge at a node, the node ID, and the link along ' +
                    'which that discharge will flow.')
        else:
            # make a dummy
            def lossfunc(Qw, dummyn, dummyl, dummygrid):
                return float(Qw)
            self._lossfunc = lossfunc

    def accumulate_flow(self, update_flow_director=True):
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
        if self.flow_director.to_n_receivers == "one":

            # step 3. Run depression finder if passed
            # Depression finder reaccumulates flow at the end of its routine.
            if self.depression_finder_provided is not None:
                # prevent internal flow rerouting (which ignores loss), and
                # do it (more slowly) here instead
                self.depression_finder.map_depressions()

            # step 2. Get r
            r = self._grid["node"]["flow__receiver_node"]

            # step 2. Stack, D, delta construction
            nd = flow_accum_bw._make_number_of_donors_array(r)
            delta = flow_accum_bw._make_delta_array(nd)
            D = flow_accum_bw._make_array_of_donors(r, delta)
            s = flow_accum_bw.make_ordered_node_array(r)
            link = self._grid.at_node['flow__link_to_receiver_node']

            # put theese in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid["node"]["flow__data_structure_delta"][:] = delta[1:]
            self._grid["link"]["flow__data_structure_D"][: len(D)] = D
            self._grid["node"]["flow__upstream_node_order"][:] = s

            # step 4. Accumulate (to one or to N depending on direction
            # method. )
            a, q = flow_accum_bw.find_drainage_area_and_discharge_lossy(
                s, r, link, self._lossfunc, self._grid, self.node_cell_area,
                self._grid.at_node["water__unit_flux_in"]
            )
            self._grid["node"]["drainage_area"][:] = a
            self._grid["node"]["surface_water__discharge"][:] = q

        else:
            # step 2. Get r and p
            r = self._grid["node"]["flow__receiver_node"]
            p = self._grid["node"]["flow__receiver_proportions"]
            link = self._grid.at_node['flow__link_to_receiver_node']

            # step 2. Stack, D, delta construction
            nd = flow_accum_to_n._make_number_of_donors_array_to_n(r, p)
            delta = flow_accum_to_n._make_delta_array_to_n(nd)
            D = flow_accum_to_n._make_array_of_donors_to_n(r, p, delta)
            s = flow_accum_to_n.make_ordered_node_array_to_n(r, p)

            # put theese in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid["node"]["flow__data_structure_delta"][:] = delta[1:]

            if self._is_raster:
                tempD = BAD_INDEX_VALUE * np.ones(
                    (self._grid.number_of_links * 2))
                tempD[: len(D)] = D
                self._grid["link"]["flow__data_structure_D"][
                    :] = tempD.reshape((self._grid.number_of_links, 2))
            else:
                self._grid["link"]["flow__data_structure_D"][: len(D)] = D
            self._grid["node"]["flow__upstream_node_order"][:] = s

            # step 3. Run depression finder if passed
            # at present this must go at the end.

            # step 4. Accumulate (to one or to N depending on dir method. )
            a, q = flow_accum_to_n.find_drainage_area_and_discharge_to_n_lossy(
                s, r, link, p, self._lossfunc, self._grid, self.node_cell_area,
                self._grid.at_node["water__unit_flux_in"]
            )
            # store drainage area and discharge.
            self._grid["node"]["drainage_area"][:] = a
            self._grid["node"]["surface_water__discharge"][:] = q

            # at the moment, this is where the depression finder needs to live.
            # at the moment, no depression finders work with to-many
            # if self.depression_finder_provided is not None:
            #     self.depression_finder.map_depressions()

        return (a, q)


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
