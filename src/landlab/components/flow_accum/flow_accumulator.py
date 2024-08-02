#!/usr/env/python

"""
flow_accumulator.py: Component to accumulate flow and calculate drainage area.

Provides the FlowAccumulator component which accumulates flow and calculates
drainage area. FlowAccumulator supports multiple methods for calculating flow
direction. Optionally a depression finding component can be specified and flow
directing, depression finding, and flow routing can all be accomplished
together.
"""


import numpy as np

from landlab import Component  # for type tests
from landlab import FieldError
from landlab import NetworkModelGrid
from landlab import RasterModelGrid
from landlab import VoronoiDelaunayGrid
from landlab.components.flow_accum import flow_accum_bw
from landlab.components.flow_accum import flow_accum_to_n
from landlab.core.messages import warning_message
from landlab.core.utils import as_id_array
from landlab.utils.return_array import return_array_at_node

from ..depression_finder.floodstatus import FloodStatus


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
       direction: *'flow__link_direction'*

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
    >>> from landlab.grid.raster_funcs import neighbor_to_arrow

    >>> mg = RasterModelGrid((3, 3))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field(
    ...     "topographic__elevation",
    ...     mg.node_x + mg.node_y,
    ...     at="node",
    ... )

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
    ...     mg, "topographic__elevation", flow_director="FlowDirectorSteepest"
    ... )

    Second, we can pass just the method name as a string to the argument
    `flow_director`:

    >>> fa = FlowAccumulator(mg, "topographic__elevation", flow_director="Steepest")

    Third, we can import a FlowDirector component from Landlab and pass it to
    `flow_director`:

    >>> from landlab.components import FlowDirectorSteepest
    >>> fa = FlowAccumulator(
    ...     mg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )

    Finally, we can instantiate a FlowDirector component and pass this
    instantiated version to `flow_director`. You might want to do this if you
    used a FlowDirector in order to set up something before starting a
    time loop and then want to use the same flow director within the loop.

    >>> fd = FlowDirectorSteepest(mg, "topographic__elevation")
    >>> fa = FlowAccumulator(
    ...     mg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )

    Now let's look at what FlowAccumulator does. Even before we run
    FlowAccumulator it has the property `surface_values` that stores the values
    of the surface over which flow is directed and accumulated.

    >>> fa.surface_values
    array([0., 1., 2., 1., 2., 3., 2., 3., 4.])

    Now let's make a more complicated elevation grid for the next examples.

    >>> mg = RasterModelGrid((5, 4))
    >>> topographic__elevation = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 21.0, 10.0, 0.0],
    ...     [0.0, 31.0, 20.0, 0.0],
    ...     [0.0, 32.0, 30.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> _ = mg.add_field("topographic__elevation", topographic__elevation, at="node")
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa = FlowAccumulator(
    ...     mg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )
    >>> fa.run_one_step()
    >>> neighbor_to_arrow(mg.at_node["flow__receiver_node"].reshape(mg.shape))
    array([['○', '○', '○', '○'],
           ['○', '↑', '↑', '○'],
           ['○', '→', '↑', '○'],
           ['○', '→', '↑', '○'],
           ['○', '○', '○', '○']], dtype='<U1')
    >>> mg.at_node["flow__receiver_node"].reshape(mg.shape)
    array([[ 0,  1,  2,  3],
           [ 4,  1,  2,  7],
           [ 8, 10,  6, 11],
           [12, 14, 10, 15],
           [16, 17, 18, 19]])
    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0., 1., 5., 0.],
           [0., 1., 5., 0.],
           [0., 1., 4., 0.],
           [0., 1., 2., 0.],
           [0., 0., 0., 0.]])

    Now let's change the cell area (100.) and the runoff rates:

    >>> mg = RasterModelGrid((5, 4), xy_spacing=(10.0, 10))

    Put the data back into the new grid.

    >>> _ = mg.add_field("topographic__elevation", topographic__elevation, at="node")
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fa = FlowAccumulator(
    ...     mg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )
    >>> runoff_rate = np.arange(mg.number_of_nodes, dtype=float)
    >>> rnff = mg.add_field("water__unit_flux_in", runoff_rate, at="node", clobber=True)
    >>> fa.run_one_step()
    >>> mg.at_node["surface_water__discharge"].reshape(mg.shape)
    array([[   0.,  500.,  5200.,    0.],
           [   0.,  500.,  5200.,    0.],
           [   0.,  900.,  4600.,    0.],
           [   0., 1300.,  2700.,    0.],
           [   0.,    0.,     0.,    0.]])

    The flow accumulator will happily work with a negative runoff rate, which
    could be used to allow, e.g., transmission losses:

    >>> runoff_rate.fill(1.0)
    >>> fa.run_one_step()
    >>> mg.at_node["surface_water__discharge"].reshape(mg.shape)
    array([[  0., 100., 500.,   0.],
           [  0., 100., 500.,   0.],
           [  0., 100., 400.,   0.],
           [  0., 100., 200.,   0.],
           [  0.,   0.,   0.,   0.]])
    >>> runoff_rate[:8] = -0.5
    >>> fa.run_one_step()
    >>> mg.at_node["surface_water__discharge"].reshape(mg.shape)
    array([[  0.,   0., 350.,   0.],
           [  0.,   0., 350.,   0.],
           [  0., 100., 400.,   0.],
           [  0., 100., 200.,   0.],
           [  0.,   0.,   0.,   0.]])

    The drainage area array is unaffected, as you would expect:

    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[  0., 100., 500.,   0.],
           [  0., 100., 500.,   0.],
           [  0., 100., 400.,   0.],
           [  0., 100., 200.,   0.],
           [  0.,   0.,   0.,   0.]])

    The FlowAccumulator component will work for both raster grids and irregular
    grids. For the example we will use a Hexagonal Model Grid, a special type
    of Voroni Grid that has regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid((5, 3), xy_of_lower_left=(-1.0, 0.0))
    >>> _ = hmg.add_field(
    ...     "topographic__elevation",
    ...     hmg.node_x + np.round(hmg.node_y),
    ...     at="node",
    ... )
    >>> fa = FlowAccumulator(
    ...     hmg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )
    >>> fa.surface_values
    array([0. , 1. , 2. ,
           0.5, 1.5, 2.5, 3.5,
           1. , 2. , 3. , 4. , 5. ,
           2.5, 3.5, 4.5, 5.5,
           3. , 4. , 5. ])

    If the FlowDirector you want to use takes keyword arguments and you want
    to specify it using a string or uninstantiated FlowDirector class, include
    those keyword arguments when you create FlowAccumulator.

    For example, in the case of a raster grid, FlowDirectorMFD can use only
    orthogonal links, or it can use both orthogonal and diagonal links.

    >>> mg = RasterModelGrid((5, 5))
    >>> topographic__elevation = mg.node_y + mg.node_x
    >>> _ = mg.add_field("topographic__elevation", topographic__elevation, at="node")
    >>> fa = FlowAccumulator(
    ...     mg, "topographic__elevation", flow_director="MFD", diagonals=True
    ... )
    >>> fa.run_one_step()
    >>> mg.at_node["flow__receiver_node"]
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
    >>> mg.at_node["drainage_area"].round(4).reshape(mg.shape)
    array([[1.4117, 2.065 , 1.3254, 0.4038, 0.    ],
           [2.065 , 3.4081, 2.5754, 1.3787, 0.    ],
           [1.3254, 2.5754, 2.1716, 1.2929, 0.    ],
           [0.4038, 1.3787, 1.2929, 1.    , 0.    ],
           [0.    , 0.    , 0.    , 0.    , 0.    ]])

    It may seem odd that there are no round numbers in the drainage area field.
    This is because flow is directed to all downhill boundary nodes and
    partitioned based on slope.

    To check that flow is conserved, sum along all boundary nodes.

    >>> round(sum(mg.at_node["drainage_area"][mg.boundary_nodes]), 4)
    9.0

    This should be the same as the number of core nodes --- as boundary nodes
    in landlab do not have area.

    >>> len(mg.core_nodes)
    9

    Next, let's set the dx spacing such that each cell has an area of one.

    >>> dx = (2.0 / (3.0**0.5)) ** 0.5
    >>> hmg = HexModelGrid((5, 3), spacing=dx, xy_of_lower_left=(-1.0745, 0.0))
    >>> _ = hmg.add_field(
    ...     "topographic__elevation",
    ...     hmg.node_x**2 + np.round(hmg.node_y) ** 2,
    ...     at="node",
    ... )
    >>> fa = FlowAccumulator(
    ...     hmg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )
    >>> fa.run_one_step()
    >>> hmg.at_node["flow__receiver_node"]
    array([ 0,  1,  2,
            3,  0,  1,  6,
            7,  3,  4,  5, 11,
           12,  8,  9, 15,
           16, 17, 18])
    >>> np.round(hmg.at_node["drainage_area"])
    array([3., 2., 0.,
           2., 3., 2., 0.,
           0., 2., 2., 1., 0.,
           0., 1., 1., 0.,
           0., 0., 0.])

    Now let's change the cell area (100.) and the runoff rates:

    >>> hmg = HexModelGrid((5, 3), spacing=dx * 10.0, xy_of_lower_left=(-10.745, 0.0))

    Put the data back into the new grid.

    >>> _ = hmg.add_field(
    ...     "topographic__elevation",
    ...     hmg.node_x**2 + np.round(hmg.node_y) ** 2,
    ...     at="node",
    ... )
    >>> fa = FlowAccumulator(
    ...     hmg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )
    >>> fa.run_one_step()
    >>> np.round(hmg.at_node["surface_water__discharge"])
    array([500.,   0.,   0.,
           200., 500., 200.,   0.,
             0., 200., 200., 100.,   0.,
             0., 100., 100.,   0.,
             0.,   0.,   0.])

    Next, let's see what happens to a raster grid when there is a depression.

    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field("topographic__elevation", mg.node_x.copy(), at="node")
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)

    This model grid has a depression in the center.

    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)
    array([[0.    , 0.5   , 1.    , 1.5   , 2.    , 2.5   , 3.    ],
           [0.005 , 0.505 , 1.005 , 1.505 , 2.005 , 2.505 , 3.005 ],
           [0.01  , 0.51  , 0.101 , 0.151 , 0.201 , 2.51  , 3.01  ],
           [0.015 , 0.515 , 0.1015, 0.1515, 0.2015, 2.515 , 3.015 ],
           [0.02  , 0.52  , 0.102 , 0.152 , 0.202 , 2.52  , 3.02  ],
           [0.025 , 0.525 , 1.025 , 1.525 , 2.025 , 2.525 , 3.025 ],
           [0.03  , 0.53  , 1.03  , 1.53  , 2.03  , 2.53  , 3.03  ]])
    >>> fa = FlowAccumulator(
    ...     mg, "topographic__elevation", flow_director=FlowDirectorSteepest
    ... )
    >>> fa.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg.at_node["flow__receiver_node"].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14, 16, 16, 17, 18, 20],
           [21, 21, 16, 23, 24, 25, 27],
           [28, 28, 23, 30, 31, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 5.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 3.  , 0.75, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 2.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.5 , 0.25, 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ]])

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
    >>> mg.at_node["flow__receiver_node"].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 11, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 39, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
           [5.25, 5.25, 0.25, 0.25, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 5.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.75, 2.25, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 0.5 , 0.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.5 , 0.25, 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ]])

    Now the flow is routed correctly. The depression finder has properties that
    including whether there is a lake at the node, which lake is at each node,
    the outlet node of each lake, and the area of each lake.

    >>> df_4.lake_at_node.reshape(mg.shape)
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]])
    >>> df_4.lake_map.reshape(mg.shape)
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
    array([2.25])

    Alternatively, we can initialize a flow accumulator with a depression
    finder specified. Calling run_one_step() will run both the accumulator
    and the depression finder with one call. For this example, we will pass the
    class DepressionFinderAndRouter to the parameter `depression_finder`.

    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field("topographic__elevation", mg.node_x.copy(), at="node")
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fa = FlowAccumulator(
    ...     mg,
    ...     "topographic__elevation",
    ...     flow_director="FlowDirectorD8",
    ...     depression_finder=DepressionFinderAndRouter,
    ... )
    >>> fa.run_one_step()

    This has the same effect of first calling the accumulator and then calling
    the depression finder.

    >>> mg.at_node["flow__receiver_node"].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
           [5.25, 5.25, 0.25, 0.25, 0.25, 0.25, 0.  ],
           [0.25, 0.25, 5.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.75, 2.25, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 0.5 , 0.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ]])

    The depression finder is stored as part of the flow accumulator, so its
    properties can be accessed through the depression finder.

    >>> fa.depression_finder.lake_at_node.reshape(mg.shape)
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]])
    >>> fa.depression_finder.lake_map.reshape(mg.shape)
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
    array([2.25])

    Finally, note that the DepressionFinderAndRouter takes a keyword argument
    *routing* ('D8', default; 'D4') that sets how connectivity is set between
    nodes. Similar to our ability to pass keyword arguments to the FlowDirector
    through FlowAccumulator, we can pass this keyword argument to the
    DepressionFinderAndRouter component.

    >>> fa = FlowAccumulator(
    ...     mg,
    ...     "topographic__elevation",
    ...     flow_director=FlowDirectorSteepest,
    ...     depression_finder=DepressionFinderAndRouter,
    ...     routing="D4",
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
    >>> area = nmg.add_ones("cell_area_at_node", at="node")
    >>> z = nmg.add_field(
    ...     "topographic__elevation",
    ...     nmg.x_of_node + nmg.y_of_node,
    ...     at="node",
    ... )
    >>> fa = FlowAccumulator(nmg)
    >>> fa.run_one_step()
    >>> nmg.at_node["flow__receiver_node"]
    array([0, 0, 2, 1])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Braun, J., Willett, S. (2013). A very efficient O(n), implicit and parallel
    method to solve the stream power equation governing fluvial incision and
    landscape evolution. Geomorphology  180-181(C), 170-179.
    https://dx.doi.org/10.1016/j.geomorph.2012.10.008

    """

    _name = "FlowAccumulator"

    _unit_agnostic = True

    _info = {
        "drainage_area": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "flow__data_structure_delta": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Node array containing the elements delta[1:] of the data "
                "structure 'delta' used for construction of the downstream-to-upstream "
                "node array"
            ),
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water__unit_flux_in": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m/s",
            "mapping": "node",
            "doc": (
                "External volume water per area per time input to each node "
                "(e.g., rainfall rate)"
            ),
        },
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        flow_director="FlowDirectorSteepest",
        runoff_rate=None,
        depression_finder=None,
        **kwargs,
    ):
        """Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and
        compatability for the flow_director and depression_finder
        keyword arguments, tests the argument of runoff_rate, and
        initializes new fields.
        """
        super().__init__(grid)
        # Keep a local reference to the grid

        # Grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        self._is_Network = isinstance(self._grid, NetworkModelGrid)
        self._kwargs = kwargs

        # STEP 1: Testing of input values, supplied either in function call or
        # as part of the grid.
        self._test_water_inputs(grid, runoff_rate)

        # save elevations and node_cell_area to class properites.
        self._surface = surface
        self._surface_values = return_array_at_node(grid, surface)

        if self._is_Network:
            try:
                node_cell_area = self._grid.at_node["cell_area_at_node"]
            except FieldError as exc:
                raise FieldError(
                    "In order for the FlowAccumulator to work, the "
                    "grid must have an at-node field called "
                    "cell_area_at_node."
                ) from exc
        else:
            node_cell_area = self._grid.cell_area_at_node.copy()
            node_cell_area[self._grid.closed_boundary_nodes] = 0.0

        self._node_cell_area = node_cell_area

        # STEP 2:
        # This component will track the following variables.
        # Attempt to create each, if they already exist, assign the existing
        # version to the local copy.

        #   - drainage area at each node
        #   - receiver of each node
        #   - delta array

        self.initialize_output_fields()

        self._drainage_area = grid.at_node["drainage_area"]
        self._discharges = grid.at_node["surface_water__discharge"]

        self._upstream_ordered_nodes = grid.at_node["flow__upstream_node_order"]
        if np.all(self._upstream_ordered_nodes == 0):
            self._upstream_ordered_nodes.fill(self._grid.BAD_INDEX)

        self._delta_structure = grid.at_node["flow__data_structure_delta"]
        if np.all(self._delta_structure == 0):
            self._delta_structure[:] = self._grid.BAD_INDEX

        self._D_structure = self._grid.BAD_INDEX * grid.ones(at="link", dtype=int)
        self._nodes_not_in_stack = True

        # STEP 3:
        # identify Flow Director method, save name, import and initialize the
        # correct flow director component if necessary; same with
        # lake/depression handler, if specified.
        self._add_director(flow_director)
        self._add_depression_finder(depression_finder)

        if len(self._kwargs) > 0:
            kwdstr = " ".join(list(self._kwargs.keys()))
            raise ValueError(f"Extra kwargs passed to FlowAccumulator:{kwdstr}")

    @property
    def surface_values(self):
        """Values of the surface over which flow is accumulated."""
        return self._surface_values

    @property
    def flow_director(self):
        """The FlowDirector used internally."""
        return self._flow_director

    @property
    def depression_finder(self):
        """The DepressionFinder used internally."""
        return self._depression_finder

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
        ...     "topographic__elevation",
        ...     mg.node_x + mg.node_y,
        ...     at="node",
        ... )
        >>> fa = FlowAccumulator(mg, "topographic__elevation")
        >>> fa.run_one_step()
        >>> fa.link_order_upstream()
        array([ 5, 14, 23,  6, 15, 24,  7, 16, 25])

        This also works for route-to-many methods

        >>> mg = RasterModelGrid((5, 5))
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> np.flipud(
        ...     mg.add_field(
        ...         "topographic__elevation",
        ...         mg.node_x + mg.node_y,
        ...         at="node",
        ...     ).reshape(mg.shape)
        ... )
        array([[4., 5., 6., 7., 8.],
               [3., 4., 5., 6., 7.],
               [2., 3., 4., 5., 6.],
               [1., 2., 3., 4., 5.],
               [0., 1., 2., 3., 4.]])
        >>> fa = FlowAccumulator(mg, "topographic__elevation", flow_director="MFD")
        >>> fa.run_one_step()
        >>> link_order = fa.link_order_upstream()
        >>> link_order  # doctest: +SKIP
        array([ 5, 14, 10,  6, 11,  7, 23, 19, 15, 20, 16, 28, 24, 29, 25])
        >>> link_order[0]
        5
        >>> sorted(link_order[1:4])
        [6, 10, 14]
        >>> sorted(link_order[4:9])
        [7, 11, 15, 19, 23]
        >>> sorted(link_order[9:13])
        [16, 20, 24, 28]
        >>> sorted(link_order[13:])
        [25, 29]
        >>> np.all(sorted(link_order) == mg.active_links)
        True
        """
        downstream_links = self._grid["node"]["flow__link_to_receiver_node"][
            self._upstream_ordered_nodes
        ]
        out = downstream_links.flatten()
        return out[out != self._grid.BAD_INDEX]

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
        ...     "topographic__elevation",
        ...     mg.node_x + mg.node_y,
        ...     at="node",
        ... )
        >>> fa = FlowAccumulator(mg, "topographic__elevation")
        >>> fa.run_one_step()
        >>> assert_array_equal(fa.headwater_nodes(), np.array([16, 17, 18]))
        """
        delta = np.concatenate(([0], self._delta_structure))
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
                grid.add_ones("water__unit_flux_in", at="node", dtype=float)
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
            if p_k in self._kwargs:
                kw[p_k] = self._kwargs.pop(p_k)

        # flow director is provided as a string.
        if isinstance(flow_director, str):
            if flow_director[:12] == "FlowDirector":
                flow_director = flow_director[12:]

            from landlab.components.flow_director import FlowDirectorD8
            from landlab.components.flow_director import FlowDirectorDINF
            from landlab.components.flow_director import FlowDirectorMFD
            from landlab.components.flow_director import FlowDirectorSteepest

            DIRECTOR_METHODS = {
                "D4": FlowDirectorSteepest,
                "Steepest": FlowDirectorSteepest,
                "D8": FlowDirectorD8,
                "MFD": FlowDirectorMFD,
                "DINF": FlowDirectorDINF,
            }

            try:
                FlowDirector = DIRECTOR_METHODS[flow_director]
            except KeyError as exc:
                raise ValueError(
                    "String provided in flow_director is not a "
                    "valid method or component name. The following"
                    "components are valid imputs:\n" + str(PERMITTED_DIRECTORS)
                ) from exc
            self._flow_director = FlowDirector(self._grid, self._surface, **kw)
        # flow director is provided as an instantiated flow director
        elif isinstance(flow_director, Component):
            if flow_director._name in PERMITTED_DIRECTORS:
                self._flow_director = flow_director
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
                self._flow_director = FlowDirector(self._grid, self._surface, **kw)
            else:
                raise ValueError(
                    "String provided in flow_director is not a "
                    "valid method or component name. The following"
                    "components are valid imputs:\n" + str(PERMITTED_DIRECTORS)
                )

        # save method as attribute
        self._method = self._flow_director._method

    def _add_depression_finder(self, depression_finder):
        """Test and add the depression finder component."""
        PERMITTED_DEPRESSION_FINDERS = ["DepressionFinderAndRouter", "LakeMapperBarnes"]

        # now do a similar thing for the depression finder.
        self._depression_finder_provided = depression_finder
        if self._depression_finder_provided is not None:
            # collect potential kwargs to pass to depression_finder
            # instantiation
            potential_kwargs = [
                "routing",
                "pits",
                "reroute_flow",
                "surface",
                "method",
                "fill_flat",
                "fill_surface",
                "redirect_flow_steepest_descent",
                "reaccumulate_flow",
                "ignore_overfill",
                "track_lakes",
            ]
            kw = {}
            for p_k in potential_kwargs:
                if p_k in self._kwargs:
                    kw[p_k] = self._kwargs.pop(p_k)

            # NEED TO TEST WHICH FLOWDIRECTOR WAS PROVIDED.
            if self._flow_director._name in ("FlowDirectorMFD", "FlowDirectorDINF"):
                raise NotImplementedError(
                    "The depression finder only works with route "
                    "to one FlowDirectors such as "
                    "FlowDirectorSteepest and  FlowDirectorD8. "
                    "Provide a different FlowDirector."
                )

            # depression finder is provided as a string.
            if isinstance(self._depression_finder_provided, str):
                from landlab.components.depression_finder.lake_mapper import (
                    DepressionFinderAndRouter,
                )
                from landlab.components.lake_fill.lake_fill_barnes import (
                    LakeMapperBarnes,
                )

                DEPRESSION_METHODS = {
                    "DepressionFinderAndRouter": DepressionFinderAndRouter,
                    "LakeMapperBarnes": LakeMapperBarnes,
                }

                try:
                    DepressionFinder = DEPRESSION_METHODS[
                        self._depression_finder_provided
                    ]
                except KeyError as exc:
                    raise ValueError(
                        "Component provided in depression_finder "
                        "is not a valid component. The following "
                        "components are valid imputs: "
                        f"{', '.join(repr(x) for x in PERMITTED_DEPRESSION_FINDERS)}."
                    ) from exc

                self._depression_finder = DepressionFinder(self._grid, **kw)
            # flow director is provided as an instantiated depression finder
            elif isinstance(self._depression_finder_provided, Component):
                if (
                    self._depression_finder_provided._name
                    in PERMITTED_DEPRESSION_FINDERS
                ):
                    self._depression_finder = self._depression_finder_provided
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

            # depression_finder is provided as an uninstantiated depression finder
            else:
                if (
                    self._depression_finder_provided._name
                    in PERMITTED_DEPRESSION_FINDERS
                ):
                    DepressionFinder = self._depression_finder_provided
                    self._depression_finder = DepressionFinder(self._grid, **kw)
                else:
                    raise ValueError(
                        "Component provided in depression_finder "
                        "is not a valid component. The following "
                        "components are valid imputs:\n"
                        + str(PERMITTED_DEPRESSION_FINDERS)
                    )

            # Make sure direction methods are consistent between the director
            # and the depression handler
            if isinstance(self._grid, RasterModelGrid):
                flow_director_method = self.flow_director_raster_method()
                depression_finder_method = (
                    self.depression_handler_raster_direction_method()
                )
                if flow_director_method != depression_finder_method:
                    message = (
                        "Incompatibility between flow-director routing method\n"
                        + "which is "
                        + flow_director_method
                        + ", and depression-handler method,\n"
                        + "which is "
                        + depression_finder_method
                    )
                    raise ValueError(warning_message(message))
        else:
            self._depression_finder = None

    def flow_director_raster_method(self):
        """Return 'D8' or 'D4' depending on the direction method used.

        (Note: only call this function for a raster gird;
        does not handle multiple-flow directors)
        """
        assert isinstance(self._grid, RasterModelGrid)
        if self._flow_director._name in ("FlowDirectorD8"):
            return "D8"
        else:
            return "D4"

    def depression_handler_raster_direction_method(self):
        """Return 'D8' or 'D4' depending on the direction method used.

        (Note: only call this function for a raster gird;
        does not handle multiple-flow directors)
        """
        assert isinstance(self._grid, RasterModelGrid)
        if self._depression_finder._name in ("DepressionFinderAndRouter"):
            return self._depression_finder._routing
        elif self._depression_finder._name in ("LakeMapperBarnes"):
            if (
                self._depression_finder._allneighbors.size
                > self.grid.adjacent_nodes_at_node.size
            ):
                return "D8"
            else:
                return "D4"
        else:
            raise ValueError("Depression finder type not recognized.")

    def pits_present(self):
        return np.any(self._grid.at_node["flow__sink_flag"][self._grid.core_nodes])

    def flooded_nodes_present(self):
        # flooded node status may not exist if no depression finder was used.
        if "flood_status_code" in self._grid.at_node:
            return np.all(
                self._grid.at_node["flood_status_code"] == FloodStatus.UNFLOODED
            )
        else:
            return False

    def accumulate_flow(self, update_flow_director=True, update_depression_finder=True):
        """Function to make FlowAccumulator calculate drainage area and
        discharge.

        Running accumulate_flow() results in the following to occur:

            1. Flow directions are updated (unless update_flow_director is set
               as False). This incluldes checking for updated boundary
               conditions.
            2. The depression finder, if present is updated (unless
               update_depression_finder is set as False).
            3. Intermediate steps that analyse the drainage network topology
               and create datastructures for efficient drainage area and
               discharge calculations.
            4. Calculation of drainage area and discharge.
            5. Return of drainage area and discharge.

        Parameters
        ----------
        update_flow_director : optional, bool
            Whether to update the flow director. Default is True.
        update_depression_finder : optional, bool
            Whether to update the depression finder, if present.
            Default is True.

        Returns
        -------
        drainage_area : array
            At node array which points to the field
            grid.at_node["drainage_area"].
        surface_water__discharge
            At node array which points to the field
            grid.at_node["surface_water__discharge"].
        """
        # set a couple of aliases
        a = self._grid["node"]["drainage_area"]
        q = self._grid["node"]["surface_water__discharge"]

        # step 1. Find flow directions by specified method
        if update_flow_director:
            self._flow_director.run_one_step()

        # further steps vary depending on how many recievers are present
        # one set of steps is for route to one (D8, Steepest/D4)

        # step 2. Get r
        r = as_id_array(self._grid["node"]["flow__receiver_node"])

        if self._flow_director._to_n_receivers == "one":
            # step 2b. Run depression finder if passed
            # Depression finder reaccumulates flow at the end of its routine.
            # At the moment, no depression finders work with to-many, so it
            # lives here
            if (
                self._depression_finder_provided is not None
                and update_depression_finder
            ):
                # only update depression finder if requested AND if there
                # are pits, or there were flooded nodes from last timestep.
                if self.pits_present or self.flooded_nodes_present:
                    self._depression_finder.update()

                # if FlowDirectorSteepest is used, update the link directions
                if self._flow_director._name == "FlowDirectorSteepest":
                    self._flow_director._determine_link_directions()

            # step 3. Stack, D, delta construction
            nd = as_id_array(flow_accum_bw._make_number_of_donors_array(r))
            delta = as_id_array(flow_accum_bw._make_delta_array(nd))
            D = as_id_array(flow_accum_bw._make_array_of_donors(r, delta))
            s = as_id_array(flow_accum_bw.make_ordered_node_array(r, nd, delta, D))

            # put these in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid.at_node["flow__data_structure_delta"][:] = as_id_array(delta[1:])
            self._D_structure = as_id_array(D)
            self._grid.at_node["flow__upstream_node_order"][:] = as_id_array(s)

            # step 4. Accumulate (to one or to N depending on direction method)
            a[:], q[:] = self._accumulate_A_Q_to_one(s, r)

        else:
            # Get p
            p = self._grid["node"]["flow__receiver_proportions"]

            # step 3. Stack, D, delta construction
            nd = as_id_array(flow_accum_to_n._make_number_of_donors_array_to_n(r, p))
            delta = as_id_array(flow_accum_to_n._make_delta_array_to_n(nd))
            D = as_id_array(flow_accum_to_n._make_array_of_donors_to_n(r, p, delta))
            s = as_id_array(
                flow_accum_to_n.make_ordered_node_array_to_n(r, p, nd, delta, D)
            )

            # put theese in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid["node"]["flow__data_structure_delta"][:] = delta[1:]
            self._D_structure = D

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
            s, r, self._node_cell_area, self._grid.at_node["water__unit_flux_in"]
        )
        return (a, q)

    def _accumulate_A_Q_to_n(self, s, r, p):
        """Accumulate area and discharge for a route-to-many scheme.

        Note this can be overridden in inherited components.
        """
        a, q = flow_accum_to_n.find_drainage_area_and_discharge_to_n(
            s, r, p, self._node_cell_area, self._grid.at_node["water__unit_flux_in"]
        )
        return (a, q)

    def run_one_step(self):
        """Accumulate flow and save to the model grid.

        1. Flow directions are updated. This incluldes checking for updated
           boundary conditions.
        2. The depression finder, if present is updated.
        3. Intermediate steps that analyse the drainage network topology
           and create datastructures for efficient drainage area and
           discharge calculations.
        4. Calculation of drainage area and discharge.
        5. Return of drainage area and discharge.

        An alternative to run_one_step() is accumulate_flow() which does the
        same things but also returns the drainage area and discharge.
        accumulate_flow() additionally provides the ability to turn off updating
        the flow director or the depression finder.
        """
        self.accumulate_flow()


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
