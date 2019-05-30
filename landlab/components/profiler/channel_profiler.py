# coding: utf8
# ! /usr/env/python
"""channel_profiler.py component to create channel profiles."""
from collections import OrderedDict

import numpy as np

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt

from landlab import RasterModelGrid
from landlab.components.profiler.base_profiler import _BaseProfiler
from landlab.utils.flow__distance import calculate_flow__distance


class ChannelProfiler(_BaseProfiler):
    """Extract and plot the channel profiles in drainage networks.

    The ChannelProfiler extracts channel networks from a landlab grid.

    In order to extract channel networks, the flow connectivity across the grid
    must already be identified. This is typically done with the FlowAccumulator
    component. However, this component does not require that the
    FlowAccumulator was used. Instead it expects that the following at-node
    grid fields will be present:

        'flow__receiver_node'
        'flow__link_to_receiver_node'

    The ChannelProfiler can work on grids that have used route-to-one or
    route-to-multiple flow directing.

    To understand how this component works it is useful to define the following
    terms: *watershed*, *outlet*, *headwater node*, *segment*.

    A *watershed* is all model grid nodes that drain to a single node, called
    the *outlet*. Channels nodes are identified as nodes that have a
    ``channel_definition_field`` value greater than or equal to the
    ``minimum_channel_threshold``. This ``channel_definition_field`` is often
    the drainage area (this component's default). We use a flexible field
    rather than only drainage area to support alternative bases for channel
    extraction.

    The default behaviour of this component is to use an exclusive definition
    of a *watershed*. That is, the two largest watersheds are defined as the
    watersheds upstream of the two nodes on the model grid boundary with the
    largest values of the ``channel_definition_field`` rather than potentially
    nested watersheds. Nested watersheds are supported through the use of the
    ``outlet_nodes`` keyword argument.

    Consider the following grid with 10 columns and 7 rows. In this grid is one
    watershed with an outlet node indicated by ``o``. Here ``X`` indicates
    nodes that are not part of the channel network (based on the
    ``channel_definition_field``) and ``.`` indicates nodes that are part of
    the network.

    In this and the following examples, we will use only D4 connectivity. The
    ChannelProfiler, however, knows nothing of connectivity other than what is
    implied by the two required grid fields.
    ::

        X X X X X X X X X X
        X . X X X X X X X X
        X . . X X X . . . X
        X X . . X X . X X X
        X X X . . . . X X X
        X X X . X X X X X X
        X X X o X X X X X X

    This component can extract the channel network from one or more watersheds.
    This option is specified with the keyword argument
    ``number_of_watersheds``.

    The *headwater nodes*, shown as ``@`` are nodes that have no upstream
    nodes with sufficient area to be classified as a channel.
    ::

        X X X X X X X X X X
        X @ X X X X X X X X
        X . . X X X . . @ X
        X X . . X X . X X X
        X X X . . . . X X X
        X X X . X X X X X X
        X X X o X X X X X X

    For each watershed, the ChannelProfiler either extracts the largest channel
    (again, based on the ``channel_definition_field``) or all channel segments
    with sufficent values in the ``channel_definition_field``.

    Default behavior of this component is to extract only the largest channel
    in the single largest watershed. This would extract the following channel
    segment.
    ::

        X X X X X X X X X X
        X . X X X X X X X X
        X . . X X X * * * X
        X X . . X X * X X X
        X X X * * * * X X X
        X X X * X X X X X X
        X X X * X X X X X X

    This component verifies that all watershed outlets have a value in the
    ``channel_definition_field`` of at least ``minimum_outlet_threshold``
    (default is 0 units). If no watersheds exist that meet this criteria, an
    error is raised.

    If a user knows exactly which node or nodes they want to use as the outlet
    nodes, then this can be specified using the ``outlet_nodes`` keyword
    argument. Otherwise the ``number_of_watersheds`` (default 1) nodes with the
    largest value in the ``channel_definition_field`` will be selected as the
    outlet nodes from the model grid boundary nodes. Setting
    ``number_of_watersheds`` to ``None`` results in selecting all nodes at the
    model grid boundary that meet the criteria for an outlet based on the
    ``channel_definition_field`` and the ``minimum_outlet_threshold``.

    The node IDs and distances upstream of the channel network are stored in a
    data structure called ``network_structure``. It is a dictionary with keys
    indicating the outlet node.

    For each watershed outlet, the value in the ``network_structure`` is itself
    a dictionary with keys that are a segment ID tuple of the
    ``(dowstream, upstream)`` nodes IDs of each channel segment.

    For our simple example, these are the node IDs:
    ::

            X  X  X  X  X  X  X  X  X  X
            X 51  X  X  X  X  X  X  X  X
            X 41 42  X  X  X 46 47 48  X
            X  X 32 33  X  X 36  X  X  X
            X  X  X 23 24 25 26  X  X  X
            X  X  X 13  X  X  X  X  X  X
            X  X  X  3  X  X  X  X  X  X

    So for our main channel only example, the outlet has an ID of 3, the
    downstream end of the channel segment is 3, and the upstream end is 48.

    The value associated with the segment ID tuple ``(3, 48)`` is itself a
    dictionary. It has three key-value pairs. First, ``"ids"`` contains a list
    of the segment node ids ordered from downstream to upstream. It includes
    the endpoints. Second, ``"distances"`` contains a list of distances
    upstream that mirrors the list in ``"ids"``. Finally, ``"color"`` is an
    RGBA tuple indicating the color for the segment.

    By default a unique color will be assigned to each watershed. To change the
    color, a user can change values stored in ``network_structure``.
    Additionally, a ``cmap`` keyword argument can provide some user control
    over the color at the instantiation of the component.

    In the main channel only example, the datastructure will look as follows:

    . code-block:: python

        {3: {
            (3, 48) : {
                "ids": [3, 13, 23, 24, 25, 26, 36, 46, 47, 48],
                "distances": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
                "color": (1, 0, 1, 1),
                }
            }
        }

    Three channel segments are idendified if ``main_channel_only=False``.
    ::

        X X X X X X X X X X     X X X X X X X X X X     X X X X X X X X X X
        X . X X X X X X X X     X . X X X X X X X X     X * X X X X X X X X
        X . . X X X . . . X     X . . X X X * * * X     X * * X X X . . . X
        X X . . X X . X X X     X X . . X X * X X X     X X * * X X . X X X
        X X X * . . . X X X     X X X * * * * X X X     X X X * . . . X X X
        X X X * X X X X X X     X X X . X X X X X X     X X X . X X X X X X
        X X X * X X X X X X     X X X . X X X X X X     X X X . X X X X X X

    The data structure associated with this set of three segments is

        . code-block:: python

            {3: {
                (3, 23) : {
                    "ids": [3, 13, 23],
                    "distances": [0, 1, 2]
                    "color": (1, 0, 1, 1),
                    },
                (23, 48) : {
                    "ids": [23, 24, 25, 26, 36, 46, 47, 48],
                    "distances": [2, 3, 4, 5, 6, 7, 8, 9],
                    "color": (1, 0, 1, 1),
                    },
                (23, 51) : {
                    "ids": [23, 33, 32, 42, 41, 51],
                    "distances": [2, 3, 4, 5, 6, 7]
                    "color": (1, 0, 1, 1),
                    },
                }
            }

    Note that the distances upstream are relative to the outlet, not the
    downstream end of the stream segment.

    Next consider a model grid with two watersheds.
    ::

        X X X X X X X X X X
        X . . X X X . X X X
        X . X X . X . X X X
        o . . . . X . X X X
        X X X X X X . X X X
        X X X . . . . X X X
        X X X X X X . . . X
        X X X X X X X X o X

    And the following node IDs.
    ::

        X   X   X   X   X   X   X   X   X   X
        X  61  62   X   X   X  66   X   X   X
        X  51   X   X  54   X  56   X   X   X
       40  41  42  43  44   X  46   X   X   X
        X   X   X   X   X   X  36   X   X   X
        X   X   X  23  24  25  26   X   X   X
        X   X   X   X   X   X  16  17  18   X
        X   X   X   X   X   X   X   X   8   X


    The data structure for ``number_of_watersheds=2`` and
    ``main_channel_only=False`` will be as follows. Note that each watershed
    has been assigned a different color tuple value.

    . code-block:: python

        {8: {
            (8, 26) : {
                "ids": [8, 18, 17, 16, 26],
                "distances": [0, 1, 2, 3, 4]
                "color": (1, 0, 1, 1),
                },
            (26, 23) : {
                "ids": [26, 25, 24, 23],
                "distances": [4, 5, 6, 7],
                "color": (1, 0, 1, 1),
                },
            (26, 66) : {
                "ids": [26, 36, 46, 56, 66],
                "distances": [4, 5, 6, 7, 8]
                "color": (1, 0, 1, 1),
                },
            },
         40: {
            (40, 41) : {
                "ids": [40, 41],
                "distances": [0, 1]
                "color": (1, 0, 0, 1),
                },
            (41, 54) : {
                "ids": [41, 42, 43, 44, 54],
                "distances": [2, 3, 4, 5, 6],
                "color": (1, 0, 0, 1),
                },
            (41, 62) : {
                "ids": [41, 51, 61, 62],
                "distances": [1, 2, 3, 4]
                "color": (1, 0, 0, 1),
                },
            }
        }

    Examples
    ---------

    Start by importing necessary modules

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, ChannelProfiler

    Create the second example grid we showed above. Note that in order to do
    this we need to enter the elevations starting from the lower left so the
    elevation order may seem upside-down    .

    >>> z = np.array([ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,
    ...                0,  0,  0,  0,  0,  0,  4,  3,  2,  0,
    ...                0,  0,  0,  8,  7,  6,  5,  0,  0,  0,
    ...                0,  0,  0,  0,  0,  0,  6,  0,  0,  0,
    ...                1,  3,  4,  5,  6,  0,  7,  0,  0,  0,
    ...                0,  4,  0,  0,  7,  0,  8,  0,  0,  0,
    ...                0,  5,  6,  0,  0,  0,  9,  0,  0,  0,
    ...                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,], dtype=float)

    >>> mg = RasterModelGrid((8, 10))
    >>> z = mg.add_field("node", "topographic__elevation", z)
    >>> mg.set_nodata_nodes_to_closed(z, 0)
    >>> fa = FlowAccumulator(mg, flow_director='D4')
    >>> fa.run_one_step()
    >>> fa.drainage_area.reshape(mg.shape)
    array([[  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  11.,   0.],
           [  0.,   0.,   0.,   0.,   0.,   0.,   9.,  10.,  11.,   0.],
           [  0.,   0.,   0.,   1.,   2.,   3.,   8.,   0.,   0.,   0.],
           [  0.,   0.,   0.,   0.,   0.,   0.,   4.,   0.,   0.,   0.],
           [  8.,   8.,   4.,   3.,   2.,   0.,   3.,   0.,   0.,   0.],
           [  0.,   3.,   0.,   0.,   1.,   0.,   2.,   0.,   0.,   0.],
           [  0.,   2.,   1.,   0.,   0.,   0.,   1.,   0.,   0.,   0.],
           [  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.]])

    >>> profiler = ChannelProfiler(
    ...     mg,
    ...     number_of_watersheds=2,
    ...     minimum_channel_threshold=0,
    ...     main_channel_only=False)
    >>> profiler.run_one_step()

    Todo

    >>> profiler.network_structure.keys()
    odict_keys([40, 8])

    Todo

    >>> profiler.network_structure[40].keys()
    odict_keys([(40, 41), (41, 54), (41, 62)])

    Todo

    >>> profiler.network_structure[40][(40, 41)]["ids"]
    array([40, 41])
    >>> profiler.network_structure[40][(40, 41)]["distances"]
    array([ 0.,  1.])
    >>> np.round(profiler.network_structure[40][(40, 41)]["color"], decimals=2)
    array([ 0.27,  0.  ,  0.33,  1.  ])

    >>> profiler.network_structure[40][(41, 54)]["ids"]
    array([41, 42, 43, 44, 54])
    >>> profiler.network_structure[40][(41, 54)]["distances"]
    array([ 1.,  2.,  3.,  4.,  5.])
    >>> np.round(profiler.network_structure[40][(41, 54)]["color"], decimals=2)
    array([ 0.27,  0.  ,  0.33,  1.  ])

    >>> profiler.network_structure[40][(41, 62)]["ids"]
    array([41, 51, 61, 62])
    >>> profiler.network_structure[40][(41, 62)]["distances"]
    array([ 1.,  2.,  3.,  4.])
    >>> np.round(profiler.network_structure[40][(41, 62)]["color"], decimals=2)
    array([ 0.27,  0.  ,  0.33,  1.  ])

    Todo

    >>> profiler.network_structure[8].keys()
    odict_keys([(8, 26), (26, 23), (26, 66)])

    Todo

    >>> profiler.network_structure[8][(8, 26)]["ids"]
    array([ 8, 18, 17, 16, 26])
    >>> profiler.network_structure[8][(8, 26)]["distances"]
    array([ 0.,  1.,  2.,  3.,  4.])
    >>> np.round(profiler.network_structure[8][(8, 26)]["color"], decimals=2)
    array([ 0.13,  0.57,  0.55,  1.  ])

    >>> profiler.network_structure[8][(26, 23)]["ids"]
    array([26, 25, 24, 23])
    >>> profiler.network_structure[8][(26, 23)]["distances"]
    array([ 4.,  5.,  6.,  7.])
    >>> np.round(profiler.network_structure[8][(26, 23)]["color"], decimals=2)
    array([ 0.13,  0.57,  0.55,  1.  ])

    >>> profiler.network_structure[8][(26, 66)]["ids"]
    array([26, 36, 46, 56, 66])
    >>> profiler.network_structure[8][(26, 66)]["distances"]
    array([ 4.,  5.,  6.,  7.,  8.])
    >>> np.round(profiler.network_structure[8][(26, 66)]["color"], decimals=2)
    array([ 0.13,  0.57,  0.55,  1.  ])

    The ChannelProfiler is designed to be flexible, and by careful combination
    of its instantiation variables can be used to extract many useful forms of
    profile. In these examples, we will use the default
    ``channel_definition_field``, the drainage area.

    To illustrate, lets start by creating a landscape model.

    >>> from landlab.components import FastscapeEroder
    >>> mg = RasterModelGrid((40, 60), xy_spacing=100)
    >>> z = mg.add_zeros('topographic__elevation', at='node')
    >>> fa = FlowAccumulator(mg)
    >>> sp = FastscapeEroder(mg, K_sp=.0001)
    >>> dt = 100
    >>> for i in range(200):
    ...     fa.run_one_step()
    ...     sp.run_one_step(dt=dt)
    ...     z[mg.core_nodes] += 0.001 * dt

    Some options:

    Default: Extract a the single biggest channel draining to the model grid
    boundary traced back all the way to the watershed divide.

    >>> profiler = ChannelProfiler(mg)

    Extract the largest channel draining to each of the four largest outlet
    nodes on the model grid boundary traced back all the way to the watershed
    divide.

    >>> profiler = ChannelProfiler(mg, number_of_watersheds=4)

    Extract the single largest channel draining to node 84. Note that the
    keyword argument ``outlet_nodes`` must be an iterable.

    >>> profiler = ChannelProfiler(mg, outlet_nodes=[84])

    Extract the largest channel draining to each of the four largest outlet
    nodes on the model grid boundary traced back to nodes with
    ``channel_definition_field`` values of 500.

    >>> profiler = ChannelProfiler(
    ...     mg,
    ...     number_of_watersheds=4,
    ...     minimum_channel_threshold=500)

    Extract a the single biggest channel draining to the model grid boundary
    based on the field ``surface_water__discharge`` traced back to discharge
    values of 500.

    >>> profiler = ChannelProfiler(mg,
    ...     channel_definition_field='surface_water__discharge',
    ...     minimum_channel_threshold=500)

    Extract the single largest channel within *all* watersheds with an outlet
    with ``channel_definition_field`` greater than 1e3. Trace the channels
    up to the point in each watershed in which the channels have values in the
    ``channel_definition_field`` of 500.

    >>> profiler = ChannelProfiler(
    ...     mg,
    ...     number_of_watersheds=None,
    ...     minimum_outlet_threshold=1e3,
    ...     minimum_channel_threshold=500)

    Extract two trunk channels beginning at the given nodes, traced up to a
    a minimum ``channel_definition_field`` value of of 500. Note that
    ``number_of_watersheds`` must match the size of ``outlet_nodes``.

    >>> profiler = ChannelProfiler(
    ...     mg,
    ...     outlet_nodes=[24, 45],
    ...     number_of_watersheds=2,
    ...     minimum_channel_threshold=500)

    Extract every possible channel (not just the largest one), leading from the
    two highest model grid boundary nodes traced back to a
    ``channel_definition_field`` threshold of 500.
    >>> profiler = ChannelProfiler(mg,
    ...     main_channel_only=False,
    ...     minimum_channel_threshold=500)
    """

    _name = "ChannelProfiler"

    _input_var_names = (
        "topographic__elevation",
        "drainage_area",
        "flow__receiver_node",
        "flow__link_to_receiver_node",
    )

    _output_var_names = ()

    _var_units = {
        "topographic__elevation": "m",
        "flow__receiver_node": "-",
        "drainage_area": "m**2",
        "flow__link_to_receiver_node": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "drainage_area": "node",
        "flow__link_to_receiver_node": "node",
    }
    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current node)",
        "drainage_area": "Upstream accumulated surface area contributing to the node's discharge",
        "flow__link_to_receiver_node": "Node array containing ID of link that leads from each node to its receiver, or BAD_INDEX_VALUE if no link",
    }

    def __init__(
        self,
        grid,
        channel_definition_field="drainage_area",
        number_of_watersheds=1,
        minimum_outlet_threshold=0,
        main_channel_only=True,
        outlet_nodes=None,
        minimum_channel_threshold=0,
        cmap="viridis",
    ):
        """
        Parameters
        ----------
        grid : Landlab Model Grid instance, required
        channel_definition_field : field name as string
            Name of field used to identify the outlet and headwater nodes of the
            channel network. Default is "drainage_area".
        minimum_outlet_threshold : float, optional
            Minimum value of the *channel_definition_field* to define a
            watershed outlet. Default is 0.0.
        minimum_channel_threshold : float, optional
            Value to use for the minimum drainage area associated with a
            plotted channel segment. Default values 0.
        number_of_watersheds : int, optional
            Total number of watersheds to plot. Default value is 1. If value is
            greater than 1 and outlet_nodes is not specified, then the
            number_of_watersheds largest watersheds is based on the drainage
            area at the model grid boundary. If given as None, then all grid
            cells on the domain boundary with a stopping field (typically
            drainage area) greater than the minimum_outlet_threshold in area are used.
        main_channel_only : Boolean, optional
            Flag to determine if only the main channel should be plotted, or if
            all stream segments with drainage area less than threshold should
            be plotted. Default value is True.
        outlet_nodes : length number_of_watersheds iterable, optional
            Length number_of_watersheds iterable containing the node IDs of
            nodes to start the channel profiles from. If not provided, the
            default is the number_of_watersheds node IDs on the model grid
            boundary with the largest terminal drainage area.
        cmap : str
            A valid matplotlib cmap string. Default is "viridis".

        """
        super(ChannelProfiler, self).__init__(grid)

        self._grid = grid
        self._cmap = plt.get_cmap(cmap)
        if channel_definition_field in grid.at_node:
            self._channel_definition_field = grid.at_node[channel_definition_field]
        else:
            msg = "Required field {name} not present. This field is required by the ChannelProfiler to define the start and stop of channel networks."
            raise ValueError(msg)

        if "flow__receiver_node" in grid.at_node:
            self._flow_receiver = grid.at_node["flow__receiver_node"]
        else:
            msg = "flow__receiver_node is a required field to run a ChannelProfiler."
            raise ValueError(msg)

        if "flow__link_to_receiver_node" in grid.at_node:
            self._link_to_flow_receiver = grid.at_node["flow__link_to_receiver_node"]
        else:
            msg = "flow__link_to_receiver_node is a required field to run a ChannelProfiler."
            raise ValueError(msg)

        self._main_channel_only = main_channel_only
        self.minimum_channel_threshold = minimum_channel_threshold

        # verify that the number of starting nodes is the specified number of channels
        if outlet_nodes is not None:
            if (number_of_watersheds is not None) and (
                len(outlet_nodes) is not number_of_watersheds
            ):
                msg = "Length of outlet_nodes must equal the" "number_of_watersheds!"
                raise ValueError(msg)
        else:
            large_outlet_ids = grid.boundary_nodes[
                np.argsort(self._channel_definition_field[grid.boundary_nodes])
            ]
            if number_of_watersheds is None:
                big_enough_watersheds = (
                    self._channel_definition_field[large_outlet_ids]
                    > minimum_outlet_threshold
                )
                outlet_nodes = large_outlet_ids[big_enough_watersheds]
            else:
                outlet_nodes = large_outlet_ids[-number_of_watersheds:]

        starting_da = self._channel_definition_field[outlet_nodes]
        outlet_nodes = np.asarray(outlet_nodes)
        if (
            np.any(starting_da < self.minimum_channel_threshold)
            or outlet_nodes.size == 0
        ):
            msg = (
                "The number of watersheds requested by the ChannelProfiler is "
                "greater than the number in the domain with channel_definition_field"
                " area. {vals}".format(vals=starting_da)
            )
            raise ValueError(msg)

        self._outlet_nodes = outlet_nodes

    @property
    def network_structure(self):
        """
        OrderedDict defining the channel network.

        The node IDs and distances upstream of the channel network are stored
        in ``network_structure``. It is a dictionary with keys of the outlet
        node ID.

        For each watershed outlet, the value in the ``network_structure`` is
        itself a dictionary with keys that are a segment ID tuple of the
        ``(dowstream, upstream)`` nodes IDs of each channel segment.

        The value associated with the segment ID tuple
        ``(dowstream, upstream)`` is itself a dictionary. It has three
        key-value pairs. First, ``"ids"`` contains a list of the segment node
        IDs ordered from downstream to upstream. It includes the endpoints.
        Second, ``"distances"`` contains a list of distances upstream that
        mirrors the list in ``"ids"``. Finally, ``"color"`` is an RGBA tuple
        indicating the color for the segment.
        """
        return self._net_struct

    @property
    def distance_along_profile(self):
        """ """
        return self._distance_along_profile

    @property
    def network_ids(self):
        """ """
        return self._net_ids

    def _get_channel_segment(self, i):
        """Get channel segment and return additional nodes to process.

        Parameters
        ----------
        i : int, required
            Node id of start of channel segment.

        Returns
        ----------
        channel_segment : list
            Node IDs of the nodes in the current channel segment.
        nodes_to_process, list
            List of nodes to add to the processing queue. These nodes are those
            that drain to the upper end of this channel segment. If
            main_channel_only = False this will be an empty list.
        """
        j = i
        channel_segment = []
        channel_upstream = True

        # add the reciever of j to the channel segment if it is not j.
        recieving_node = self._flow_receiver[j]
        if recieving_node != j:
            channel_segment.append(recieving_node)

        while channel_upstream:

            # add the new node to the channel segment
            channel_segment.append(j)

            # get supplying nodes
            supplying_nodes = np.where(self._flow_receiver == j)[0]

            # remove supplying nodes that are the outlet node
            supplying_nodes = supplying_nodes[np.where(supplying_nodes != i)]

            # if only adding the biggest channel, continue upstream choosing the
            # largest node until no more nodes remain.
            if (self._main_channel_only) and (len(supplying_nodes) > 0):
                max_drainage = np.argmax(
                    self._channel_definition_field[supplying_nodes]
                )

                if (
                    self._channel_definition_field[supplying_nodes[max_drainage]]
                    < self.minimum_channel_threshold
                ):
                    nodes_to_process = []
                    channel_upstream = False
                else:
                    j = supplying_nodes[max_drainage]

            # if considering multiple channel segments, continue upstream until
            # there are two or more donors with sufficient discharge, then
            # break, returning those nodes as starting points.
            else:

                # get all upstream drainage areas
                upstream_das = self._channel_definition_field[supplying_nodes]

                # if no nodes upstream exceed the threshold, exit
                if np.sum(upstream_das > self.minimum_channel_threshold) == 0:
                    nodes_to_process = []
                    channel_upstream = False

                # otherwise
                else:
                    # if only one upstream node exceeds the threshold, proceed
                    # up the channel.
                    if np.sum(upstream_das > self.minimum_channel_threshold) == 1:
                        max_drainage = np.argmax(
                            self._channel_definition_field[supplying_nodes]
                        )
                        j = supplying_nodes[max_drainage]
                    # otherwise provide the multiple upstream nodes to be
                    # processed into a new channel.
                    else:
                        nodes_to_process = supplying_nodes[
                            upstream_das > self.minimum_channel_threshold
                        ]
                        channel_upstream = False

        return (channel_segment, nodes_to_process)

    def _create_profile_structure(self):
        """
        Create the profile_IDs data structure for channel network.

        The bound attribute self._profile structure is the channel segment
        datastructure. profile structure is a list of length
        number_of_watersheds. Each element of profile_structure is itself a
        list of length number of stream segments that drain to each of the
        starting nodes. Each stream segment list contains the node ids of a
        stream segment from downstream to upstream.
        """
        self._net_struct = OrderedDict()

        if self._main_channel_only:
            for i in self._outlet_nodes:
                (channel_segment, nodes_to_process) = self._get_channel_segment(i)
                segment_tuple = (channel_segment[0], channel_segment[-1])
                self._net_struct[i] = {
                    segment_tuple: {"ids": np.array(channel_segment)}
                }

        else:
            for i in self._outlet_nodes:
                channel_network = OrderedDict()
                queue = [i]
                while len(queue) > 0:
                    node_to_process = queue.pop(0)
                    (channel_segment, nodes_to_process) = self._get_channel_segment(
                        node_to_process
                    )
                    segment_tuple = (channel_segment[0], channel_segment[-1])
                    channel_network[segment_tuple] = {"ids": np.array(channel_segment)}
                    queue.extend(nodes_to_process)
                self._net_struct[i] = channel_network

        self._calculate_distances()
        self._assign_colors()
        self._create_flat_structures()

    def _create_flat_structures(self):
            """ """
            self._net_ids = []
            self._distance_along_profile = []
            self._colors = []
            for outlet_id in self._net_struct:
                seg_tuples = self._net_struct[outlet_id].keys()
                self._net_ids.append([self._net_struct[outlet_id][seg]["ids"] for seg in seg_tuples])
                self._distance_along_profile.append([self._net_struct[outlet_id][seg]["distances"] for seg in seg_tuples])
                self._colors.append([self._net_struct[outlet_id][seg]["color"] for seg in seg_tuples])

    def _assign_colors(self):
        """Assign colors"""
        num_watersheds = len(self._net_struct)
        norm = mpl.colors.Normalize(vmin=0, vmax=num_watersheds)
        mappable = cm.ScalarMappable(norm=norm, cmap=self._cmap)
        colors = [mappable.to_rgba(val) for val in range(num_watersheds)]

        for idx, outlet_id in enumerate(self._net_struct):
            for segment_tuple in self._net_struct[outlet_id]:
                self._net_struct[outlet_id][segment_tuple]["color"] = colors[idx]

    def _calculate_distances(self):
        """Get distances along the network datastructure"""
        distance_upstream = calculate_flow__distance(self._grid)
        for outlet_id in self._net_struct:
            offset = distance_upstream[outlet_id]
            for segment_tuple in self._net_struct[outlet_id]:
                ids = self._net_struct[outlet_id][segment_tuple]["ids"]
                d = distance_upstream[ids]
                self._net_struct[outlet_id][segment_tuple]["distances"] = d - offset
