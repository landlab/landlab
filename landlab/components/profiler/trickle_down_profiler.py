# ! /usr/env/python
"""trickle_down_profiler.py component to create channel profiles."""
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

from landlab.components.profiler.base_profiler import _BaseProfiler
from landlab.core.utils import as_id_array
from landlab.utils.flow__distance import calculate_flow__distance


class TrickleDownProfiler(_BaseProfiler):
    """Extract and a profile from one or more node IDs to their downstream termini.

    The TrickleDownProfiler extracts channel networks from a landlab grid.
    Unlike the ChannelProfiler which starts at one or more watershed outlets
    and works upstream until it reaches the end of the channel (based on
    a specified threshold, such as drainage area) the TrickleDownProfiler
    starts at a *starting node* and works its way downhill until it reaches
    an outlet or sink.

    In order to follow the channel network, the flow connectivity across the
    grid must already be identified. This is typically done with the
    FlowAccumulator component. However, this component does not require that the
    FlowAccumulator was used. Instead it expects that the following at-node
    grid fields will be present:
    ::

        'flow__receiver_node'
        'flow__link_to_receiver_node'

    The TrickleDownProfiler can work on grids that have used route-to-one or
    route-to-multiple flow directing.

    To understand how this component works it is useful to define the following
    terms: *outlet*, *starting node*, and *segment*.

    Consider the following grid with 10 columns and 7 rows. ``@`` represents
    the *starting node*, ``.`` represents the nodes downstream, and the
    watershed outlet node is indicated by ``o``.

    In this and the following examples, we will use only D4 connectivity. The
    ChannelProfiler, however, knows nothing of connectivity other than what is
    implied by the two required grid fields.
    ::

        X X X X X X X X X X
        X X X X X X X X X X
        X X X X X X . . @ X
        X X X X X X . X X X
        X X X . . . . X X X
        X X X . X X X X X X
        X X X o X X X X X X

    For each starting node, the TrickleDownProfiler follows the network
    downstream until it reaches the outlet or sink. One or more starting nodes
    can be used, depending on a user's needs.

    The node IDs and distances upstream of the channel network are stored in
    ``data_structure``. It is a dictionary with keys indicating the starting
    node.

    For each starting node, the value in the ``data_structure`` is itself
    a dictionary with keys that are a segment ID tuple of the
    ``(dowstream, upstream)`` nodes IDs of each channel segment.

    For our simple example, these are the node IDs:
    ::

            X  X  X  X  X  X  X  X  X  X
            X  X  X  X  X  X  X  X  X  X
            X  X  X  X  X  X 46 47 48  X
            X  X  X  X  X  X 36  X  X  X
            X  X  X 23 24 25 26  X  X  X
            X  X  X 13  X  X  X  X  X  X
            X  X  X  3  X  X  X  X  X  X

    The starting node is 48 and the outlet node is 3.

    The value associated with the segment ID tuple ``(3, 48)`` is itself a
    dictionary. It has three key-value pairs. First, ``"ids"`` contains a list
    of the segment node ids ordered from downstream to upstream. It includes
    the endpoints. Second, ``"distances"`` contains a list of distances
    upstream that mirrors the list in ``"ids"``. Finally, ``"color"`` is an
    RGBA tuple indicating the color for the segment.

    By default a unique color will be assigned to each starting node. To change
    the color, a user can change values stored in ``data_structure``.
    Additionally, a ``cmap`` keyword argument can provide some user control
    over the color at the instantiation of the component.

    For example with a starting node of 48, the data structure will look as
    follows:

    .. code-block:: python

        {
            48: {
                (3, 48): {
                    "ids": [3, 13, 23, 24, 25, 26, 36, 46, 47, 48],
                    "distances": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                    "color": (1, 0, 1, 1),
                }
            }
        }

    Note that the distances upstream are relative to the outlet.

    Next consider an example with two starting nodes, each noted with an ``@``.
    ::

        X X X X X X X X X X
        X X X X X X @ X X X
        X X X X @ X . X X X
        o . . . . X . X X X
        X X X X X X . X X X
        X X X X X X . X X X
        X X X X X X . . . X
        X X X X X X X X o X

    And the following node IDs.
    ::

        X   X   X   X   X   X   X   X   X   X
        X   X   X   X   X   X  66   X   X   X
        X   X   X   X  54   X  56   X   X   X
       40  41  42  43  44   X  46   X   X   X
        X   X   X   X   X   X  36   X   X   X
        X   X   X   X   X   X  26   X   X   X
        X   X   X   X   X   X  16  17  18   X
        X   X   X   X   X   X   X   X   8   X

    With our starting nodes of 54 and 66 our data structure will look like.

    .. code-block:: python

        {
            54: {
                (40, 54): {
                    "ids": [40, 41, 42, 43, 44, 54],
                    "distances": [0, 1, 3, 4, 5, 6],
                    "color": [0.27, 0.0, 0.33, 1.0],
                },
            },
            66: {
                (8, 66): {
                    "ids": [8, 18, 17, 16, 26, 36, 46, 56, 66],
                    "distances": [0, 1, 2, 3, 4, 5, 6, 7, 8],
                    "color": [0.13, 0.57, 0.55, 1.0],
                },
            },
        }

    Examples
    --------

    Start by importing necessary modules

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, TrickleDownProfiler

    Create the second example grid we showed above. Note that in order to do
    this we need to enter the elevations starting from the lower left so the
    elevation order may seem upside-down. In addition, in this example,
    elevation is only provided along the profiles. The third line of code below
    sets all nodes with a value of zero to closed, such that these nodes are
    igored.
    >>> z = np.array(
    ...     [
    ...         [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    ...         [0, 0, 0, 0, 0, 0, 4, 3, 2, 0],
    ...         [0, 0, 0, 8, 7, 6, 5, 0, 0, 0],
    ...         [0, 0, 0, 0, 0, 0, 6, 0, 0, 0],
    ...         [1, 3, 4, 5, 6, 0, 7, 0, 0, 0],
    ...         [0, 4, 0, 0, 7, 0, 8, 0, 0, 0],
    ...         [0, 5, 6, 0, 0, 0, 9, 0, 0, 0],
    ...         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ...     ],
    ...     dtype=float,
    ... )

    >>> mg = RasterModelGrid((8, 10))
    >>> z = mg.add_field("topographic__elevation", z, at="node")
    >>> mg.set_nodata_nodes_to_closed(z, 0)
    >>> fa = FlowAccumulator(mg, flow_director="D4")
    >>> fa.run_one_step()
    >>> fa.node_drainage_area.reshape(mg.shape)
    array([[  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  11.,   0.],
           [  0.,   0.,   0.,   0.,   0.,   0.,   9.,  10.,  11.,   0.],
           [  0.,   0.,   0.,   1.,   2.,   3.,   8.,   0.,   0.,   0.],
           [  0.,   0.,   0.,   0.,   0.,   0.,   4.,   0.,   0.,   0.],
           [  8.,   8.,   4.,   3.,   2.,   0.,   3.,   0.,   0.,   0.],
           [  0.,   3.,   0.,   0.,   1.,   0.,   2.,   0.,   0.,   0.],
           [  0.,   2.,   1.,   0.,   0.,   0.,   1.,   0.,   0.,   0.],
           [  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.]])

    >>> profiler = TrickleDownProfiler(mg, starting_nodes=[54, 66])
    >>> profiler.run_one_step()

    The keys of the property ``data_structure`` are the IDs of the two outlet
    nodes.

    >>> profiler.data_structure.keys()
    odict_keys([54, 66])

    Within the data structure, the value at key 54, is a dictionary of the
    one segment, each specified by a ``(dowstream, upstream)`` tuple:

    >>> profiler.data_structure[54].keys()
    dict_keys([(40, 54)])

    The value of the segment between nodes 40 and 54 has the following
    components:

    >>> profiler.data_structure[54][(40, 54)]["ids"]
    array([40, 41, 42, 43, 44, 54])
    >>> profiler.data_structure[54][(40, 54)]["distances"]
    array([0., 1., 2., 3., 4., 5.])
    >>> np.round(profiler.data_structure[54][(40, 54)]["color"], decimals=2)
    array([0.27, 0.  , 0.33, 1.  ])

    The rest of the ``profile_structure`` encodes information about the second
    profile which starts at node 66.

    >>> profiler.data_structure[66].keys()
    dict_keys([(8, 66)])

    >>> profiler.data_structure[66][(8, 66)]["ids"]
    array([ 8, 18, 17, 16, 26, 36, 46, 56, 66])
    >>> profiler.data_structure[66][(8, 66)]["distances"]
    array([0., 1., 2., 3., 4., 5., 6., 7., 8.])
    >>> np.round(profiler.data_structure[66][(8, 66)]["color"], decimals=2)
    array([0.13, 0.57, 0.55, 1.  ])


    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "TrickleDownProfiler"

    _unit_agnostic = True

    _info = {
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
    }

    def __init__(
        self,
        grid,
        starting_nodes=None,
        cmap="viridis",
    ):
        """
        Parameters
        ----------
        grid : Landlab Model Grid instance
        starting_nodes : iterable
        cmap : str, optional
            A valid matplotlib cmap string. Default is "viridis".

        """
        super().__init__(grid)

        self._cmap = plt.colormaps[cmap]

        self._flow_receiver = grid.at_node["flow__receiver_node"]
        self._starting_nodes = starting_nodes

    @property
    def data_structure(self):
        """OrderedDict defining the trickle down network.

        The IDs and upstream distance of the channel network nodes are stored
        in ``data_structure``. It is a dictionary with keys of the outlet node
        ID.

        For each starting node, the value in the ``data_structure`` is
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
        return self._data_struct

    def _create_profile_structure(self):
        """Create the profile_IDs data structure for channel network.

        The bound attribute self._profile structure is the channel segment
        datastructure. Profile structure is a list of length
        starting_nodes. Each element of profile_structure is itself a
        list of length number of stream segments that drain to each of the
        starting nodes. Each stream segment list contains the node ids of a
        stream segment from downstream to upstream.
        """
        self._data_struct = OrderedDict()

        for i in self._starting_nodes:
            channel_segment = []
            current_node = i

            # march downstream
            while self._flow_receiver[current_node] != current_node:
                channel_segment.append(current_node)
                current_node = self._flow_receiver[current_node]
            channel_segment.append(current_node)

            channel_segment.reverse()
            segment_tuple = (current_node, i)
            self._data_struct[i] = {
                segment_tuple: {"ids": as_id_array(channel_segment)}
            }

        self._calculate_distances()
        self.assign_colors()
        self._create_flat_structures()

    def _create_flat_structures(self):
        """Create expected flattened structures for ids, distances, and colors."""
        self._nodes = []

        self._distance_along_profile = []
        self._colors = []

        for outlet_id in self._data_struct:
            seg_tuples = self._data_struct[outlet_id].keys()
            self._nodes.extend(
                [self._data_struct[outlet_id][seg]["ids"] for seg in seg_tuples]
            )
            self._distance_along_profile.extend(
                [self._data_struct[outlet_id][seg]["distances"] for seg in seg_tuples]
            )
            self._colors.extend(
                [self._data_struct[outlet_id][seg]["color"] for seg in seg_tuples]
            )

    def assign_colors(self, color_mapping=None):
        """Assign a unique color for each starting node.

        Parameters
        ----------
        color_mapping : str
            Color map name.
        """

        if color_mapping is None:
            num_watersheds = len(self._data_struct)
            norm = mpl.colors.Normalize(vmin=0, vmax=num_watersheds)
            mappable = cm.ScalarMappable(norm=norm, cmap=self._cmap)
            color_mapping = {
                outlet_id: mappable.to_rgba(idx)
                for idx, outlet_id in enumerate(self._data_struct)
            }

        for outlet_id in self._data_struct:
            for segment_tuple in self._data_struct[outlet_id]:
                self._data_struct[outlet_id][segment_tuple]["color"] = color_mapping[
                    outlet_id
                ]

    def _calculate_distances(self):
        """Get distances along the network data structure."""
        distance_upstream = calculate_flow__distance(self._grid)
        for outlet_id in self._data_struct:
            for segment_tuple in self._data_struct[outlet_id]:
                ids = self._data_struct[outlet_id][segment_tuple]["ids"]
                d = distance_upstream[ids]
                self._data_struct[outlet_id][segment_tuple]["distances"] = d
