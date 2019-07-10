#! /usr/env/python
r"""
Landlab's Continuous-Time Stochastic (CTS) cellular automata modeling package.

Overview
--------

A CellLab CTS model implements a particular type of cellular
automaton (CA): a continuous-time stochastic CA. The approach is based on that
of Narteau et al. (2002, 2009) and Rozier and Narteau (2014). Like a normal
CA, the domain consists of a lattice of cells, each of which has a discrete
state. Unlike a conventional CA, the updating process is stochastic, and takes
place in continuous rather than discrete time. Any given pair (or "doublet")
of adjacent cell states has a certain specified probability of transition to a
different pair of states. The transition probability is given in the form of an
average *transition rate*, :math:`\lambda` (with dimensions of 1/T); the actual
time of transition is a random variable drawn from an exponential probability
distribution with mean :math:`1/\lambda`.

Subclasses
----------

Landlab provides for several different lattice and connection types:

-  RasterCTS: regular raster grid with transitions between horizontal and
   vertical cell pairs
-  OrientedRasterCTS: like a RasterLCA, but different transition rates can
   be assigned to vertical and horizontal pairs. This property of
   orientation can be used, for example, to implement rules representing
   gravitational attraction, or flow of a fluid with a particular
   direction.
-  RasterD8CTS: like a RasterLCA, but includes diagonal as well as vertical
   and horizontal cell pairs.
-  OrientedRasterD8CTS: as above but orientation also matters.
-  HexCTS: hexagonal grid
-  OrientedHexCTS: hexagonal grid, with transition rates allowed to vary
   according to orientation.

Encoding of "states"
--------------------
As in any traditional cellular automaton model, a LandlabCellularAutomaton
contains a grid of cells ("nodes" in Landlab parlance), each of which is has a
discrete state. States are represented by integers (0, 1, ... N).

In addition, every active link has an *orientation code* and a *link state
code*. The orientation code represents the orientation of the link in space: is
it "vertical" (aligned with the y axis), "horizontal" (aligned with x), or in
some other orientation? The number of possible orientations depends on the
subclass. The base class has only one orientation code (0) (meaning
"orientation doesn't matter), but this is overridden in some of the subclasses.
For example, the OrientedRasterLCA has two orientation codes (0 and 1, for
vertical and horizontal), while the OrientedHexLCA has three (representing the
three axes in a hex-cell / triagonal grid).

Each active link also has a *link state code*. The *state* of a link refers to
its particular combination of nodes and its orientation. For example, link
state 1 refers to a link in which the tail-node has state 0, the head-node has
state 1, and the orientation code is 0. The number of possible link states is
equal to R N^2, where R is the number of orientations (1 to 3, depending on the
subclass) and N is the number of possible node states. The simplest possible
Landlab CA model would have just one orientation code and two possible cell
states, so that there are four unique link states. These would be represented
by the tuples of (tail-node state, head-node state, orientation) as follows::

    link state 0 = (0, 0, 0)
    link state 1 = (0, 1, 0)
    link state 2 = (1, 0, 0)
    link state 3 = (1, 1, 0)

Main data structures
--------------------
node_state : 1d array of int (x number of nodes in grid)
    Node-based grid of node-state codes. This is the grid of cell (sic) states.

link_state_dict : dictionary
    Keys are 3-element tuples that represent the cell-state pairs and
    orientation code for each possible link type; values are the corresponding
    link-state codes. Allows you to look up the link-state code corresponding
    to a particular pair of adjacent nodes with a particular orientation.

node_pair : list (x number of possible link states)
    List of 3-element tuples representing all the various link states. Allows
    you to look up the node states and orientation corresponding to a
    particular link-state ID.

priority_queue : PriorityQueue object containing event records
    Queue containing all future transition events, sorted by time of occurrence
    (from soonest to latest).

next_update : 1d array (x number of links)
    Time (in the future) at which the link will undergo its next transition.
    You might notice that the update time for every scheduled transition is
    also stored with each event in the event queue. Why store it twice?
    Because a scheduled event might be invalidated after the event has been
    scheduled (because another transition has changed one of a link's two
    nodes, for example). The way to tell whether a scheduled event is still
    valid is to compare its time with the corresponding transition time in the
    *next_update* array. If they are different, the event is discarded.

link_orientation : 1d array of int8 (x number of links)
    Orientation code for each link.

link_state : 1d array of int (x number of links)
    State code for each link.

n_trn : 1d array of int (x number of possible link states)
    Number of transitions ("trn" stands for "transition") from a given link
    state.

trn_to : 1d array of ints (x # transitions)
    Stores the link-state code(s) to which a particular transition ID can
    transition.

trn_rate : 1d array of floats (# transitions)
    Rate associated with each link-state transition.


Created GT Sep 2014, starting from link_cap.py.
"""
from __future__ import print_function

import numpy as np
import pylab as plt

import landlab
from landlab.ca.cfuncs import (
    PriorityQueue,
    get_next_event_new,
    push_transitions_to_event_queue,
    run_cts_new,
)

_NEVER = 1e50

_DEBUG = False

_CORE = landlab.grid.base.CORE_NODE


class Transition(object):

    """A transition from one state to another.

    Represents a transition from one state ("from_state") to another
    ("to_state") at a link. The transition probability is represented by a rate
    parameter "rate", with dimensions of 1/T. The probability distribution of
    time until the transition event occurs is exponentional with mean 1/rate.
    The optional name parameter allows the caller to assign a name to any given
    transition.

    Note that from_state and to_state can now be either integer IDs for the
    standardised ordering of the link states (as before), or tuples explicitly
    describing the node state at each end, and the orientation.
    Orientation is 0: horizontal, L-R; 1: vertical, bottom-top.
    For such a tuple, order is (left/bottom, right/top, orientation).

    Transition() constructor sets 3 required properties and 2 optional
    properties for a transition from one cell pair to another.

    Parameters
    ----------
    from_state : int
        Code for the starting state of the cell pair (link)
    to_state : int
        Code for the new state of the cell pair (link)
    rate : float
        Average rate at which this transition occurs (dimension of 1/time)
    name : string (optional)
        Name for this transition
    swap_properties : bool (optional)
        Flag: should properties be exchanged between the two cells?
    """

    def __init__(
        self,
        from_state,
        to_state,
        rate,
        name=None,
        swap_properties=False,
        prop_update_fn=None,
    ):
        """
        Transition() constructor sets 3 required properties and 2 optional
        properties for a transition from one cell pair to another.

        Parameters
        ----------
        from_state : int
            Code for the starting state of the cell pair (link)
        to_state : int
            Code for the new state of the cell pair (link)
        rate : float
            Average rate at which this transition occurs (dimension of 1/time)
        name : string (optional)
            Name for this transition
        swap_properties : bool (optional)
            Flag: should properties be exchanged between the two cells?
        """
        self.from_state = from_state
        self.to_state = to_state
        self.rate = rate
        self.name = name
        self.swap_properties = swap_properties
        self.prop_update_fn = prop_update_fn


class CAPlotter(object):

    """Handle display of a CellLab-CTS grid.

    CAPlotter() constructor keeps a reference to the CA model, and
    optionally a colormap to be used with plots.

    Parameters
    ----------
    ca : LandlabCellularAutomaton object
        Reference to a CA model
    cmap : Matplotlib colormap, optional
        Colormap to be used in plotting

    Examples
    --------
    >>> from landlab import RasterModelGrid, HexModelGrid
    >>> from landlab.ca.celllab_cts import Transition
    >>> from landlab.ca.raster_cts import RasterCTS
    >>> import numpy as np
    >>> grid = RasterModelGrid((3, 5))
    >>> nsd = {0 : 'zero', 1 : 'one'}
    >>> trn_list = []
    >>> trn_list.append(Transition((0, 1, 0), (1, 1, 0), 1.0))
    >>> ins = np.arange(15) % 2
    >>> ca = RasterCTS(grid, nsd, trn_list, ins)
    >>> cap = CAPlotter(ca)
    >>> cap.gridtype
    'rast'
    >>> cap._cmap.name
    'jet'

    >>> from landlab.ca.hex_cts import HexCTS
    >>> import matplotlib
    >>> grid = HexModelGrid(3, 3)
    >>> ins = np.zeros(grid.number_of_nodes, dtype=int)
    >>> ca = HexCTS(grid, nsd, trn_list, ins)
    >>> cap = CAPlotter(ca, cmap=matplotlib.cm.pink)
    >>> cap.gridtype
    'hex'
    >>> cap._cmap.name
    'pink'
    """

    def __init__(self, ca, cmap=None, **kwds):
        """
        CAPlotter() constructor keeps a reference to the CA model, and
        optionally a colormap to be used with plots.

        Parameters
        ----------
        ca : LandlabCellularAutomaton object
            Reference to a CA model
        cmap : Matplotlib colormap, optional
            Colormap to be used in plotting
        """
        import matplotlib

        # Set the colormap; default to matplotlib's "jet" colormap
        if cmap is None:
            self._cmap = matplotlib.cm.jet
        else:
            self._cmap = cmap

        # Keep a reference to the CA model
        self.ca = ca

        # Initialize the plot and remember the grid type
        plt.ion()
        plt.figure(1)
        if type(ca.grid) is landlab.grid.hex.HexModelGrid:
            self.gridtype = "hex"
        else:
            self.gridtype = "rast"

    def update_plot(self):
        """Plot the current node state grid."""
        plt.clf()
        if self.gridtype == "rast":
            nsr = self.ca.grid.node_vector_to_raster(self.ca.node_state)
            plt.imshow(nsr, interpolation="None", origin="lower", cmap=self._cmap)
        else:
            self.ca.grid.hexplot(self.ca.node_state, color_map=self._cmap)

        plt.draw()
        plt.pause(0.001)

    def finalize(self):
        """Wrap up plotting.

        Wrap up plotting by switching off interactive model and showing the
        plot.
        """
        plt.ioff()
        plt.show()


class CellLabCTSModel(object):

    """Link-type (or doublet-type) cellular automaton model.

    A CellLabCTSModel implements a link-type (or doublet-type) cellular
    automaton model. A link connects a pair of cells. Each cell has a state
    (represented by an integer code), and each link also has a state that is
    determined by the states of the cell pair.

    Parameters
    ----------
    model_grid : Landlab ModelGrid object
        Reference to the model's grid
    node_state_dict : dict
        Keys are node-state codes, values are the names associated with
        these codes
    transition_list : list of Transition objects
        List of all possible transitions in the model
    initial_node_states : array of ints (x number of nodes in grid)
        Starting values for node-state grid
    prop_data : array (x number of nodes in grid), optional
        Array of properties associated with each node/cell
    prop_reset_value : number or object, optional
        Default or initial value for a node/cell property (e.g., 0.0).
        Must be same type as *prop_data*.
    """

    def __init__(
        self,
        model_grid,
        node_state_dict,
        transition_list,
        initial_node_states,
        prop_data=None,
        prop_reset_value=None,
        seed=0,
    ):
        """Initialize the CA model.

        Parameters
        ----------
        model_grid : Landlab ModelGrid object
            Reference to the model's grid
        node_state_dict : dict
            Keys are node-state codes, values are the names associated with
            these codes
        transition_list : list of Transition objects
            List of all possible transitions in the model
        initial_node_states : array of ints (x number of nodes in grid)
            Starting values for node-state grid
        prop_data : array (x number of nodes in grid), optional
            Array of properties associated with each node/cell
        prop_reset_value : number or object, optional
            Default or initial value for a node/cell property (e.g., 0.0).
            Must be same type as *prop_data*.
        seed : int, optional
            Seed for random number generation.
        """

        # Keep a copy of the model grid
        self.grid = model_grid

        # Initialize random number generation
        np.random.seed(seed)

        # Create an array that knows which links are connected to a boundary
        # node
        self.bnd_lnk = np.zeros(self.grid.number_of_links, dtype=np.int8)
        for link_id in range(self.grid.number_of_links):
            if (
                self.grid.status_at_node[self.grid.node_at_link_tail[link_id]] != _CORE
                or self.grid.status_at_node[self.grid.node_at_link_head[link_id]]
                != _CORE
            ):
                self.bnd_lnk[link_id] = True

        # Set up the initial node-state grid
        self.set_node_state_grid(initial_node_states)

        # Current simulation time starts out at zero
        self.current_time = 0.0

        # Figure out how many states there are, and make sure the input data
        # are self consistent.
        #   There are 2 x (N^2) link states, where N is the number of node
        # states. For example, if there are just two node states, 0 and 1, then
        # the possible oriented link pairs are listed below:
        #   0-0 0-1 1-0 1-1  0 0 1 1
        #                    0 1 0 1
        self.num_node_states = len(node_state_dict)
        self.num_node_states_sq = self.num_node_states * self.num_node_states
        self.num_link_states = self.number_of_orientations * self.num_node_states_sq

        assert type(transition_list) is list, "transition_list must be a list!"
        assert transition_list, "Transition list must contain at least one transition"
        last_type = None
        for t in transition_list:
            # TODO: make orientation optional for cases where
            # self.number_of_orientations = 1
            if isinstance(t.from_state, tuple) and isinstance(t.to_state, tuple):
                this_type = tuple
            else:
                this_type = int

            if this_type is tuple:
                # added to allow from and to states to be tuples, not just ids
                for i in t.from_state[:-1]:
                    assert (
                        i < self.num_node_states
                    ), "Transition from_state out of range"
                for i in t.to_state[:-1]:
                    assert i < self.num_node_states, "Transition to_state out of range"
                assert (
                    t.from_state[-1] < self.number_of_orientations
                ), "Encoding for orientation in from_state must be < number of orientations."
                assert (
                    t.to_state[-1] < self.number_of_orientations
                ), "Encoding for orientation in to_state must be < number of orientations."
            else:
                assert (
                    t.from_state < self.num_link_states
                ), "Transition from_state out of range"
                assert (
                    t.to_state < self.num_link_states
                ), "Transition to_state out of range"

            assert (
                last_type == this_type or last_type is None
            ), "All transition types must be either int IDs, or all tuples."
            # this test to ensure all entries are either IDs, or tuples, not
            # mixed
            last_type = this_type

        # Create priority queue for events and next_update array for links
        self.next_update = self.grid.add_zeros("link", "next_update_time")
        self.priority_queue = PriorityQueue()
        self.next_trn_id = -np.ones(self.grid.number_of_links, dtype=np.int)

        # Assign link types from node types
        self.create_link_state_dict_and_pair_list()

        # DEJH adds: convert transition_list to IDs if necessary
        # This is the new part that allows Transition from_ and to_ types
        # to be specified either as ints, or as tuples.
        transition_list_as_ID = transition_list[:]
        if type(transition_list[0].from_state) == tuple:
            # (then they all are..., because of the assertions in __init__)
            for i in range(len(transition_list)):
                transition_list_as_ID[i].from_state = self.link_state_dict[
                    transition_list[i].from_state
                ]
                transition_list_as_ID[i].to_state = self.link_state_dict[
                    transition_list[i].to_state
                ]

        # Set up the information needed to determine the orientation of links
        # in the lattice. The default method just creates an array of zeros
        # (all orientations considered the same), but this will be overridden
        # in subclasses that do use orientation.
        self.setup_array_of_orientation_codes()

        # Using the grid of node states, figure out all the link states
        self.assign_link_states_from_node_types()

        # Create transition data for links
        self.setup_transition_data(transition_list_as_ID)

        # Put the various transitions on the event queue
        self.push_transitions_to_event_queue()

        # In order to keep track of cell "properties", we create an array of
        # indices that refer to locations in the caller's code where properties
        # are tracked.
        self.propid = np.arange(self.grid.number_of_nodes)
        if prop_data is None:
            self.prop_data = np.zeros(self.grid.number_of_nodes)
            self.prop_reset_value = 0.0
        else:
            self.prop_data = prop_data
            self.prop_reset_value = prop_reset_value

    def set_node_state_grid(self, node_states):
        """Set the grid of node-state codes to node_states.

        Sets the grid of node-state codes to node_states. Also checks
        to make sure node_states is in the proper format, which is to
        say, it's a Numpy array of the same length as the number of nodes in
        the grid.

        **Creates**:

        *  self.node_state : 1D array of ints (x number of nodes in grid)
           The node-state array

        Parameters
        ----------
        node_states : 1D array of ints (x number of nodes in grid)

        Notes
        -----
        The node-state array is attached to the grid as a field with the name
        'node_state'.
        """
        assert (
            type(node_states) is np.ndarray
        ), "initial_node_states must be a Numpy array"
        assert (
            len(node_states) == self.grid.number_of_nodes
        ), "length of initial_node_states must equal number of nodes in grid"
        self.grid.at_node["node_state"] = node_states
        self.node_state = node_states

    def create_link_state_dict_and_pair_list(self):
        """Create a dict of link-state to node-state.

        Creates a dictionary that can be used as a lookup table to find out
        which link state corresponds to a particular pair of node states. The
        dictionary keys are 3-element tuples, each of which represents the
        state of the TAIL node, the HEAD node, and the orientation of the link.
        The values are integer codes representing the link state numbers.

        Notes
        -----
        Performance note: making self.node_pair a tuple does not appear to
        change time to lookup values in update_node_states. Changing it to a
        2D array of int actually slows it down.
        """
        self.link_state_dict = {}
        self.node_pair = []
        k = 0
        for orientation in range(self.number_of_orientations):
            for tail_state in range(self.num_node_states):
                for head_state in range(self.num_node_states):
                    self.link_state_dict[(tail_state, head_state, orientation)] = k
                    self.node_pair.append((tail_state, head_state, orientation))
                    k += 1

    def setup_array_of_orientation_codes(self):
        """Create array of active link orientation codes.

        Creates and configures an array that contain the orientation code for
        each active link (and corresponding cell pair).

        **creates**:

        * ``self.link_orientation`` : 1D numpy array

        Notes
        -----

        The setup varies depending on the type of LCA. The default is
        non-oriented, in which case we just have an array of zeros. Subclasses
        will override this method to handle lattices in which orientation
        matters (for example, vertical vs. horizontal in an OrientedRasterLCA).
        """
        self.link_orientation = np.zeros(self.grid.number_of_links, dtype=np.int8)

    def assign_link_states_from_node_types(self):
        """Assign link-state code for each link.

        Takes lists/arrays of "tail" and "head" node IDs for each link, and a
        dictionary that associates pairs of node states (represented as a
        3-element tuple, comprising the TAIL state, FROM state, and
        orientation) to link states.

        **creates**:

        * ``self.link_state`` : 1D numpy array
        """
        self.link_state = np.zeros(self.grid.number_of_links, dtype=int)

        for i in self.grid.active_links:
            orientation = self.link_orientation[i]
            node_pair = (
                self.node_state[self.grid.node_at_link_tail[i]],
                self.node_state[self.grid.node_at_link_head[i]],
                orientation,
            )
            self.link_state[i] = self.link_state_dict[node_pair]

    def setup_transition_data(self, xn_list):
        """Create transition data arrays."""

        # First, create an array that stores the number of possible transitions
        # out of each state.
        n_xn = np.zeros(self.num_link_states, dtype=int)
        for xn in xn_list:
            n_xn[xn.from_state] += 1
        self.n_trn = np.zeros(self.num_link_states, dtype=np.int)

        # Now, create arrays to hold the "to state" and transition rate for each
        # transition. These arrays are dimensioned N x M where N is the number
        # of states, and M is the maximum number of transitions from a single
        # state (for example if state 3 could transition either to state 1 or
        # state 4, and the other states only had one or zero possible
        # transitions, then the maximum would be 2).
        max_transitions = np.max(n_xn)
        self.trn_id = np.zeros((self.num_link_states, max_transitions), dtype=np.int)
        num_transitions = len(xn_list)
        self.trn_to = np.zeros(num_transitions, dtype=np.int)
        self.trn_rate = np.zeros(num_transitions)
        self.trn_propswap = np.zeros(num_transitions, dtype=np.int8)
        self.trn_prop_update_fn = np.zeros(num_transitions, dtype=object)

        for trn in range(num_transitions):
            self.trn_to[trn] = xn_list[trn].to_state
            self.trn_rate[trn] = xn_list[trn].rate
            self.trn_propswap[trn] = xn_list[trn].swap_properties
            if xn_list[trn].prop_update_fn is not None:
                self.trn_prop_update_fn[trn] = xn_list[trn].prop_update_fn
                self._use_propswap_or_callback = True
            from_state = xn_list[trn].from_state
            self.trn_id[from_state, self.n_trn[from_state]] = trn
            self.n_trn[from_state] += 1

    def push_transitions_to_event_queue(self):
        """
        Initializes the event queue by creating transition events for each
        cell pair that has one or more potential transitions and pushing these
        onto the queue. Also records scheduled transition times in the
        self.next_update array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.ca.celllab_cts import Transition
        >>> from landlab.ca.oriented_raster_cts import OrientedRasterCTS
        >>> import numpy as np
        >>> grid = RasterModelGrid((3, 5))
        >>> nsd = {0 : 'zero', 1 : 'one'}
        >>> trn_list = []
        >>> trn_list.append(Transition((0, 1, 0), (1, 0, 0), 1.0))
        >>> trn_list.append(Transition((1, 0, 0), (0, 1, 0), 2.0))
        >>> trn_list.append(Transition((0, 1, 1), (1, 0, 1), 3.0))
        >>> trn_list.append(Transition((0, 1, 1), (1, 1, 1), 4.0))
        >>> ins = np.arange(15) % 2
        >>> cts = OrientedRasterCTS(grid, nsd, trn_list, ins)
        >>> ev0 = cts.priority_queue._queue[0]
        >>> np.round(100 * ev0[0])
        12.0
        >>> ev0[2]  # this is the link ID
        16
        >>> ev6 = cts.priority_queue._queue[6]
        >>> np.round(100 * ev6[0])
        27.0
        >>> ev6[2]  # this is the link ID
        6
        >>> cts.next_trn_id[ev0[2]]  # ID of the transition to occur at this link
        3
        >>> cts.next_trn_id[cts.grid.active_links]
        array([-1,  2, -1,  1,  0,  1,  0,  2, -1,  3])
        """
        push_transitions_to_event_queue(
            self.grid.number_of_active_links,
            self.grid.active_links,
            self.n_trn,
            self.link_state,
            self.trn_id,
            self.trn_rate,
            self.next_update,
            self.next_trn_id,
            self.priority_queue,
        )

    def update_link_state_new(self, link, new_link_state, current_time):
        """
        Implements a link transition by updating the current state of the link
        and (if appropriate) choosing the next transition event and pushing it
        on to the event queue.

        Parameters
        ----------
        link : int
            ID of the link to update
        new_link_state : int
            Code for the new state
        current_time : float
            Current time in simulation
        """

        # If the link connects to a boundary, we might have a different state
        # than the one we planned
        if self.bnd_lnk[link]:
            fns = self.node_state[self.grid.node_at_link_tail[link]]
            tns = self.node_state[self.grid.node_at_link_head[link]]
            orientation = self.link_orientation[link]
            new_link_state = int(
                orientation * self.num_node_states_sq + fns * self.num_node_states + tns
            )

        self.link_state[link] = new_link_state
        if self.n_trn[new_link_state] > 0:
            (event_time, trn_id) = get_next_event_new(
                link,
                new_link_state,
                current_time,
                self.n_trn,
                self.trn_id,
                self.trn_rate,
            )
            self.priority_queue.push(link, event_time)
            self.next_update[link] = event_time
            self.next_trn_id[link] = trn_id
        else:
            self.next_update[link] = _NEVER
            self.next_trn_id[link] = -1

    def update_component_data(self, new_node_state_array):
        """Update all component data.

        Call this method to update all data held by the component, if, for
        example, another component or boundary conditions modify the node
        statuses outside the component between run steps.

        This method updates all necessary properties, including both node and
        link states.

        *new_node_state_array* is the updated list of node states, which must
        still all be compatible with the state list originally supplied to
        this component.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.ca.celllab_cts import Transition
        >>> from landlab.ca.raster_cts import RasterCTS
        >>> import numpy as np
        >>> grid = RasterModelGrid((3, 5))
        >>> nsd = {0 : 'zero', 1 : 'one'}
        >>> trn_list = []
        >>> trn_list.append(Transition((0, 1, 0), (1, 1, 0), 1.0))
        >>> ins = np.zeros(15, dtype=np.int)
        >>> ca = RasterCTS(grid, nsd, trn_list, ins)
        >>> list(ca.node_state[6:9])
        [0, 0, 0]
        >>> list(ca.link_state[9:13])
        [0, 0, 0, 0]
        >>> len(ca.priority_queue._queue)  # there are no transitions
        0
        >>> nns = np.arange(15) % 2        # make a new node-state grid...
        >>> ca.update_component_data(nns)  # ...and assign it
        >>> list(ca.node_state[6:9])
        [0, 1, 0]
        >>> list(ca.link_state[9:13])
        [2, 1, 2, 1]
        >>> len(ca.priority_queue._queue)  # now there are 5 transitions
        5
        """
        self.set_node_state_grid(new_node_state_array)
        self.assign_link_states_from_node_types()
        self.push_transitions_to_event_queue()

    # @profile
    def run(
        self, run_to, node_state_grid=None, plot_each_transition=False, plotter=None
    ):
        """Run the model forward for a specified period of time.

        Parameters
        ----------
        run_to : float
            Time to run to, starting from self.current_time
        node_state_grid : 1D array of ints (x number of nodes) (optional)
            Node states (if given, replaces model's current node state grid)
        plot_each_transition : bool (optional)
            Option to display the grid after each transition
        plotter : CAPlotter object (optional)
            Needed if caller wants to plot after every transition

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.ca.celllab_cts import Transition
        >>> from landlab.ca.oriented_raster_cts import OrientedRasterCTS
        >>> import numpy as np
        >>> grid = RasterModelGrid((3, 5))
        >>> nsd = {0 : 'zero', 1 : 'one'}
        >>> trn_list = []
        >>> trn_list.append(Transition((0, 1, 0), (1, 0, 0), 1.0))
        >>> trn_list.append(Transition((1, 0, 0), (0, 1, 0), 2.0))
        >>> trn_list.append(Transition((0, 1, 1), (1, 0, 1), 3.0))
        >>> trn_list.append(Transition((0, 1, 1), (1, 1, 1), 4.0))
        >>> ins = np.arange(15) % 2
        >>> cts = OrientedRasterCTS(grid, nsd, trn_list, ins)
        """
        if node_state_grid is not None:
            self.set_node_state_grid(node_state_grid)

        self.current_time = run_cts_new(
            run_to,
            self.current_time,
            self.priority_queue,
            self.next_update,
            self.grid.node_at_link_tail,
            self.grid.node_at_link_head,
            self.node_state,
            self.next_trn_id,
            self.trn_to,
            self.grid.status_at_node,
            self.num_node_states,
            self.num_node_states_sq,
            self.bnd_lnk,
            self.link_orientation,
            self.link_state,
            self.n_trn,
            self.trn_id,
            self.trn_rate,
            self.grid.links_at_node,
            self.grid.active_link_dirs_at_node,
            self.trn_propswap,
            self.propid,
            self.prop_data,
            self.prop_reset_value,
            self.trn_prop_update_fn,
            self,
            plot_each_transition,
            plotter,
        )
