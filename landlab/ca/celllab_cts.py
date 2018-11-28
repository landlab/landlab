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
node_state : 1d array (x number of nodes in grid)
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

priority_queue : heap containing event records
    Queue containing all future transition events, sorted by time of occurrence
    (from soonest to latest).

next_update : 1d array (x number of active links)
    Time (in the future) at which the link will undergo its next transition.
    You might notice that the update time for every scheduled transition is
    also stored in each Event object in the event queue. Why store it twice?
    Because a scheduled event might be invalidated after the event has been
    scheduled (because another transition has changed one of a link's two
    nodes, for example). The way to tell whether a scheduled event is still
    valid is to compare its time with the corresponding transition time in the
    *next_update* array. If they are different, the event is discarded.

link_orientation : 1d array of ints (x number of active links)
    Orientation code for each link.

link_state : 1d array of ints (x number of active links)
    State code for each link.

n_xn : 1d array of ints (x number of possible link states)
    Number of transitions ("xn" stands for "transition") from a given link
    state.

xn_to : 2d array of ints (# possible link states x max. # transitions)
    Stores the link-state code(s) to which a particular link state can
    transition. "max. # transitions" means the maximum number of transitions
    from a single state. For example, if each link state is associated with one
    and only one transition, then the maximum is 1, but if there is at least
    one link state that can have either of two different transitions, then the
    maximum would be two.

xn_rate : 2d array of floats (# possible link states x max. # transitions)
    Rate associated with each link-state transition.


Created GT Sep 2014, starting from link_cap.py.
"""
from __future__ import print_function

import numpy as np
import pylab as plt

# X from _heapq import heappush
# X from _heapq import heappop
# X from _heapq import heapify
import landlab

# X _USE_CYTHON = False

_CYTEST = True

_RUN_NEW = True

_TESTING = True

if _TESTING:
    from .cfuncs import PriorityQueue

# Xif _USE_CYTHON:
#    from .cfuncs import (update_link_states_and_transitions,
#                         run_cts, run_cts_lean, PriorityQueue)

# TODO:
#   - REMOVE TEST CLAUSE (MAKE IT NORMAL BEHAVIOR)
#   - REMOVE OBSOLUTE, UNUSED METHODS
#   - RENAME "_new" METHODS (REMOVE THE "_new")
#   - UPDATE CALLING SYNTAX IN CELLLAB_CTS AND FRIENDS
if _CYTEST:
    from landlab.ca.cfuncs import (
        update_node_states,
        push_transitions_to_event_queue_new,
        do_transition_new,
        update_link_states_and_transitions_new,
        run_cts_new,
        get_next_event_new,
    )

_NEVER = 1e50

_DEBUG = False

_TEST = False

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


# Xif _USE_CYTHON:
#    from .cfuncs import Event
# else:
#    class Event(object):
#        """Transition event at a link.
#
#        Represents a transition event at a link. The transition occurs at a
#        given link and a given time, and it involves a transition into the
#        state xn_to (an integer code representing the new link state;
#        "xn" is shorthand for "transition").
#
#        The class overrides the __lt__ (less than operator) method so that when
#        Event() objects are placed in a PriorityQueue, the earliest event is
#        given the highest priority (i.e., placed at the top of the queue).
#
#        Event() constructor sets 3 required properties and one optional
#        property.
#
#        Parameters
#        ----------
#        time : float
#            Time at which the event is scheduled to occur
#        link : int
#            ID of the link at which event occurs
#        xn_to : int
#            New state to which this cell pair (link) will transition
#        propswap : bool (optional)
#            Flag: does this event involve an exchange of properties between
#            the two cells?
#
#        Examples
#        --------
#        >>> from landlab.ca.celllab_cts import Event
#        >>> e1 = Event( 10.0, 1, 2)
#        >>> e2 = Event( 2.0, 3, 1)
#        >>> e1 < e2
#        False
#        >>> e2 < e1
#        True
#        """
#
#        def __init__(self, time, link, xn_to, propswap=False, prop_update_fn=None):
#            """
#            Event() constructor sets 3 required properties and one optional
#            property.
#
#            Parameters
#            ----------
#            time : float
#                Time at which the event is scheduled to occur
#            link : int
#                ID of the link at which event occurs
#            xn_to : int
#                New state to which this cell pair (link) will transition
#            propswap : bool (optional)
#                Flag: does this event involve an exchange of properties between
#                the two cells?
#            """
#            self.time = time
#            self.link = link
#            self.xn_to = xn_to
#            self.propswap = propswap
#            self.prop_update_fn = prop_update_fn
#
#        def __lt__(self, other):
#            """
#            Overridden less-than operator: returns true if the event on the left
#            has an earlier scheduled time than the event on the right
#            """
#            return self.time < other.time


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
        # Are we calling this from a subclass __init__? If so, then the
        # variable self.number_of_orientations should already be defined.
        try:
            self.number_of_orientations == 1
        except AttributeError:
            # if self.number_of_orientations not already defined
            self.number_of_orientations = 1

        # Keep a copy of the model grid; remember how many active links in it
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
        # X self.event_queue = []
        # X heapify(self.event_queue)
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
        if _CYTEST:
            push_transitions_to_event_queue_new(
                len(self.grid.active_links),
                self.grid.active_links,
                self.n_trn,
                self.link_state,
                self.trn_id,
                self.trn_rate,
                self.next_update,
                self.next_trn_id,
                self.priority_queue,
            )
        elif _RUN_NEW:
            self.push_transitions_to_event_queue_new()
        else:
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

        # Determine and remember whether we will handle property swaps and/or
        # callbacks in this model.
        if np.amax(self.xn_propswap) > 0:
            self._use_propswap_or_callback = True
        else:
            self._use_propswap_or_callback = False

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

        if False and _DEBUG:
            print()
            print("create_link_state_dict_and_pair_list(): dict is:")
            print((self.link_state_dict))
            print("  and the pair list is:")
            print((self.node_pair))

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

        if False and _DEBUG:
            print()
            print("assign_link_states_from_node_types(): the link state array is:")
            print((self.link_state))

    def setup_transition_data(self, xn_list):
        """Create transition data arrays."""

        # First, create an array that stores the number of possible transitions
        # out of each state.
        self.n_xn = np.zeros(self.num_link_states, dtype=int)
        for xn in xn_list:
            self.n_xn[xn.from_state] += 1
        self.n_trn = np.zeros(self.num_link_states, dtype=np.int)

        # Now, create arrays to hold the "to state" and transition rate for each
        # transition. These arrays are dimensioned N x M where N is the number
        # of states, and M is the maximum number of transitions from a single
        # state (for example if state 3 could transition either to state 1 or
        # state 4, and the other states only had one or zero possible
        # transitions, then the maximum would be 2).
        max_transitions = np.max(self.n_xn)
        self.trn_id = np.zeros((self.num_link_states, max_transitions), dtype=np.int)
        num_transitions = len(xn_list)
        self.trn_to = np.zeros(num_transitions, dtype=np.int)
        self.trn_rate = np.zeros(num_transitions)
        self.trn_propswap = np.zeros(num_transitions, dtype=np.int8)
        self.trn_prop_update_fn = np.zeros(num_transitions, dtype=object)
        # OLD
        self.xn_to = np.zeros((self.num_link_states, max_transitions), dtype=int)
        self.xn_rate = np.zeros((self.num_link_states, max_transitions))
        self.xn_propswap = np.zeros(
            (self.num_link_states, max_transitions), dtype=np.int8
        )
        self.xn_prop_update_fn = np.zeros(
            (self.num_link_states, max_transitions), dtype=object
        )

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

        # OLD
        # Populate the "to" and "rate" arrays
        # reset this and then re-do (inefficient but should work)
        self.n_xn[:] = 0
        for xn in xn_list:
            from_state = xn.from_state
            self.xn_to[from_state][self.n_xn[from_state]] = xn.to_state
            self.xn_rate[from_state][self.n_xn[from_state]] = xn.rate
            self.xn_propswap[from_state][self.n_xn[from_state]] = xn.swap_properties
            if xn.prop_update_fn is not None:
                self.xn_prop_update_fn[from_state][
                    self.n_xn[from_state]
                ] = xn.prop_update_fn
                self._use_propswap_or_callback = True
            self.n_xn[from_state] += 1

        if False and _DEBUG:
            print()
            print("setup_transition_data():")
            print(("  n_xn", self.n_xn))
            print(("  to:", self.xn_to))
            print(("  rate:", self.xn_rate))

    def current_link_state(self, link_id):
        """Get the current state of a link.

        Used to determine whether the link state at link *link_id* has changed
        due to an independent change in the node-state grid. Returns the
        current state of the link based on the states of its two end nodes;
        this can be compared to the entry in self.link_state to determine
        whether the state has changed.

        Parameters
        ----------
        link_id : int
            ID of the active link to test

        Returns
        -------
        int
            New link state code

        Notes
        -----
        Vectorizing this might yield some speed.
        """

        # Find out the states of the two nodes, and the orientation
        tail_node_state = self.node_state[self.grid.node_at_link_tail[link_id]]
        head_node_state = self.node_state[self.grid.node_at_link_head[link_id]]
        orientation = self.link_orientation[link_id]

        # Return the corresponding state code.
        return (
            orientation * self.num_node_states_sq
            + tail_node_state * self.num_node_states
            + head_node_state
        )

    def update_link_states_and_transitions(self, current_time):
        """
        Following an "external" change to the node state grid, updates link
        states where necessary and creates any needed events.

        Notes
        -----
        **Algorithm**::

            FOR each active link:
                if the actual node pair is different from the link's code:
                    change the link state to be correct
                    schedule an event
        """
        if _CYTEST:
            update_link_states_and_transitions_new(
                self.grid.active_links,
                self.node_state,
                self.grid.node_at_link_tail,
                self.grid.node_at_link_head,
                self.link_orientation,
                self.bnd_lnk,
                self.link_state,
                self.n_trn,
                self.priority_queue,
                self.next_update,
                self.next_trn_id,
                self.trn_id,
                self.trn_rate,
                self.num_node_states,
                self.num_node_states_sq,
                current_time,
            )
        # X       elif _USE_CYTHON:
        #            update_link_states_and_transitions(self.grid.active_links,
        #                                               self.node_state,
        #                                               self.grid.node_at_link_tail,
        #                                               self.grid.node_at_link_head,
        #                                               self.link_orientation,
        #                                               self.bnd_lnk,
        #                                               self.link_state,
        #                                               self.n_xn,
        #                                               self.event_queue,
        #                                               self.next_update,
        #                                               self.xn_to, self.xn_rate,
        #                                               self.num_node_states,
        #                                               self.num_node_states_sq,
        #                                               current_time,
        #                                               self.xn_propswap,
        #                                               self.xn_prop_update_fn)
        else:
            for i in self.grid.active_links:
                current_state = self.current_link_state(i)
                if current_state != self.link_state[i]:
                    self.update_link_state(i, current_state, current_time)

    def update_link_states_and_transitions_new(self, current_time):
        """
        Following an "external" change to the node state grid, updates link
        states where necessary and creates any needed events.

        Notes
        -----
        **Algorithm**::

            FOR each active link:
                if the actual node pair is different from the link's code:
                    change the link state to be correct
                    schedule an event
        """
        if _CYTEST:
            update_link_states_and_transitions_new(
                self.grid.active_links,
                self.node_state,
                self.grid.node_at_link_tail,
                self.grid.node_at_link_head,
                self.link_orientation,
                self.bnd_lnk,
                self.link_state,
                self.n_trn,
                self.priority_queue,
                self.next_update,
                self.next_trn_id,
                self.trn_id,
                self.trn_rate,
                self.num_node_states,
                self.num_node_states_sq,
                current_time,
            )
        else:
            for i in self.grid.active_links:
                current_state = self.current_link_state(i)
                if current_state != self.link_state[i]:
                    self.update_link_state_new(i, current_state, current_time)

    # X    def get_next_event(self, link, current_state, current_time):
    #        """Get the next event for a link.
    #
    #        Returns the next event for link with ID "link", which is in state
    #        "current state".
    #
    #        Parameters
    #        ----------
    #        link : int
    #            ID of the link
    #        current_state : int
    #            Current state code for the link
    #        current_time : float
    #            Current time in simulation (i.e., time of event just processed)
    #
    #        Returns
    #        -------
    #        Event object
    #            The returned Event object contains the time, link ID, and type of
    #            the next transition event at this link.
    #
    #        Notes
    #        -----
    #        If there is only one potential transition out of the current state, a
    #        time for the transition is selected at random from an exponential
    #        distribution with rate parameter appropriate for this transition.
    #
    #        If there are more than one potential transitions, a transition time is
    #        chosen for each, and the smallest of these applied.
    #
    #        Assumes that there is at least one potential transition from the
    #        current state.
    #        """
    #        assert (self.n_xn[current_state] > 0), \
    #            'must have at least one potential transition'
    #
    #        # Find next event time for each potential transition
    #        if self.n_xn[current_state] == 1:
    #            xn_to = self.xn_to[current_state][0]
    #            propswap = self.xn_propswap[current_state][0]
    #            next_time = np.random.exponential(
    #                1.0 / self.xn_rate[current_state][0])
    #            if _DEBUG:
    #                print('next_time for 1 xn is ' + str(next_time))
    #            prop_update_fn = self.xn_prop_update_fn[current_state][0]
    #        else:
    #            next_time = _NEVER
    #            xn_to = None
    #            propswap = False
    #            for i in range(self.n_xn[current_state]):
    #                this_next = np.random.exponential(
    #                    1.0 / self.xn_rate[current_state][i])
    #                if _DEBUG:
    #                    print('this_next for >1 xn is ' + str(this_next))
    #                if this_next < next_time:
    #                    next_time = this_next
    #                    xn_to = self.xn_to[current_state][i]
    #                    propswap = self.xn_propswap[current_state][i]
    #                    prop_update_fn = self.xn_prop_update_fn[current_state][i]
    #
    #        # Create and setup event, and return it
    #        my_event = Event(next_time + current_time, link,
    #                         xn_to, propswap, prop_update_fn)
    #
    #        if _DEBUG:
    #            print('get_next_event():')
    #            print(('  current state:', current_state))
    #            print(('  node pair:', self.node_pair[current_state]))
    #            print(('  next_time:', my_event.time))
    #            print(('  link:', my_event.link))
    #            print(('  xn_to:', my_event.xn_to))
    #            print(('  nxn:', self.n_xn[current_state]))
    #            print(('  rate:', self.xn_rate[current_state][:]))
    #            print(('  propswap:', my_event.propswap))
    #
    #        return my_event

    def get_next_event_new(self, link, current_state, current_time):
        """Get the next event for a link.

        Returns the next event for link with ID "link", which is in state
        "current state".

        Parameters
        ----------
        link : int
            ID of the link
        current_state : int
            Current state code for the link
        current_time : float
            Current time in simulation (i.e., time of event just processed)

        Returns
        -------
        Event object
            The returned Event object contains the time, link ID, and type of
            the next transition event at this link.

        Notes
        -----
        If there is only one potential transition out of the current state, a
        time for the transition is selected at random from an exponential
        distribution with rate parameter appropriate for this transition.

        If there are more than one potential transitions, a transition time is
        chosen for each, and the smallest of these applied.

        Assumes that there is at least one potential transition from the
        current state.
        """
        assert (
            self.n_xn[current_state] > 0
        ), "must have at least one potential transition"

        # Find next event time for each potential transition: new version
        if self.n_trn[current_state] == 1:
            trn_id = self.trn_id[current_state, 0]
            next_time = np.random.exponential(1.0 / self.trn_rate[trn_id])
            if _DEBUG:
                print("rand (" + str(self.trn_rate[trn_id]) + ": " + str(next_time))
        else:
            next_time = _NEVER
            trn_id = -1
            for i in range(self.n_trn[current_state]):
                this_next = np.random.exponential(
                    1.0 / self.trn_rate[self.trn_id[current_state][i]]
                )
                if _DEBUG:
                    print(
                        "rand2 ("
                        + str(self.trn_rate[self.trn_id[current_state][i]])
                        + ": "
                        + str(this_next)
                    )
                if this_next < next_time:
                    if _DEBUG:
                        print(" using it")
                    next_time = this_next
                    trn_id = self.trn_id[current_state, i]

        if _DEBUG:
            print("get_next_event_new():")
            print(("  current state:", current_state))
            print(("  node pair:", self.node_pair[current_state]))
            print(("  next_time:", next_time))
            print(("  link:", link))
            print(("  trn_id:", trn_id))

        return (next_time + current_time, trn_id)

    # X    def push_transitions_to_event_queue(self):
    #        """
    #        Initializes the event queue by creating transition events for each
    #        cell pair that has one or more potential transitions and pushing these
    #        onto the queue. Also records scheduled transition times in the
    #        self.next_update array.
    #
    #        Examples
    #        --------
    #        >>> from landlab import RasterModelGrid
    #        >>> from landlab.ca.celllab_cts import Transition
    #        >>> from landlab.ca.oriented_raster_cts import OrientedRasterCTS
    #        >>> import numpy as np
    #        >>> grid = RasterModelGrid((3, 5))
    #        >>> nsd = {0 : 'zero', 1 : 'one'}
    #        >>> trn_list = []
    #        >>> trn_list.append(Transition((0, 1, 0), (1, 0, 0), 1.0))
    #        >>> trn_list.append(Transition((1, 0, 0), (0, 1, 0), 2.0))
    #        >>> trn_list.append(Transition((0, 1, 1), (1, 0, 1), 3.0))
    #        >>> trn_list.append(Transition((0, 1, 1), (1, 1, 1), 4.0))
    #        >>> ins = np.arange(15) % 2
    #        >>> cts = OrientedRasterCTS(grid, nsd, trn_list, ins)
    #        """
    #        if False and _DEBUG:
    #            print(('push_transitions_to_event_queue():',
    #                   self.num_link_states, self.n_xn))
    #
    #        for i in self.grid.active_links:
    #
    #            if self.n_xn[self.link_state[i]] > 0:
    #                event = self.get_next_event(i, self.link_state[i], 0.0)
    #                heappush(self.event_queue, event)
    #                self.next_update[i] = event.time
    #
    #            else:
    #                self.next_update[i] = _NEVER
    #
    #        if False and _DEBUG:
    #            print('  push_transitions_to_event_queue(): events in queue are now:')
    #            for e in self.event_queue:
    #                print('    next_time:', e.time, 'link:',
    #                      e.link, 'xn_to:', e.xn_to)

    def push_transitions_to_event_queue_new(self):
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
        if False and _DEBUG:
            print(
                (
                    "push_transitions_to_event_queue_new():",
                    self.num_link_states,
                    self.n_xn,
                )
            )

        for i in self.grid.active_links:

            if self.n_xn[self.link_state[i]] > 0:
                (ev_time, trn_id) = self.get_next_event_new(i, self.link_state[i], 0.0)
                self.priority_queue.push(i, ev_time)
                self.next_update[i] = ev_time
                self.next_trn_id[i] = trn_id

            else:
                self.next_update[i] = _NEVER

    # @profile
    def update_node_states(self, tail_node, head_node, new_link_state):
        """Update the states of the two nodes in the given link.

        Parameters
        ----------
        tail_node : int
            ID of the tail node of the link (cell pair) in question
        head_node : int
            ID of the head node of the link (cell pair) in question
        new_link_state : int
            Link state code for the new cell pair

        Returns
        -------
        (bool, bool)
            Flags indicating whether the tail node and head node, respectively,
            have changed state
        """
        # Remember the previous state of each node so we can detect whether the
        # state has changed
        old_tail_node_state = self.node_state[tail_node]
        old_head_node_state = self.node_state[head_node]

        # Change to the new states
        if self.grid.status_at_node[tail_node] == _CORE:
            self.node_state[tail_node] = self.node_pair[new_link_state][0]
        if self.grid.status_at_node[head_node] == _CORE:
            self.node_state[head_node] = self.node_pair[new_link_state][1]

        return (
            self.node_state[tail_node] != old_tail_node_state,
            self.node_state[head_node] != old_head_node_state,
        )

    # X    def update_link_state(self, link, new_link_state, current_time):
    #        """
    #        Implements a link transition by updating the current state of the link
    #        and (if appropriate) choosing the next transition event and pushing it
    #        on to the event queue.
    #
    #        Parameters
    #        ----------
    #        link : int
    #            ID of the link to update
    #        new_link_state : int
    #            Code for the new state
    #        current_time : float
    #            Current time in simulation
    #        """
    #        if _DEBUG:
    #            print('update_link_state() link ' + str(link) + ' to state ' + str(new_link_state) + ' at time ' + str(current_time))
    #
    #        # If the link connects to a boundary, we might have a different state
    #        # than the one we planned
    #        if self.bnd_lnk[link]:
    #            fns = self.node_state[self.grid.node_at_link_tail[link]]
    #            tns = self.node_state[self.grid.node_at_link_head[link]]
    #            orientation = self.link_orientation[link]
    #            new_link_state = orientation * self.num_node_states_sq + \
    #                fns * self.num_node_states + tns
    #
    #        self.link_state[link] = new_link_state
    #        if self.n_xn[new_link_state] > 0:
    #            event = self.get_next_event(link, new_link_state, current_time)
    #            heappush(self.event_queue, event)
    #            self.next_update[link] = event.time
    #        else:
    #            self.next_update[link] = _NEVER

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
        if _DEBUG:
            print(
                "update_link_state_new() link "
                + str(link)
                + " to state "
                + str(new_link_state)
                + " at time "
                + str(current_time)
            )

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
            if _CYTEST:
                (event_time, trn_id) = get_next_event_new(
                    link,
                    new_link_state,
                    current_time,
                    self.n_trn,
                    self.trn_id,
                    self.trn_rate,
                )
            else:
                (event_time, trn_id) = self.get_next_event_new(
                    link, new_link_state, current_time
                )
            self.priority_queue.push(link, event_time)
            if _DEBUG:
                print(
                    "Pushed event at "
                    + str(link)
                    + " for time "
                    + str(event_time)
                    + " id "
                    + str(trn_id)
                )
            self.next_update[link] = event_time
            self.next_trn_id[link] = trn_id
        else:
            self.next_update[link] = _NEVER
            self.next_trn_id[link] = -1

    def do_transition(
        self, event, current_time, plot_each_transition=False, plotter=None
    ):
        """Transition state.

        Implements a state transition.

        Parameters
        ----------
        event : Event object
            Event object containing the data for the current transition event
        current_time : float
            Current time in simulation
        plot_each_transition : bool (optional)
            True if caller wants to show a plot of the grid after this
            transition
        plotter : CAPlotter object
            Sent if caller wants a plot after this transition

        Notes
        -----
        First checks that the transition is still valid by comparing the
        link's next_update time with the corresponding update time in the
        event object.

        If the transition is valid, we:

        1. Update the states of the two nodes attached to the link
        2. Update the link's state, choose its next transition, and push
           it on the event queue.
        3. Update the states of the other links attached to the two nodes,
           choose their next transitions, and push them on the event queue.
        """
        if _DEBUG:
            print()
            print(
                "py do_transition() for link",
                event.link,
                "time",
                event.time,
                " cur time ",
                current_time,
            )
            print(np.count_nonzero(self.node_state != 0))
            print(self.link_state[event.link])

        # We'll process the event if its update time matches the one we have
        # recorded for the link in question. If not, it means that the link has
        # changed state since the event was pushed onto the event queue, and
        # in that case we'll ignore it.
        if event.time == self.next_update[event.link]:

            if _DEBUG:
                print("  event time =", event.time)

            tail_node = self.grid.node_at_link_tail[event.link]
            head_node = self.grid.node_at_link_head[event.link]
            tail_changed, head_changed = self.update_node_states(
                tail_node, head_node, event.xn_to
            )
            self.update_link_state(event.link, event.xn_to, event.time)

            # Next, when the state of one of the link's nodes changes, we have
            # to update the states of the OTHER links attached to it. This
            # could happen to one or both nodes.
            if tail_changed:

                if _DEBUG:
                    print(" tail node has changed state, so updating its links")
                    print(" links at node " + str(tail_node) + " are:")
                    print(self.grid.links_at_node[tail_node, :])
                    print(self.grid.active_link_dirs_at_node[tail_node, :])

                for i in range(self.grid.links_at_node.shape[1]):

                    link = self.grid.links_at_node[tail_node, i]
                    dir_code = self.grid.active_link_dirs_at_node[tail_node, i]
                    if _DEBUG:
                        print("tail checking link", link)
                        print(
                            "  dir code "
                            + str(dir_code)
                            + " event link "
                            + str(event.link)
                        )
                    if dir_code != 0 and link != event.link:

                        this_link_fromnode = self.grid.node_at_link_tail[link]
                        this_link_tonode = self.grid.node_at_link_head[link]
                        orientation = self.link_orientation[link]
                        current_pair = (
                            self.node_state[this_link_fromnode],
                            self.node_state[this_link_tonode],
                            orientation,
                        )
                        new_link_state = self.link_state_dict[current_pair]
                        new_link_state2 = (
                            orientation * self.num_node_states_sq
                            + self.node_state[this_link_fromnode] * self.num_node_states
                            + self.node_state[this_link_tonode]
                        )
                        assert new_link_state == new_link_state2, "oops"
                        self.update_link_state(link, new_link_state, event.time)

            if head_changed:

                if _DEBUG:
                    print(" head node has changed state, so updating its links")

                for i in range(self.grid.links_at_node.shape[1]):

                    link = self.grid.links_at_node[head_node, i]
                    dir_code = self.grid.active_link_dirs_at_node[head_node, i]
                    if _DEBUG:
                        print("head checking link", link)
                    if dir_code != 0 and link != event.link:
                        this_link_fromnode = self.grid.node_at_link_tail[link]
                        this_link_tonode = self.grid.node_at_link_head[link]
                        orientation = self.link_orientation[link]
                        current_pair = (
                            self.node_state[this_link_fromnode],
                            self.node_state[this_link_tonode],
                            orientation,
                        )
                        new_link_state = self.link_state_dict[current_pair]
                        new_link_state2 = (
                            orientation * self.num_node_states_sq
                            + self.node_state[this_link_fromnode] * self.num_node_states
                            + self.node_state[this_link_tonode]
                        )
                        assert new_link_state == new_link_state2, "oops"
                        self.update_link_state(link, new_link_state, event.time)

            # If requested, display a plot of the grid
            if plot_each_transition and (plotter is not None):
                plotter.update_plot()

            # If this event involves an exchange of properties (i.e., the
            # event involves motion of an object that posses properties we
            # want to track), implement the swap.
            #   If the event requires a call to a user-defined callback
            # function, we handle that here too.
            if event.propswap:
                tmp = self.propid[tail_node]
                self.propid[tail_node] = self.propid[head_node]
                self.propid[head_node] = tmp
                if self.grid.status_at_node[tail_node] != _CORE:
                    self.prop_data[self.propid[tail_node]] = self.prop_reset_value
                if self.grid.status_at_node[head_node] != _CORE:
                    self.prop_data[self.propid[head_node]] = self.prop_reset_value
                if event.prop_update_fn is not None:
                    event.prop_update_fn(self, tail_node, head_node, event.time)

        elif _DEBUG:
            print(
                "  event time is",
                event.time,
                "but update time is",
                self.next_update[event.link],
                "so event will be ignored",
            )

    def do_transition_new(
        self,
        event_link,
        event_time,
        current_time,
        plot_each_transition=False,
        plotter=None,
    ):
        """Transition state.

        Implements a state transition.

        Parameters
        ----------
        event : Event object
            Event object containing the data for the current transition event
        current_time : float
            Current time in simulation
        plot_each_transition : bool (optional)
            True if caller wants to show a plot of the grid after this
            transition
        plotter : CAPlotter object
            Sent if caller wants a plot after this transition

        Notes
        -----
        First checks that the transition is still valid by comparing the
        link's next_update time with the corresponding update time in the
        event object.

        If the transition is valid, we:

        1. Update the states of the two nodes attached to the link
        2. Update the link's state, choose its next transition, and push
           it on the event queue.
        3. Update the states of the other links attached to the two nodes,
           choose their next transitions, and push them on the event queue.

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
        >>> (tm, idx, link) = cts.priority_queue.pop()
        >>> np.round(100 * tm)
        12.0
        >>> idx
        6
        >>> link
        16
        >>> cts.grid.node_at_link_tail[link]
        8
        >>> cts.grid.node_at_link_head[link]
        13
        >>> cts.next_trn_id[link]
        3
        >>> cts.trn_to[cts.next_trn_id[link]]
        7
        >>> cts.do_transition_new(link, tm, 0.0)
        >>> cts.node_state[8]
        1
        >>> cts.node_state[13]
        1
        >>> cts.link_state[16]
        7
        >>> cts.next_update[16] == _NEVER
        True
        >>> cts.next_trn_id[16]
        -1
        """
        if _DEBUG:
            print()
            print(
                "py do_transition_new() for link",
                event_link,
                "time",
                event_time,
                " cur time ",
                current_time,
            )

        # We'll process the event if its update time matches the one we have
        # recorded for the link in question. If not, it means that the link has
        # changed state since the event was pushed onto the event queue, and
        # in that case we'll ignore it.
        if event_time == self.next_update[event_link]:

            if _DEBUG:
                print("  event time =", event_time)

            tail_node = self.grid.node_at_link_tail[event_link]
            head_node = self.grid.node_at_link_head[event_link]
            trn_id = self.next_trn_id[event_link]
            trn_to = self.trn_to[trn_id]
            if _CYTEST:
                old_tail_node_state = self.node_state[tail_node]
                old_head_node_state = self.node_state[head_node]
                update_node_states(
                    self.node_state,
                    self.grid.status_at_node,
                    tail_node,
                    head_node,
                    trn_to,
                    self.num_node_states,
                )
                tail_changed = old_tail_node_state != self.node_state[tail_node]
                head_changed = old_head_node_state != self.node_state[head_node]
            else:
                tail_changed, head_changed = self.update_node_states(
                    tail_node, head_node, trn_to
                )
            self.update_link_state_new(event_link, trn_to, event_time)

            # Next, when the state of one of the link's nodes changes, we have
            # to update the states of the OTHER links attached to it. This
            # could happen to one or both nodes.
            if tail_changed:

                if _DEBUG:
                    print(" tail node has changed state, so updating its links")
                    print(" links at node " + str(tail_node) + " are:")
                    print(self.grid.links_at_node[tail_node, :])
                    print(self.grid.active_link_dirs_at_node[tail_node, :])

                for i in range(self.grid.links_at_node.shape[1]):

                    link = self.grid.links_at_node[tail_node, i]
                    dir_code = self.grid.active_link_dirs_at_node[tail_node, i]
                    if _DEBUG:
                        print("tail checking link", link)
                        print(
                            "  dir code "
                            + str(dir_code)
                            + " event link "
                            + str(event_link)
                        )
                    if dir_code != 0 and link != event_link:

                        this_link_fromnode = self.grid.node_at_link_tail[link]
                        this_link_tonode = self.grid.node_at_link_head[link]
                        orientation = self.link_orientation[link]
                        current_pair = (
                            self.node_state[this_link_fromnode],
                            self.node_state[this_link_tonode],
                            orientation,
                        )
                        new_link_state = self.link_state_dict[current_pair]
                        new_link_state2 = (
                            orientation * self.num_node_states_sq
                            + self.node_state[this_link_fromnode] * self.num_node_states
                            + self.node_state[this_link_tonode]
                        )
                        assert new_link_state == new_link_state2, "oops"
                        self.update_link_state_new(link, new_link_state, event_time)

            if head_changed:

                if _DEBUG:
                    print(" head node has changed state, so updating its links")

                for i in range(self.grid.links_at_node.shape[1]):

                    link = self.grid.links_at_node[head_node, i]
                    dir_code = self.grid.active_link_dirs_at_node[head_node, i]
                    if _DEBUG:
                        print("head checking link", link)
                    if dir_code != 0 and link != event_link:
                        this_link_fromnode = self.grid.node_at_link_tail[link]
                        this_link_tonode = self.grid.node_at_link_head[link]
                        orientation = self.link_orientation[link]
                        current_pair = (
                            self.node_state[this_link_fromnode],
                            self.node_state[this_link_tonode],
                            orientation,
                        )
                        new_link_state = self.link_state_dict[current_pair]
                        new_link_state2 = (
                            orientation * self.num_node_states_sq
                            + self.node_state[this_link_fromnode] * self.num_node_states
                            + self.node_state[this_link_tonode]
                        )
                        assert new_link_state == new_link_state2, "oops"
                        self.update_link_state_new(link, new_link_state, event_time)

            # If requested, display a plot of the grid
            if plot_each_transition and (plotter is not None):
                plotter.update_plot()

            # If this event involves an exchange of properties (i.e., the
            # event involves motion of an object that posses properties we
            # want to track), implement the swap.
            #   If the event requires a call to a user-defined callback
            # function, we handle that here too.
            if self.trn_propswap[trn_id]:
                tmp = self.propid[tail_node]
                self.propid[tail_node] = self.propid[head_node]
                self.propid[head_node] = tmp
                if self.grid.status_at_node[tail_node] != _CORE:
                    self.prop_data[self.propid[tail_node]] = self.prop_reset_value
                if self.grid.status_at_node[head_node] != _CORE:
                    self.prop_data[self.propid[head_node]] = self.prop_reset_value
                if self.trn_prop_update_fn[trn_id] is not None:
                    self.trn_prop_update_fn[trn_id](
                        self, tail_node, head_node, event_time
                    )

            if False and _DEBUG:
                n = self.grid.number_of_nodes
                for r in range(self.grid.number_of_node_rows):
                    for c in range(self.grid.number_of_node_columns):
                        n -= 1
                        print("{0:.0f}".format(self.node_state[n]), end=" ")
                    print()
                if self.propid is not None:
                    print()
                    n = self.grid.number_of_nodes
                    for r in range(self.grid.number_of_node_rows):
                        for c in range(self.grid.number_of_node_columns):
                            n -= 1
                            print("{0:2.0f}".format(self.propid[n]), end=" ")
                        print()

        elif _DEBUG:
            print(
                "  event time is",
                event_time,
                "but update time is",
                self.next_update[event_link],
                "so event will be ignored",
            )

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
        """
        self.set_node_state_grid(new_node_state_array)
        self.assign_link_states_from_node_types()
        self.push_transitions_to_event_queue_new()

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

        # X        if plot_each_transition or self._use_propswap_or_callback:
        #            lean_run = False
        #        else:
        #            lean_run = True

        # X        if _USE_CYTHON and not lean_run:
        #            self.current_time = run_cts(run_to, self.current_time,
        #               plot_each_transition,
        #               plotter,
        #               self.event_queue,
        #               self.next_update,
        #               self.grid.node_at_link_tail,
        #               self.grid.node_at_link_head,
        #               self.node_state,
        #               self.link_state,
        #               self.grid.status_at_node,
        #               self.link_orientation,
        #               self.propid,
        #               self.prop_data,
        #               self.n_xn,
        #               self.xn_to,
        #               self.xn_rate,
        #               self.grid.links_at_node,
        #               self.grid.active_link_dirs_at_node,
        #               self.num_node_states,
        #               self.num_node_states_sq,
        #               self.prop_reset_value,
        #               self.xn_propswap,
        #               self.xn_prop_update_fn,
        #               self.bnd_lnk,
        #               self)
        #        elif _USE_CYTHON and lean_run:
        #            self.current_time = run_cts_lean(run_to, self.current_time,
        #               self.event_queue,
        #               self.next_update,
        #               self.grid.node_at_link_tail,
        #               self.grid.node_at_link_head,
        #               self.node_state,
        #               self.link_state,
        #               self.grid.status_at_node,
        #               self.link_orientation,
        #               self.n_xn,
        #               self.xn_to,
        #               self.xn_rate,
        #               self.grid.links_at_node,
        #               self.grid.active_link_dirs_at_node,
        #               self.num_node_states,
        #               self.num_node_states_sq,
        #               self.bnd_lnk)

        #         if _RUN_NEW:

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

    #         else:
    #
    #            # Continue until we've run out of either time or events
    #            while self.current_time < run_to and self.event_queue:
    #
    #                if _DEBUG:
    #                    print('Current Time = ', self.current_time)
    #
    #                # Is there an event scheduled to occur within this run?
    #                if self.event_queue[0].time <= run_to:
    #
    #                    # If so, pick the next transition event from the event queue
    #                    ev = heappop(self.event_queue)
    #
    #                    if _DEBUG:
    #                        print('Event:', ev.time, ev.link, ev.xn_to)
    #
    #                    self.do_transition(ev, self.current_time,
    #                                       plot_each_transition, plotter)
    #
    #                    # Update current time
    #                    self.current_time = ev.time
    #
    #                # If there is no event scheduled for this span of time, simply
    #                # advance current_time to the end of the current run period.
    #                else:
    #                    self.current_time = run_to
    #
    #                if _DEBUG:
    #                    print(self.node_state)

    def run_new(self, run_to, plot_each_transition=False, plotter=None):
        """Test of new approach using priority queue."""

        # Continue until we've run out of either time or events
        while self.current_time < run_to and self.priority_queue._queue:

            if _DEBUG:
                print("Current Time = ", self.current_time)

            # Is there an event scheduled to occur within this run?
            if self.priority_queue._queue[0][0] <= run_to:

                # If so, pick the next transition event from the event queue
                (ev_time, ev_idx, ev_link) = self.priority_queue.pop()

                if _DEBUG:
                    print(
                        "Event:",
                        ev_time,
                        ev_link,
                        self.trn_to[self.next_trn_id[ev_link]],
                    )

                if _CYTEST:
                    do_transition_new(
                        ev_link,
                        ev_time,
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
                        plot_each_transition=False,
                        plotter=None,
                    )
                else:
                    self.do_transition_new(
                        ev_link,
                        ev_time,
                        self.current_time,
                        plot_each_transition,
                        plotter,
                    )

                # Update current time
                self.current_time = ev_time

            # If there is no event scheduled for this span of time, simply
            # advance current_time to the end of the current run period.
            else:
                self.current_time = run_to

            if _DEBUG:
                print(self.node_state)

    def transition_info_as_string(self, event):
        """Returns info about a particular event as a string, for debug."""
        link = event[2]
        tail = self.grid.node_at_link_tail[link]
        head = self.grid.node_at_link_head[link]
        new_link_state = self.trn_to[self.next_trn_id[link]]
        new_tail_state = (new_link_state / self.num_node_states) % self.num_node_states
        new_head_state = new_link_state % self.num_node_states
        info_str = (
            str(event[0])
            + " "  # sched time
            + str(self.next_update[link])
            + " "  # sched time
            + str(link)
            + " "  # link ID
            + str(tail)
            + "=>"  # tail ID
            + str(head)
            + " "  # head ID
            + str(self.link_orientation[link])
            + " "  # orientation
            + str(self.node_state[tail])
            + "=>"  # tail state
            + str(self.node_state[head])
            + " "  # head state
            + str(self.link_state[link])
            + " "  # link state
            + str(self.next_trn_id[link])
            + " "  # trn ID
            + str(new_link_state)
            + " "  # new link state
            + str(new_tail_state)
            + "=>"  # new tail state
            + str(new_head_state)  # new head state
        )
        return info_str

    def print_scheduled_transitions(self):
        """Display list of transitions in PQ, and related data, for debug."""
        print("tme tml lnk tln hdn orn tst hst tid nls nts nhs")
        for trn_event in self.priority_queue._queue:
            print(self.transition_info_as_string(trn_event))
