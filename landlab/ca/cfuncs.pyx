"""
Created on Thu Jun 30 12:40:39 2016

@author: gtucker
"""

import numpy as np

cimport cython
cimport numpy as np

import sys  # for debug

from _heapq import heappop
from _heapq import heappush

from landlab.grid.nodestatus import NodeStatus


cdef double _NEVER = 1.0e50

cdef int _CORE = NodeStatus.CORE

DTYPE = np.double
ctypedef np.double_t DTYPE_t

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t

DTYPE_INT8 = np.int8
ctypedef np.int8_t DTYPE_INT8_t

DTYPE_UINT8 = np.uint8
ctypedef np.uint8_t DTYPE_UINT8_t

cdef char _DEBUG = 0


cdef class PriorityQueue:
    """
    Implements a priority queue.
    """
    cdef public object _queue
    cdef public int _index

    def __init__(self):
        self._queue = []
        self._index = 0

    def push(self, int item, double priority):
        heappush(self._queue, (priority, self._index, item))
        self._index += 1

    def pop(self):
        assert len(self._queue) > 0, "Q is empty"
        return heappop(self._queue)


cdef class Event:
    """
    Represents a transition event at a link. The transition occurs at a given
    link and a given time, and it involves a transition into the state xn_to
    (an integer code representing the new link state; "xn" is shorthand for
    "transition").

    The class overrides the __lt__ (less than operator) method so that when
    Event() objects are placed in a PriorityQueue, the earliest event is
    given the highest priority (i.e., placed at the top of the queue).

    Event() constructor sets 3 required properties and one optional
    property.

    Parameters
    ----------
    time : float
        Time at which the event is scheduled to occur
    link : int
        ID of the link at which event occurs
    xn_to : int
        New state to which this cell pair (link) will transition
    propswap : bool (optional)
        Flag: does this event involve an exchange of properties between
        the two cells?

    Examples
    --------
    >>> from landlab.ca.celllab_cts import Event
    >>> e1 = Event( 10.0, 1, 2)
    >>> e2 = Event( 2.0, 3, 1)
    >>> e1 < e2
    False
    >>> e2 < e1
    True
    """
    cdef public double time
    cdef public int link
    cdef public int xn_to
    cdef public char propswap
    cdef public object prop_update_fn

    def __init__(self, double time, int link, int xn_to, object propswap=False,
                 object prop_update_fn=None):
        """
        Event() constructor sets 3 required properties and one optional
        property.

        Parameters
        ----------
        time : float
            Time at which the event is scheduled to occur
        link : int
            ID of the link at which event occurs
        xn_to : int
            New state to which this cell pair (link) will transition
        propswap : bool (optional)
            Flag: does this event involve an exchange of properties between
            the two cells?
        """
        self.time = time
        self.link = link
        self.xn_to = xn_to
        self.propswap = propswap
        self.prop_update_fn = prop_update_fn

    def __richcmp__(self, Event other, int op):
        """
        Overridden less-than operator: returns true if the event on the left
        has an earlier scheduled time than the event on the right
        """
        return self.time < other.time


@cython.boundscheck(True)
@cython.wraparound(False)
cdef int current_link_state(
    DTYPE_INT_t link_id,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
):
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
    node_state : array of int
        State codes of nodes
    node_at_link_tail : array of int
        ID of node at link tails
    node_at_link_head : array of int
        ID of node at link heads
    link_orientation : array of 1-byte int
        Orientation codes: 0, 1, or (with hex) 2
    num_nodes_states : int
        Total number of possible node states
    num_node_states_sq : int
        Square of number of node states (precomputed for speed)

    Returns
    -------
    int
        New link state code
    """
    cdef int tail_node_state, head_node_state
    cdef char orientation

    # Find out the states of the two nodes, and the orientation
    tail_node_state = node_state[node_at_link_tail[link_id]]
    head_node_state = node_state[node_at_link_head[link_id]]
    orientation = link_orientation[link_id]

    # Return the corresponding state code.
    return (orientation * num_node_states_sq +
            tail_node_state * num_node_states + head_node_state)


@cython.boundscheck(True)
@cython.wraparound(False)
cpdef update_link_states_and_transitions(
    np.ndarray[DTYPE_INT_t, ndim=1] active_links,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    event_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT8_t, ndim=2] xn_propswap,
    xn_prop_update_fn,
):
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
    cdef int current_state
    cdef int i, j

    for j in range(len(active_links)):
        i = active_links[j]
        current_state = current_link_state(
            i,
            node_state,
            node_at_link_tail,
            node_at_link_head,
            link_orientation,
            num_node_states,
            num_node_states_sq,
        )
        if current_state != link_state[i]:
            update_link_state(
                i,
                current_state,
                current_time,
                bnd_lnk,
                node_state,
                node_at_link_tail,
                node_at_link_head,
                link_orientation,
                num_node_states,
                num_node_states_sq,
                link_state,
                n_xn,
                event_queue,
                next_update,
                xn_to,
                xn_rate,
                xn_propswap,
                xn_prop_update_fn,
            )


@cython.boundscheck(True)
@cython.wraparound(False)
cpdef update_link_states_and_transitions_new(
    np.ndarray[DTYPE_INT_t, ndim=1] active_links,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_trn,
    priority_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] next_trn_id,
    np.ndarray[DTYPE_INT_t, ndim=2] trn_id,
    np.ndarray[DTYPE_t, ndim=1] trn_rate,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    DTYPE_t current_time,
):
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
    cdef int current_state
    cdef int i, j

    for j in range(len(active_links)):
        i = active_links[j]
        current_state = current_link_state(
            i,
            node_state,
            node_at_link_tail,
            node_at_link_head,
            link_orientation,
            num_node_states,
            num_node_states_sq,
        )
        if current_state != link_state[i]:
            update_link_state_new(
                i,
                current_state,
                current_time,
                bnd_lnk,
                node_state,
                node_at_link_tail,
                node_at_link_head,
                link_orientation,
                num_node_states,
                num_node_states_sq,
                link_state,
                n_trn,
                priority_queue,
                next_update,
                next_trn_id,
                trn_id,
                trn_rate,
            )


@cython.boundscheck(True)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef update_node_states(
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    DTYPE_INT_t tail_node,
    DTYPE_INT_t head_node,
    DTYPE_INT_t new_link_state,
    DTYPE_INT_t num_states,
):
    """Update the states of 2 nodes that underwent a transition."""
    if _DEBUG:
        print(("UNS", tail_node, head_node, new_link_state, num_states))
    # Change to the new states
    if status_at_node[tail_node] == _CORE:
        # assume integer division!!
        node_state[tail_node] = (new_link_state / num_states) % num_states
    if status_at_node[head_node] == _CORE:
        node_state[head_node] = new_link_state % num_states
    if _DEBUG:
        print(("UNS new tail state: ", node_state[tail_node]))
        print(("UNS new head state: ", node_state[head_node]))


@cython.boundscheck(True)
@cython.wraparound(False)
cpdef get_next_event(
    DTYPE_INT_t link, DTYPE_INT_t current_state,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    np.ndarray[DTYPE_INT8_t, ndim=2] xn_propswap,
    xn_prop_update_fn,
):
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
    (see celllab_cts.py for other parameters)

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
    cdef int my_xn_to
    cdef int i
    cdef char propswap
    cdef double next_time, this_next
    cdef Event my_event

    # Find next event time for each potential transition
    if n_xn[current_state] == 1:
        my_xn_to = xn_to[current_state, 0]
        propswap = xn_propswap[current_state, 0]
        next_time = np.random.exponential(1.0 / xn_rate[current_state, 0])
        # next_time = -(1.0 / xn_rate[current_state, 0]) * log(1.0 - rand())
        prop_update_fn = xn_prop_update_fn[current_state, 0]
    else:
        next_time = _NEVER
        my_xn_to = 0
        propswap = 0
        for i in range(n_xn[current_state]):
            this_next = np.random.exponential(1.0 / xn_rate[current_state, i])
            # this_next = -(1.0 / xn_rate[current_state, i]) * log(1.0 - rand())
            if this_next < next_time:
                next_time = this_next
                my_xn_to = xn_to[current_state, i]
                propswap = xn_propswap[current_state, i]
                prop_update_fn = xn_prop_update_fn[current_state, i]

    # Create and setup event, and return it
    my_event = Event(next_time + current_time, link, my_xn_to, propswap, prop_update_fn)

    return my_event


@cython.boundscheck(True)
@cython.wraparound(False)
cpdef get_next_event_new(
    DTYPE_INT_t link, DTYPE_INT_t current_state,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT_t, ndim=1] n_trn,
    np.ndarray[DTYPE_INT_t, ndim=2] trn_id,
    np.ndarray[DTYPE_t, ndim=1] trn_rate,
):
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
    (see celllab_cts.py for other parameters)

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
    cdef int this_trn_id
    cdef int i
    cdef double next_time, this_next

    # Find next event time for each potential transition
    if n_trn[current_state] == 1:
        this_trn_id = trn_id[current_state, 0]
        next_time = np.random.exponential(1.0 / trn_rate[this_trn_id])
    else:
        next_time = _NEVER
        this_trn_id = -1
        for i in range(n_trn[current_state]):
            this_next = np.random.exponential(1.0 / trn_rate[trn_id[current_state][i]])
            if this_next < next_time:
                next_time = this_next
                this_trn_id = trn_id[current_state, i]

    return (next_time + current_time, this_trn_id)


cpdef push_transitions_to_event_queue(
    int number_of_active_links,
    np.ndarray[DTYPE_INT_t, ndim=1] active_links,
    np.ndarray[DTYPE_INT_t, ndim=1] n_trn,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=2] trn_id,
    np.ndarray[DTYPE_t, ndim=1] trn_rate,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] next_trn_id,
    PriorityQueue priority_queue,
):
    """
    Initializes the event queue by creating transition events for each
    cell pair that has one or more potential transitions and pushing these
    onto the queue. Also records scheduled transition times in the
    self.next_update array.
    """
    for j in range(number_of_active_links):

        i = active_links[j]
        if n_trn[link_state[i]] > 0:
            (ev_time, this_trn_id) = get_next_event_new(
                i, link_state[i], 0.0, n_trn, trn_id, trn_rate
            )
            priority_queue.push(i, ev_time)
            next_update[i] = ev_time
            next_trn_id[i] = this_trn_id

        else:
            next_update[i] = _NEVER


@cython.boundscheck(True)
@cython.wraparound(False)
cdef void update_link_state(
    DTYPE_INT_t link, DTYPE_INT_t new_link_state,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn, event_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    np.ndarray[DTYPE_INT8_t, ndim=2] xn_propswap,
    np.ndarray[object, ndim=2] xn_prop_update_fn,
):
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
    (see celllab_cts.py for other parameters)
    """
    cdef int fns, tns
    cdef int orientation
    cdef Event event

    # If the link connects to a boundary, we might have a different state
    # than the one we planned
    if bnd_lnk[link]:
        fns = node_state[node_at_link_tail[link]]
        tns = node_state[node_at_link_head[link]]
        orientation = link_orientation[link]
        new_link_state = orientation * num_node_states_sq + \
            fns * num_node_states + tns

    link_state[link] = new_link_state
    if n_xn[new_link_state] > 0:
        event = get_next_event(
            link,
            new_link_state,
            current_time,
            n_xn,
            xn_to,
            xn_rate,
            xn_propswap,
            xn_prop_update_fn,
        )
        heappush(event_queue, event)
        next_update[link] = event.time
    else:
        next_update[link] = _NEVER


@cython.boundscheck(True)
@cython.wraparound(False)
cdef void update_link_state_new(
    DTYPE_INT_t link, DTYPE_INT_t new_link_state,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_trn,
    PriorityQueue priority_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] next_trn_id,
    np.ndarray[DTYPE_INT_t, ndim=2] trn_id,
    np.ndarray[DTYPE_t, ndim=1] trn_rate,
):
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
    (see celllab_cts.py for other parameters)
    """
    cdef int fns, tns
    cdef int this_trn_id
    cdef int orientation

    if _DEBUG:
        print(("ULSN", link, link_state[link], new_link_state, current_time))

    # If the link connects to a boundary, we might have a different state
    # than the one we planned
    if bnd_lnk[link]:
        fns = node_state[node_at_link_tail[link]]
        tns = node_state[node_at_link_head[link]]
        orientation = link_orientation[link]
        new_link_state = orientation * num_node_states_sq + \
            fns * num_node_states + tns
        if _DEBUG:
            print((" bnd True", new_link_state))

    link_state[link] = new_link_state
    if n_trn[new_link_state] > 0:
        (event_time, this_trn_id) = get_next_event_new(
            link, new_link_state, current_time, n_trn, trn_id, trn_rate
        )
        priority_queue.push(link, event_time)
        next_update[link] = event_time
        next_trn_id[link] = this_trn_id
    else:
        next_update[link] = _NEVER
        next_trn_id[link] = -1


@cython.boundscheck(True)
@cython.wraparound(False)
cdef void do_transition(
    Event event,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT_t, ndim=1] propid,
    object prop_data,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    DTYPE_INT_t prop_reset_value,
    np.ndarray[DTYPE_INT8_t, ndim=2] xn_propswap,
    xn_prop_update_fn,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    event_queue,
    this_cts_model,
    plot_each_transition=False,
    plotter=None,
):
    """Transition state.

    Implements a state transition.

    Parameters
    ----------
    event : Event object
        Event object containing the data for the current transition event
    plot_each_transition : bool (optional)
        True if caller wants to show a plot of the grid after this
        transition
    plotter : CAPlotter object
        Sent if caller wants a plot after this transition
    (see celllab_cts.py for other parameters)

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
    cdef int tail_node, head_node  # IDs of tail and head nodes at link
    cdef int old_tail_node_state
    cdef int old_head_node_state
    cdef int this_link_tail_node   # Tail ID for an adjacent link
    cdef int this_link_head_node   # Head ID for an adjacent link
    cdef int link                  # ID of a link
    cdef int new_link_state        # New link state after transition
    cdef int tmp                   # Used to exchange property IDs
    cdef char dir_code             # Direction code for link at node
    cdef char orientation          # Orientation code for link
    cdef int i

    # We'll process the event if its update time matches the one we have
    # recorded for the link in question. If not, it means that the link has
    # changed state since the event was pushed onto the event queue, and
    # in that case we'll ignore it.
    if event.time == next_update[event.link]:

        tail_node = node_at_link_tail[event.link]
        head_node = node_at_link_head[event.link]

        # Remember the previous state of each node so we can detect whether the
        # state has changed
        old_tail_node_state = node_state[tail_node]
        old_head_node_state = node_state[head_node]

        update_node_states(
            node_state,
            status_at_node,
            tail_node,
            head_node,
            event.xn_to,
            num_node_states,
        )
        update_link_state(
            event.link,
            event.xn_to,
            event.time,
            bnd_lnk,
            node_state,
            node_at_link_tail,
            node_at_link_head,
            link_orientation,
            num_node_states,
            num_node_states_sq,
            link_state,
            n_xn,
            event_queue,
            next_update,
            xn_to,
            xn_rate,
            xn_propswap,
            xn_prop_update_fn,
        )

        # Next, when the state of one of the link's nodes changes, we have
        # to update the states of the OTHER links attached to it. This
        # could happen to one or both nodes.
        if node_state[tail_node] != old_tail_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[tail_node, i]
                dir_code = active_link_dirs_at_node[tail_node, i]

                if dir_code != 0 and link != event.link:

                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node]
                    )
                    update_link_state(
                        link,
                        new_link_state,
                        event.time,
                        bnd_lnk,
                        node_state,
                        node_at_link_tail,
                        node_at_link_head,
                        link_orientation,
                        num_node_states,
                        num_node_states_sq,
                        link_state,
                        n_xn,
                        event_queue,
                        next_update,
                        xn_to,
                        xn_rate,
                        xn_propswap,
                        xn_prop_update_fn,
                    )

        if node_state[head_node] != old_head_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[head_node, i]
                dir_code = active_link_dirs_at_node[head_node, i]

                if dir_code != 0 and link != event.link:
                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node]
                    )
                    update_link_state(
                        link,
                        new_link_state,
                        event.time,
                        bnd_lnk,
                        node_state,
                        node_at_link_tail,
                        node_at_link_head,
                        link_orientation,
                        num_node_states,
                        num_node_states_sq,
                        link_state,
                        n_xn,
                        event_queue,
                        next_update,
                        xn_to,
                        xn_rate,
                        xn_propswap,
                        xn_prop_update_fn,
                    )

        # If requested, display a plot of the grid
        if plot_each_transition and (plotter is not None):
            plotter.update_plot()

        # If this event involves an exchange of properties (i.e., the
        # event involves motion of an object that posses properties we
        # want to track), implement the swap.
        #   If the event requires a call to a user-defined callback
        # function, we handle that here too.
        if event.propswap:
            tmp = propid[tail_node]
            propid[tail_node] = propid[head_node]
            propid[head_node] = tmp
            if status_at_node[tail_node] != _CORE:
                prop_data[propid[tail_node]] = prop_reset_value
            if status_at_node[head_node] != _CORE:
                prop_data[propid[head_node]] = prop_reset_value
            if event.prop_update_fn is not None:
                event.prop_update_fn(this_cts_model, tail_node, head_node, event.time)


# @cython.boundscheck(False)
# @cython.wraparound(False)
cpdef void do_transition_new(
    DTYPE_INT_t event_link,
    DTYPE_t event_time,
    PriorityQueue priority_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] next_trn_id,
    np.ndarray[DTYPE_INT_t, ndim=1] trn_to,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_trn,
    np.ndarray[DTYPE_INT_t, ndim=2] trn_id,
    np.ndarray[DTYPE_t, ndim=1] trn_rate,
    np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=1] trn_propswap,
    np.ndarray[DTYPE_INT_t, ndim=1] propid,
    object prop_data,
    DTYPE_INT_t prop_reset_value,
    object trn_prop_update_fn,
    object this_cts_model,
    plot_each_transition=False,
    plotter=None,
):
    """Transition state.

    Implements a state transition.

    Parameters
    ----------
    event : Event object
        Event object containing the data for the current transition event
    plot_each_transition : bool (optional)
        True if caller wants to show a plot of the grid after this
        transition
    plotter : CAPlotter object
        Sent if caller wants a plot after this transition
    (see celllab_cts.py for other parameters)

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
    cdef int tail_node, head_node  # IDs of tail and head nodes at link
    cdef int old_tail_node_state
    cdef int old_head_node_state
    cdef int this_trn_id
    cdef int this_trn_to
    cdef int this_link_tail_node   # Tail ID for an adjacent link
    cdef int this_link_head_node   # Head ID for an adjacent link
    cdef int link                  # ID of a link
    cdef int new_link_state        # New link state after transition
    cdef int tmp                   # Used to exchange property IDs
    cdef char dir_code             # Direction code for link at node
    cdef char orientation          # Orientation code for link
    cdef int i

    if _DEBUG:
        print(
            (
                "DTN",
                event_time,
                event_link,
                link_state[event_link],
                next_update[event_link],
            )
        )

    # We'll process the event if its update time matches the one we have
    # recorded for the link in question. If not, it means that the link has
    # changed state since the event was pushed onto the event queue, and
    # in that case we'll ignore it.
    if event_time == next_update[event_link]:

        tail_node = node_at_link_tail[event_link]
        head_node = node_at_link_head[event_link]

        # DEBUG
        if status_at_node[tail_node] == 4 or status_at_node[head_node] == 4:
            print(
                (
                    "TRN INFO: ",
                    event_time,
                    event_link,
                    link_state[event_link],
                    next_update[event_link],
                )
            )
            print("TAIL " + str(tail_node) + " " + status_at_node[tail_node])
            print("HEAD " + str(head_node) + " " + status_at_node[tail_node])
            # _DEBUG = True

        # Remember the previous state of each node so we can detect whether the
        # state has changed
        old_tail_node_state = node_state[tail_node]
        old_head_node_state = node_state[head_node]

        this_trn_id = next_trn_id[event_link]
        this_trn_to = trn_to[this_trn_id]

        if _DEBUG:
            print((this_trn_id, this_trn_to))
            print(("tail:", tail_node))
            print(("tail state:", old_tail_node_state))
            print(("head:", head_node))
            print(("head state:", old_head_node_state))

        update_node_states(
            node_state,
            status_at_node,
            tail_node,
            head_node,
            this_trn_to,
            num_node_states,
        )
        update_link_state_new(
            event_link,
            this_trn_to,
            event_time,
            bnd_lnk,
            node_state,
            node_at_link_tail,
            node_at_link_head,
            link_orientation,
            num_node_states,
            num_node_states_sq,
            link_state,
            n_trn,
            priority_queue,
            next_update,
            next_trn_id,
            trn_id,
            trn_rate,
        )

        # Next, when the state of one of the link's nodes changes, we have
        # to update the states of the OTHER links attached to it. This
        # could happen to one or both nodes.
        if node_state[tail_node] != old_tail_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[tail_node, i]
                dir_code = active_link_dirs_at_node[tail_node, i]

                if dir_code != 0 and link != event_link:

                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node]
                    )
                    update_link_state_new(
                        link,
                        new_link_state,
                        event_time,
                        bnd_lnk,
                        node_state,
                        node_at_link_tail,
                        node_at_link_head,
                        link_orientation,
                        num_node_states,
                        num_node_states_sq,
                        link_state,
                        n_trn,
                        priority_queue,
                        next_update,
                        next_trn_id,
                        trn_id,
                        trn_rate,
                    )

        if node_state[head_node] != old_head_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[head_node, i]
                dir_code = active_link_dirs_at_node[head_node, i]

                if dir_code != 0 and link != event_link:
                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node]
                    )
                    update_link_state_new(
                        link,
                        new_link_state,
                        event_time,
                        bnd_lnk,
                        node_state,
                        node_at_link_tail,
                        node_at_link_head,
                        link_orientation,
                        num_node_states,
                        num_node_states_sq,
                        link_state,
                        n_trn,
                        priority_queue,
                        next_update,
                        next_trn_id,
                        trn_id,
                        trn_rate,
                    )

        # If requested, display a plot of the grid
        if plot_each_transition and (plotter is not None):
            plotter.update_plot()

        # If this event involves an exchange of properties (i.e., the
        # event involves motion of an object that posses properties we
        # want to track), implement the swap.
        #   If the event requires a call to a user-defined callback
        # function, we handle that here too.
        if trn_propswap[this_trn_id]:
            tmp = propid[tail_node]
            propid[tail_node] = propid[head_node]
            propid[head_node] = tmp
            if status_at_node[tail_node] != _CORE:
                prop_data[propid[tail_node]] = prop_reset_value
            if status_at_node[head_node] != _CORE:
                prop_data[propid[head_node]] = prop_reset_value
            if trn_prop_update_fn[this_trn_id] != 0:
                trn_prop_update_fn[this_trn_id](
                    this_cts_model, tail_node, head_node, event_time)

cpdef double run_cts_new(
    double run_to,
    double current_time,
    PriorityQueue priority_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] next_trn_id,
    np.ndarray[DTYPE_INT_t, ndim=1] trn_to,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_trn,
    np.ndarray[DTYPE_INT_t, ndim=2] trn_id,
    np.ndarray[DTYPE_t, ndim=1] trn_rate,
    np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=1] trn_propswap,
    np.ndarray[DTYPE_INT_t, ndim=1] propid,
    object prop_data,
    DTYPE_INT_t prop_reset_value,
    trn_prop_update_fn,
    this_cts_model,
    char plot_each_transition,
    object plotter,
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
    (see celllab_cts.py for other parameters)
    """
    import sys
    cdef double ev_time
    cdef int _ev_idx
    cdef int ev_link

    # Continue until we've run out of either time or events
    while current_time < run_to and priority_queue._queue:

        if _DEBUG:
            print("current time = ", current_time)

        # Is there an event scheduled to occur within this run?
        if priority_queue._queue[0][0] <= run_to:

            # If so, pick the next transition event from the event queue
            (ev_time, _ev_idx, ev_link) = priority_queue.pop()

            # ... and execute the transition
            do_transition_new(
                ev_link,
                ev_time,
                priority_queue,
                next_update,
                node_at_link_tail,
                node_at_link_head,
                node_state,
                next_trn_id,
                trn_to,
                status_at_node,
                num_node_states,
                num_node_states_sq,
                bnd_lnk,
                link_orientation,
                link_state,
                n_trn,
                trn_id,
                trn_rate,
                links_at_node,
                active_link_dirs_at_node,
                trn_propswap,
                propid, prop_data,
                prop_reset_value,
                trn_prop_update_fn,
                this_cts_model,
                plot_each_transition,
                plotter,
            )

            # Update current time
            current_time = ev_time

        # If there is no event scheduled for this span of time, simply
        # advance current_time to the end of the current run period.
        else:
            current_time = run_to

    return current_time


cpdef double run_cts(
    double run_to,
    double current_time,
    char plot_each_transition,
    object plotter,
    object event_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT_t, ndim=1] propid,
    object prop_data,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    DTYPE_INT_t prop_reset_value,
    np.ndarray[DTYPE_INT8_t, ndim=2] xn_propswap,
    xn_prop_update_fn,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    this_cts_model,
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
    (see celllab_cts.py for other parameters)
    """
    cdef Event ev

    # Continue until we've run out of either time or events
    while current_time < run_to and event_queue:

        # Is there an event scheduled to occur within this run?
        if event_queue[0].time <= run_to:

            # If so, pick the next transition event from the event queue
            ev = heappop(event_queue)

            # ... and execute the transition
            do_transition(
                ev,
                next_update,
                node_at_link_tail,
                node_at_link_head,
                node_state,
                link_state,
                status_at_node,
                link_orientation,
                propid,
                prop_data,
                n_xn,
                xn_to,
                xn_rate,
                links_at_node,
                active_link_dirs_at_node,
                num_node_states,
                num_node_states_sq,
                prop_reset_value,
                xn_propswap,
                xn_prop_update_fn,
                bnd_lnk,
                event_queue,
                this_cts_model,
                plot_each_transition,
                plotter,
            )

            # Update current time
            current_time = ev.time

        # If there is no event scheduled for this span of time, simply
        # advance current_time to the end of the current run period.
        else:
            current_time = run_to

    return current_time


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_next_event_lean(
    DTYPE_INT_t link,
    DTYPE_INT_t current_state,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
):
    """Get the next event for a link.

    Returns the next event for link with ID "link", which is in state
    "current state". This "lean" version omits parameters related to property
    exchange and callback function.

    Parameters
    ----------
    link : int
        ID of the link
    current_state : int
        Current state code for the link
    current_time : float
        Current time in simulation (i.e., time of event just processed)
    (see celllab_cts.py for other parameters)

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
    cdef int my_xn_to
    cdef int i
    cdef double next_time, this_next
    cdef Event my_event

    # Find next event time for each potential transition
    if n_xn[current_state] == 1:
        my_xn_to = xn_to[current_state, 0]
        next_time = np.random.exponential(1.0 / xn_rate[current_state, 0])
        # next_time = -(1.0 / xn_rate[current_state, 0]) * log(1.0 - rand())
    else:
        next_time = _NEVER
        my_xn_to = 0
        for i in range(n_xn[current_state]):
            this_next = np.random.exponential(1.0 / xn_rate[current_state, i])
            # this_next = -(1.0 / xn_rate[current_state, i]) * log(1.0 - rand())
            if this_next < next_time:
                next_time = this_next
                my_xn_to = xn_to[current_state, i]

    # Create and setup event, and return it
    my_event = Event(next_time + current_time, link, my_xn_to)

    return my_event


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void update_link_state_lean(
    DTYPE_INT_t link,
    DTYPE_INT_t new_link_state,
    DTYPE_t current_time,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn, event_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
):
    """
    Implements a link transition by updating the current state of the link
    and (if appropriate) choosing the next transition event and pushing it
    on to the event queue. This "lean" version omits parameters related to
    property exchange and callback function.

    Parameters
    ----------
    link : int
        ID of the link to update
    new_link_state : int
        Code for the new state
    current_time : float
        Current time in simulation
    (see celllab_cts.py for other parameters)
    """
    cdef int fns, tns
    cdef int orientation
    cdef Event event

    # If the link connects to a boundary, we might have a different state
    # than the one we planned
    if bnd_lnk[link]:
        fns = node_state[node_at_link_tail[link]]
        tns = node_state[node_at_link_head[link]]
        orientation = link_orientation[link]
        new_link_state = orientation * num_node_states_sq + fns * num_node_states + tns

    link_state[link] = new_link_state
    if n_xn[new_link_state] > 0:
        event = get_next_event_lean(
            link, new_link_state, current_time, n_xn, xn_to, xn_rate
        )
        heappush(event_queue, event)
        next_update[link] = event.time
    else:
        next_update[link] = _NEVER


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void do_transition_lean(
    Event event,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
    object event_queue,
):
    """Transition state.

    Implements a state transition. This "lean" version omits parameters related
    to property exchange and callback function.

    Parameters
    ----------
    event : Event object
        Event object containing the data for the current transition event
    (see celllab_cts.py for other parameters)

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
    cdef int tail_node, head_node  # IDs of tail and head nodes at link
    cdef int old_tail_node_state
    cdef int old_head_node_state
    cdef int this_link_tail_node   # Tail ID for an adjacent link
    cdef int this_link_head_node   # Head ID for an adjacent link
    cdef int link                  # ID of a link
    cdef int new_link_state        # New link state after transition
    cdef char dir_code             # Direction code for link at node
    cdef char orientation          # Orientation code for link
    cdef int i

    # print 'dtl'

    # We'll process the event if its update time matches the one we have
    # recorded for the link in question. If not, it means that the link has
    # changed state since the event was pushed onto the event queue, and
    # in that case we'll ignore it.
    if event.time == next_update[event.link]:

        tail_node = node_at_link_tail[event.link]
        head_node = node_at_link_head[event.link]

        # Remember the previous state of each node so we can detect whether the
        # state has changed
        old_tail_node_state = node_state[tail_node]
        old_head_node_state = node_state[head_node]

        update_node_states(
            node_state,
            status_at_node,
            tail_node,
            head_node,
            event.xn_to,
            num_node_states,
        )
        update_link_state_lean(
            event.link,
            event.xn_to,
            event.time,
            bnd_lnk,
            node_state,
            node_at_link_tail,
            node_at_link_head,
            link_orientation,
            num_node_states,
            num_node_states_sq,
            link_state,
            n_xn,
            event_queue,
            next_update,
            xn_to,
            xn_rate,
        )

        # Next, when the state of one of the link's nodes changes, we have
        # to update the states of the OTHER links attached to it. This
        # could happen to one or both nodes.
        if node_state[tail_node] != old_tail_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[tail_node, i]
                dir_code = active_link_dirs_at_node[tail_node, i]

                if dir_code != 0 and link != event.link:

                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node]
                    )
                    update_link_state_lean(
                        link,
                        new_link_state,
                        event.time,
                        bnd_lnk,
                        node_state,
                        node_at_link_tail,
                        node_at_link_head,
                        link_orientation,
                        num_node_states,
                        num_node_states_sq,
                        link_state,
                        n_xn,
                        event_queue,
                        next_update,
                        xn_to,
                        xn_rate,
                    )

        if node_state[head_node] != old_head_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[head_node, i]
                dir_code = active_link_dirs_at_node[head_node, i]

                if dir_code != 0 and link != event.link:
                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node]
                    )
                    update_link_state_lean(
                        link,
                        new_link_state,
                        event.time,
                        bnd_lnk,
                        node_state,
                        node_at_link_tail,
                        node_at_link_head,
                        link_orientation,
                        num_node_states,
                        num_node_states_sq,
                        link_state,
                        n_xn,
                        event_queue,
                        next_update,
                        xn_to,
                        xn_rate,
                    )


@cython.boundscheck(False)
cpdef double run_cts_lean(
    double run_to,
    double current_time,
    object event_queue,
    np.ndarray[DTYPE_t, ndim=1] next_update,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
    np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
    np.ndarray[DTYPE_INT_t, ndim=1] node_state,
    np.ndarray[DTYPE_INT_t, ndim=1] link_state,
    np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
    np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
    np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
    np.ndarray[DTYPE_t, ndim=2] xn_rate,
    np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
    np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
    DTYPE_INT_t num_node_states,
    DTYPE_INT_t num_node_states_sq,
    np.ndarray[DTYPE_INT8_t, ndim=1] bnd_lnk,
):
    """Run the model forward for a specified period of time. This "lean"
    version omits parameters related to property exchange and callback fn.

    Parameters
    ----------
    run_to : float
        Time to run to, starting from self.current_time
    (see celllab_cts.py for other parameters)
    """
    cdef Event ev

    # Continue until we've run out of either time or events
    while current_time < run_to and event_queue:

        # Is there an event scheduled to occur within this run?
        if event_queue[0].time <= run_to:

            # print 'popping'

            # If so, pick the next transition event from the event queue
            ev = heappop(event_queue)

            # ... and execute the transition
            do_transition_lean(
                ev,
                next_update,
                node_at_link_tail,
                node_at_link_head,
                node_state,
                link_state,
                status_at_node,
                link_orientation,
                n_xn,
                xn_to,
                xn_rate,
                links_at_node,
                active_link_dirs_at_node,
                num_node_states,
                num_node_states_sq,
                bnd_lnk,
                event_queue,
            )

            # Update current time
            current_time = ev.time

        # If there is no event scheduled for this span of time, simply
        # advance current_time to the end of the current run period.
        else:
            current_time = run_to

    return current_time
