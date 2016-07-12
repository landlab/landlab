"""
Created on Thu Jun 30 12:40:39 2016

@author: gtucker
"""

import numpy as np
cimport numpy as np
cimport cython
from landlab import CORE_NODE
from _heapq import heappush

_NEVER = 1.0e50

cdef int _CORE = CORE_NODE

DTYPE = np.double
ctypedef np.double_t DTYPE_t

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_INT8 = np.int8
ctypedef np.int8_t DTYPE_INT8_t

_DEBUG = False



class Event():
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

    def __init__(self, time, link, xn_to, propswap=False, prop_update_fn=None):
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

    def __lt__(self, other):
        """
        Overridden less-than operator: returns true if the event on the left
        has an earlier scheduled time than the event on the right
        """
        return self.time < other.time


@cython.boundscheck(False)
def update_node_states(np.ndarray[DTYPE_INT_t, ndim=1] node_state,
                       np.ndarray[DTYPE_INT8_t, ndim=1] status_at_node,
                       DTYPE_INT_t tail_node, 
                       DTYPE_INT_t head_node,
                       DTYPE_INT_t new_link_state,
                       DTYPE_INT_t num_states):

    # Change to the new states
    if status_at_node[tail_node] == _CORE:
        node_state[tail_node] = (new_link_state // num_states) % num_states
    if status_at_node[head_node] == _CORE:
        node_state[head_node] = new_link_state % num_states


@cython.boundscheck(False)
def get_next_event(DTYPE_INT_t link, DTYPE_INT_t current_state, 
                   DTYPE_t current_time, 
                   np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
                   np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
                   np.ndarray[DTYPE_t, ndim=2] xn_rate,
                   xn_propswap, xn_prop_update_fn):
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
    cdef int my_xn_to
    cdef int i
    cdef bint propswap
    cdef double next_time, this_next

    assert (n_xn[current_state] > 0), \
        'must have at least one potential transition'

    # Find next event time for each potential transition
    if n_xn[current_state] == 1:
        my_xn_to = xn_to[current_state][0]
        propswap = xn_propswap[current_state][0]
        next_time = np.random.exponential(1.0 / xn_rate[current_state][0])
        prop_update_fn = xn_prop_update_fn[current_state][0]
    else:
        next_time = _NEVER
        my_xn_to = 0
        propswap = False
        for i in range(n_xn[current_state]):
            this_next = np.random.exponential(1.0 / xn_rate[current_state][i])
            if this_next < next_time:
                next_time = this_next
                my_xn_to = xn_to[current_state][i]
                propswap = xn_propswap[current_state][i]
                prop_update_fn = xn_prop_update_fn[current_state][i]

    # Create and setup event, and return it
    my_event = Event(next_time + current_time, link,
                     my_xn_to, propswap, prop_update_fn)

    if _DEBUG:
        print('get_next_event():')
        print(('  next_time:', my_event.time))
        print(('  link:', my_event.link))
        print(('  xn_to:', my_event.xn_to))
        print(('  propswap:', my_event.propswap))

    return my_event


@cython.boundscheck(False)
def update_link_state(DTYPE_INT_t link, DTYPE_INT_t new_link_state, 
                      DTYPE_t current_time,
                      bnd_lnk,
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
                      xn_propswap, xn_prop_update_fn):
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
    cdef int fns, tns
    cdef int orientation

    if _DEBUG:
        print('update_link_state() link ' + str(link) + ' to state ' + str(new_link_state))
    # If the link connects to a boundary, we might have a different state
    # than the one we planned
    # if self.grid.status_at_node[self.grid.link_fromnode[link]]!=_CORE or \
    #   self.grid.status_at_node[self.grid.link_tonode[link]]!=_CORE:
    if bnd_lnk[link]:
        fns = node_state[node_at_link_tail[link]]
        tns = node_state[node_at_link_head[link]]
        orientation = link_orientation[link]
        new_link_state = orientation * num_node_states_sq + \
            fns * num_node_states + tns

    link_state[link] = new_link_state
    if n_xn[new_link_state] > 0:
        event = get_next_event(link, new_link_state, current_time, n_xn, xn_to,
                               xn_rate, xn_propswap, xn_prop_update_fn)
        heappush(event_queue, event)
        next_update[link] = event.time
    else:
        next_update[link] = _NEVER


@cython.boundscheck(False)
def do_transition(event,
                  np.ndarray[DTYPE_t, ndim=1] next_update,                  
                  np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,                  
                  np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,                  
                  np.ndarray[DTYPE_INT_t, ndim=1] node_state,            
                  np.ndarray[DTYPE_INT_t, ndim=1] link_state,
                  np.ndarray[DTYPE_INT8_t, ndim=1] status_at_node,
                  np.ndarray[DTYPE_INT8_t, ndim=1] link_orientation,
                  np.ndarray[DTYPE_INT_t, ndim=1] propid,
                  prop_data,
                  np.ndarray[DTYPE_INT_t, ndim=1] n_xn,
                  np.ndarray[DTYPE_INT_t, ndim=2] xn_to,
                  np.ndarray[DTYPE_t, ndim=2] xn_rate, 
                  np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
                  np.ndarray[DTYPE_INT8_t, ndim=2] active_link_dirs_at_node,
                  DTYPE_INT_t num_node_states,
                  DTYPE_INT_t num_node_states_sq,
                  DTYPE_INT_t prop_reset_value,
                  xn_propswap,
                  xn_prop_update_fn,
                  bnd_lnk, event_queue,
                  this_cts_model,
                  plot_each_transition=False,
                  plotter=None):
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
    cdef int this_link_tail_node   # Tail ID for an adjacent link
    cdef int this_link_head_node   # Head ID for an adjacent link
    cdef int link                  # ID of a link
    cdef int new_link_state        # New link state after transition
    cdef int tmp                   # Used to exchange property IDs
    cdef char tail_changed, head_changed,  # Booleans
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

        update_node_states(node_state, status_at_node, tail_node,
                           head_node, event.xn_to, num_node_states)
        update_link_state(event.link, event.xn_to, event.time,
                          bnd_lnk, node_state,
                          node_at_link_tail,
                          node_at_link_head,
                          link_orientation, num_node_states,
                          num_node_states_sq, link_state,
                          n_xn, event_queue,
                          next_update, xn_to, xn_rate,
                          xn_propswap, xn_prop_update_fn)

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
                    current_pair = (node_state[this_link_tail_node],
                                    node_state[this_link_head_node],
                                    orientation)
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node])
                    update_link_state(link, new_link_state, event.time,
                                      bnd_lnk, node_state,
                                      node_at_link_tail,
                                      node_at_link_head,
                                      link_orientation, num_node_states,
                                      num_node_states_sq, link_state,
                                      n_xn, event_queue,
                                      next_update, xn_to, xn_rate,
                                      xn_propswap, xn_prop_update_fn)

        if node_state[head_node] != old_head_node_state:

            for i in range(links_at_node.shape[1]):

                link = links_at_node[head_node, i]
                dir_code = active_link_dirs_at_node[head_node, i]

                if dir_code != 0 and link != event.link:
                    this_link_tail_node = node_at_link_tail[link]
                    this_link_head_node = node_at_link_head[link]
                    orientation = link_orientation[link]
                    current_pair = (node_state[this_link_tail_node],
                                    node_state[this_link_head_node],
                                    orientation)
                    new_link_state = (
                        orientation * num_node_states_sq +
                        node_state[this_link_tail_node] * num_node_states +
                        node_state[this_link_head_node])
                    update_link_state(link, new_link_state, event.time,
                                      bnd_lnk, node_state,
                                      node_at_link_tail,
                                      node_at_link_head,
                                      link_orientation, num_node_states,
                                      num_node_states_sq, link_state,
                                      n_xn, event_queue,
                                      next_update, xn_to, xn_rate,
                                      xn_propswap, xn_prop_update_fn)

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
                event.prop_update_fn(
                    this_cts_model, tail_node, head_node, event.time)

