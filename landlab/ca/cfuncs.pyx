"""
Created on Thu Jun 30 12:40:39 2016

@author: gtucker
"""

import numpy as np
cimport numpy as np
cimport cython
from landlab import CORE_NODE as _CORE
from heapq import heappush

_NEVER = 1.0e50

DTYPE = np.double
ctypedef np.double_t DTYPE_t

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_INT8 = np.int8
ctypedef np.int8_t DTYPE_INT8_t


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


def update_node_states_cython(np.ndarray[DTYPE_INT_t, ndim=1] node_state,
                              np.ndarray[DTYPE_INT8_t, ndim=1] status_at_node,
                              DTYPE_INT_t tail_node, 
                              DTYPE_INT_t head_node,
                              DTYPE_INT_t new_link_state,
                              node_pair):

    # Remember the previous state of each node so we can detect whether the
    # state has changed
    old_tail_node_state = node_state[tail_node]
    old_head_node_state = node_state[head_node]

    # Change to the new states
    if status_at_node[tail_node] == _CORE:
        node_state[tail_node] = node_pair[new_link_state][0]
    if status_at_node[head_node] == _CORE:
        node_state[head_node] = node_pair[new_link_state][1]

    return node_state[tail_node] != old_tail_node_state, \
           node_state[head_node] != old_head_node_state

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
    cdef bint propswap
    cdef double next_time, this_next

    assert (n_xn[current_state] > 0), \
        'must have at least one potential transition'

    # Find next event time for each potential transition
    if n_xn[current_state] == 1:
        my_xn_to = xn_to[current_state][0]
        propswap = xn_propswap[current_state][0]
        next_time = np.random.exponential(
            1.0 / xn_rate[current_state][0])
        prop_update_fn = xn_prop_update_fn[current_state][0]
    else:
        next_time = _NEVER
        my_xn_to = 0
        propswap = False
        for i in range(n_xn[current_state]):
            this_next = np.random.exponential(
                1.0 / xn_rate[current_state][i])
            if this_next < next_time:
                next_time = this_next
                my_xn_to = xn_to[current_state][i]
                propswap = xn_propswap[current_state][i]
                prop_update_fn = xn_prop_update_fn[current_state][i]

    # Create and setup event, and return it
    my_event = Event(next_time + current_time, link,
                     my_xn_to, propswap, prop_update_fn)

    return my_event

def update_link_state(DTYPE_INT_t link, DTYPE_INT_t new_link_state, 
                      DTYPE_t current_time,
                      bnd_lnk,
                      np.ndarray[DTYPE_INT_t, ndim=1] node_state, 
                      np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
                      np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
                      np.ndarray[DTYPE_INT_t, ndim=1] link_orientation,
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

