#! /usr/env/python
"""
Link-based cellular automaton modeling tools.

Created GT Oct 2013
"""

from heapq import heappush
from heapq import heappop
from landlab import RasterModelGrid
import landlab
from landlab.components.fracture_grid.fracture_grid import *
import numpy
from landlab.plot import imshow_grid
import pylab as plt

_ROCK = 0
_SAP = 1
_NODE_STATE = { _ROCK : 'rock', \
                _SAP : 'saprolite' }
_LINK_STATE = [ (_ROCK, _ROCK), (_ROCK, _SAP), (_SAP, _ROCK), (_SAP, _SAP) ]
_NEVER = 1e12

_DEBUG = False

_TEST = False


class Transition():
    """
    Represents a transition from one state ("from_state") to another 
    ("to_state") at a link. The transition probability is represented by a rate 
    parameter "rate", with dimensions of 1/T. The probability distribution of 
    time until the transition event occurs is exponentional with mean 1/rate. 
    The optional name parameter allows the caller to assign a name to any given 
    transition.
    """
    def __init__(self, from_state, to_state, rate, name=None):
        
        self.from_state = from_state
        self.to_state = to_state
        self.rate = rate
        self.name = name
        
        
class Event():
    """
    Represents a transition event at a link. The transition occurs at a given
    link and a given time, and it involves a transition into the state xn_to
    (an integer code representing the new link state; "xn" is shorthand for
    "transition").
    
    The class overrides the __lt__ (less than operator) method so that when
    Event() objects are placed in a PriorityQueue, the earliest event is
    given the highest priority (i.e., placed at the top of the queue).
    
    Example:
        
        >>> e1 = Event( 10.0, 1, 2)
        >>> e2 = Event( 2.0, 3, 1)
        >>> e1 < e2
        False
        >>> e2 < e1
        True
    """
    def __init__(self, time, link, xn_to):

        self.time = time
        self.link = link
        self.xn_to = xn_to
        
    def __lt__(self, other):
        
        if self.time < other.time:
            return True
        else:
            return False
        
        
class LinkCellularAutomaton():
    """
    The LinkCellularAutomaton implements a link-type (or doublet-type) cellular
    automaton model. A link connects a pair of cells. Each cell has a state
    (represented by an integer code), and each link also has a state that is
    determined by the states of the cell pair.
    """
    def __init__(self, model_grid, node_state_dict, transition_list,
                 initial_node_states):
                 
        # Keep a copy of the model grid
        assert (type(model_grid) is landlab.grid.raster.RasterModelGrid), \
               'model_grid must be a Landlab RasterModelGrid'
        self.grid = model_grid
        
        # Assign the node states
        assert (type(initial_node_states) is numpy.ndarray), \
               'initial_node_states must be a Numpy array'
        assert (len(initial_node_states)==self.grid.num_nodes), \
               'length of initial_node_states must equal number of nodes in grid'
        self.node_state = initial_node_states
                 
        # Figure out how many states there are, and make sure the input data
        # are self consistent.
        self.num_node_states = len(node_state_dict)
        self.num_link_states = self.num_node_states*self.num_node_states
        
        assert (type(transition_list) is list), 'transition_list must be a list!'
        assert (transition_list), \
               'Transition list must contain at least one transition'
        for t in transition_list:
            assert (t.from_state < self.num_link_states), \
                   'Transition from_state out of range'
            assert (t.to_state < self.num_link_states), \
                   'Transition to_state out of range'
        
        # Create priority queue for events and next_update array for links
        self.event_queue = []
        self.next_update = self.grid.create_active_link_array_zeros()
    
        # Assign link types from node types
        self.create_link_state_dictionary_and_pair_list()
    
        # Using the grid of node states, figure out all the link states
        self.assign_link_states_from_node_types()
    
        # Create transition data for links
        self.setup_transition_data(transition_list)

        # Put the various transitions on the event queue
        self.push_transitions_to_event_queue()


    def create_link_state_dictionary_and_pair_list(self):
        """
        Creates a dictionary that can be used as a lookup table to find out 
        which link state corresponds to a particular pair of node states. The 
        dictionary keys are 2-element tuples that represent the states of the 
        FROM and TO nodes in the link. The values are integer codes representing 
        the link state numbers.
        """
        self.link_state_dict = {}
        self.cell_pair = []
        k=0
        for fromstate in range(self.num_node_states):
            for tostate in range(self.num_node_states):
                self.link_state_dict[(fromstate,tostate)] = k
                k+=1
                self.cell_pair.append((fromstate,tostate))
    
        if _DEBUG:
            print 
            print 'create_link_state_dictionary_and_pair_list(): dict is:'
            print self.link_state_dict
            print '  and the pair list is:'
            print self.cell_pair


    def assign_link_states_from_node_types(self):
        """
        Assigns a link-state code for each link, and returns a list of these.
        
        Takes lists/arrays of "from" and "to" node IDs for each link, and a 
        dictionary that associates pairs of node states (represented as a 
        2-element tuple) to link states.
        """
        self.link_state = numpy.zeros(self.grid.num_nodes, dtype=int)
    
        k = 0
        for i in range(self.grid.num_active_links):
            node_pair = (self.node_state[self.grid.activelink_fromnode[i]], \
                         self.node_state[self.grid.activelink_tonode[i]])
            #print 'node pair:', node_pair, 'dict:', link_state_dict[node_pair]
            self.link_state[k] = self.link_state_dict[node_pair]
            k += 1
        
        if _DEBUG:
            print 
            print 'assign_link_states_from_node_types(): the link state array is:'
            print link_state


    def setup_transition_data(self, xn_list):
        """
        Using the transition list and the number of link states, creates 
        three arrays that collectively contain data on state transitions:
            n_xn: for each link state, contains the number of transitions out of 
                  that state.
            xn_to: 2D array that records, for each link state and each
                   transition, the new state into which the link transitions.
            xn_rate: 2D array that records, for each link state and each
                     transition, the rate (1/time) of the transition.                 
        """
        # First, create an array that stores the number of possible transitions
        # out of each state.
        self.n_xn = numpy.zeros(self.num_link_states, dtype=int)
        for xn in xn_list:
            self.n_xn[xn.from_state] += 1
        
        # Now, create arrays to hold the "to state" and transition rate for each
        # transition. These arrays are dimensioned N x M where N is the number
        # of states, and M is the maximum number of transitions from a single 
        # state (for example if state 3 could transition either to state 1 or 
        # state 4, and the other states only had one or zero possible 
        # transitions, then the maximum would be 2).
        max_transitions = numpy.max(self.n_xn)
        self.xn_to = numpy.zeros((self.num_link_states, max_transitions), dtype=int)
        self.xn_rate = numpy.zeros((self.num_link_states, max_transitions))
    
        #print n_xn, xn_to, xn_rate
    
        # Populate the "to" and "rate" arrays
        self.n_xn[:] = 0  # reset this and then re-do (inefficient but should work)
        for xn in xn_list:
            #print 'from:',xn.from_state,'to:',xn.to_state,'rate:',xn.rate
            from_state = xn.from_state
            self.xn_to[from_state][self.n_xn[from_state]] = xn.to_state
            self.xn_rate[from_state][self.n_xn[from_state]] = xn.rate
            self.n_xn[from_state] += 1
    
        if _DEBUG:
            print 
            print 'setup_transition_data():'
            print '  n_xn',self.n_xn
            print '  to:',self.xn_to
            print '  rate:',self.xn_rate
    
    
    def get_next_event(self, link, current_state, current_time):
        """
        Returns the next event for link with ID "link", which is in state
        "current state".
    
        If there is only one potential transition out of the current state, a 
        time for the transition is selected at random from an exponential 
        distribution with rate parameter appropriate for this transition.
    
        If there are more than one potential transitions, a transition time is 
        chosen for each, and the smallest of these applied.
    
        Assumes that there is at least one potential transition from the current
        state.
    
        Inputs: link - ID of the link
                current_state - current state code for the link
                current_time - current time in simulation (i.e., time of event
                                just processed)
            
        Returns: an Event() object containing the time, link ID, and type of the
                next transition event at this link.
        """
        assert (self.n_xn[current_state]>0), \
               'must have at least one potential transition'
    
        #rate = self.xn_rate[current_state][0]
        #rate = xn_rate[current_state][0] * (_SURFACE - ...
    
        # Find next event time for each potential transition
        if self.n_xn[current_state]==1:
            xn = self.xn_to[current_state][0]
            next_time = numpy.random.exponential(1.0/self.xn_rate[current_state][0])
        else:
            next_time = _NEVER
            xn = None
            for i in range(self.n_xn[current_state]):
                this_next = numpy.random.exponential(1.0/self.xn_rate[current_state][i])
                if this_next < next_time:
                    next_time = this_next
                    xn = self.xn_to[current_state][i]
    
        # Create and setup event, and return it
        my_event = Event(next_time+current_time, link, xn)
    
        if _DEBUG:
            print 'get_next_event():'
            print '  next_time:',my_event.time
            print '  link:',my_event.link
            print '  xn_to:',my_event.xn_to
    
        return my_event
    
    
    def push_transitions_to_event_queue(self):
    
        if _DEBUG:
            print 'push_transitions_to_event_queue():'
        for i in range(self.num_link_states):
        
            if self.n_xn[self.link_state[i]] > 0:
                #print 'link',i,'has state',link_state[i],'and',n_xn[link_state[i]],'potential transitions'
                event = self.get_next_event(i, self.link_state[i], 0.0)
                heappush(self.event_queue, event)
                self.next_update[i] = event.time
            
            else:
                self.next_update[i] = _NEVER
            
        if _DEBUG:
            print '  push_transitions_to_event_queue(): events in queue are now:'
            for e in self.event_queue:
                print '    next_time:',e.time,'link:',e.link,'xn_to:',e.xn_to
            
            
def update_node_states(fromnode, tonode, node_state, pair, new_link_state):
    """
    Updates the states of the two nodes in the given link.
    """
    
    # Remember the previous state of each node so we can detect whether the 
    # state has changed
    old_fromnode_state = node_state[fromnode]
    old_tonode_state = node_state[tonode]
    
    # Change to the new states
    node_state[fromnode] = pair[new_link_state][0]
    node_state[tonode] = pair[new_link_state][1]
    
    if _DEBUG:
        print 'update_node_states() for',fromnode,'and',tonode
        print '  fromnode was',old_fromnode_state,'and is now',node_state[fromnode]
        print '  tonode was',old_tonode_state,'and is now',node_state[tonode]
    
    return node_state[fromnode]!=old_fromnode_state, \
           node_state[tonode]!=old_tonode_state
           
           
def update_link_state(link, link_state, new_link_state, n_xn, eq, next_update,
                      xn_to, xn_rate, current_time):
    """
    Implements a link transition by updating the current state of the link and
    (if appropriate) choosing the next transition event and pushing it on to the
    event queue.
    
    Inputs:
        link - ID of the link to update
        link_state - array of state codes for all links (MODIFIES)
        new_link_state - code for the new state
        n_xs - array of number of potential transitions from each state
        eq - event queue (MODIFIES)
        next_update - time of next update (transition) at each link (MODIFIES)
        xn_to - array of transitions "to" from each state
        xn_rate - array of transition rate for each transition
    """
    if _DEBUG:
        print
        print 'update_link_state:',xn_to
    link_state[link] = new_link_state
    if n_xn[new_link_state] > 0:
        event = get_next_event(link, new_link_state, n_xn, xn_to, xn_rate, 
                               current_time)
        heappush(eq, event)
        next_update[link] = event.time
    else:
        next_update[link] = _NEVER
            
    if _DEBUG:
        print
        print '  at link',link
        print '  state changed to',link_state[link]
        print '  update time now',next_update[link]
        
            
def do_transition(event, model_grid, node_state, next_update, pair, link_state,
                  n_xn, eq, xn_rate, node_active_links, ls_dict, xn_to, 
                  current_time):
    """
    Implements a state transition. First checks that the transition is still
    valid by comparing the link's next_update time with the corresponding update
    time in the event object.
    
    If the transition is valid, we:
        1) Update the states of the two nodes attached to the link
        2) Update the link's state, choose its next transition, and push it on
           the event queue.
        3) Update the states of the other links attached to the two nodes, 
           choose their next transitions, and push them on the event queue.
           
    Inputs:
        event - Event() object containing the transition data.
        model_grid - ModelGrid() object
        node_state - array of states for each node
        next_update - time of next update for each link
        pair - list of node-state pairs corresponding to each link state
        link_state - array of states for each link
        n_xs - array with number of transitions out of each link state
        eq - event queue
        xn_rate - array with rate of each transition
        node_active_links - list of arrays containing IDs of links connected to
            each node
        ls_dict - dictionary of link-state codes corresponding to each 
            node-state pair
        
    """
    
    if _DEBUG:
        print
        print 'do_transition() for link',event.link
        
    # We'll process the event if its update time matches the one we have 
    # recorded for the link in question. If not, it means that the link has
    # changed state since the event was pushed onto the event queue, and in that
    # case we'll ignore it.
    if event.time == next_update[event.link]:
        
        if _DEBUG:
            print '  event time =',event.time
        
        fromnode = model_grid.activelink_fromnode[event.link]
        tonode = model_grid.activelink_tonode[event.link]
        from_changed, to_changed = update_node_states(fromnode, tonode, 
                                                      node_state,
                                                      pair, event.xn_to)
        update_link_state(event.link, link_state, event.xn_to, n_xn, eq, next_update,
                          xn_to, xn_rate[link_state[event.link]], event.time)

        # Next, when the state of one of the link's nodes changes, we have to
        # update the states of the OTHER links attached to it. This could happen
        # to one or both nodes.
        if from_changed:
            
            if _DEBUG:
                print '    fromnode has changed state, so updating its links'
            
            for link in node_active_links[:,fromnode]:
                
                if link!=-1 and link!=event.link:
                    
                    fromnode = model_grid.activelink_fromnode[link]
                    tonode = model_grid.activelink_tonode[link]
                    current_pair = (node_state[fromnode], node_state[tonode])
                    new_link_state = ls_dict[current_pair]
                    update_link_state(link, link_state, new_link_state, n_xn, eq, 
                                      next_update, xn_to, xn_rate, event.time)

        if to_changed:
            
            if _DEBUG:
                print '    tonode has changed state, so updating its links'
            
            for link in node_active_links[:,tonode]:
                
                if link!=-1 and link!=event.link:
                    
                    fromnode = model_grid.activelink_fromnode[link]
                    tonode = model_grid.activelink_tonode[link]
                    current_pair = (node_state[fromnode], node_state[tonode])
                    new_link_state = ls_dict[current_pair]
                    update_link_state(link, link_state, new_link_state, n_xn, eq, 
                                      next_update, xn_to, xn_rate, event.time)

    elif _DEBUG:
        print '  event time is',event.time,'but update time is',next_update[event.link],'so event will be ignored'
    
    
def main():
    
    # INITIALIZE

    # User-defined parameters
    nr = 400
    nc = 400
    frac_spacing = 32
    plot_interval = 0.25
    next_plot = plot_interval
    run_duration = 5.0
        
    # Create grid and set up boundaries
    mg = RasterModelGrid(nr, nc, 1.0)
    mg.set_inactive_boundaries(True, True, True, True)
    node_active_links = mg.active_node_links()
    
    # Create other data structures
    # ..todo:
    #   model grid should distinguish between making array of all links
    #       vs active links.
    #   should allow data type argument.
    next_update = mg.create_active_link_array_zeros()
    node_type = mg.create_node_array_zeros()
    if _TEST:
        # for testing, make one node be type 1
        node_type[mg.ncols+1] = 1
    else:
        node_type = make_frac_grid(frac_spacing, model_grid=mg)
    #node_type[mg.ncols*2-2] = 1
    #print 'node type:',node_type
    
    # Create priority queue for events
    eq = []
    
    # Assign link types from node types
    num_node_states = len(_NODE_STATE)
    ls_dict, pair = create_link_state_dictionary_and_pair_list(num_node_states)
    #print 'ls dict:',ls_dict
    print 'from and to:',mg.activelink_fromnode, mg.activelink_tonode
    link_state = assign_link_states_from_node_types(mg.activelink_fromnode,
                                                    mg.activelink_tonode,
                                                    node_type, ls_dict)
    print 'ls:',link_state
    num_link_states = len(ls_dict)
    print 'there are',num_link_states,'link states'
    
    # Create a list of transitions
    xn_list = setup_transition_list()
    
    # Create transition data for links
    num_xn, xn_to, xn_rate = setup_transition_data(num_link_states, xn_list)
    print 'num_xn:', num_xn
    
    # Put the various transitions on the event queue
    push_transitions_to_event_queue(eq, link_state, num_xn, xn_to, xn_rate, next_update)
    
    print 'node type:', node_type
    print 'link state', link_state
    print 'next update', next_update
    print 'eq:'
    #while eq:
    #    e = heappop(eq)
    #    print e.time, e.link, e.xn_to

    #return   
            
    current_time = 0.0

    # Plot initial state
    plt.figure()
    imshow_grid(mg, node_type)
    
    # if test, keep track of # of sap cells
    tm = [0.0]
    nsap = [1]
        
        
    # RUN
    
    # Continue until we've run out of either time or events
    while current_time < run_duration and eq:
        
        print 'Current Time = ', current_time
        
        # Pick the next transition event from the event queue
        ev = heappop(eq)
        
        #print 'Event:',ev.time,ev.link,ev.xn_to
        
        do_transition(ev, mg, node_type, next_update, pair, link_state,
                      num_xn, eq, xn_rate, node_active_links, ls_dict, xn_to,
                      current_time)
        
        # Update current time
        current_time = ev.time
        
        # Monitor the number of total sap cells for testing
        tm.append(current_time)
        nsap.append(numpy.sum(node_type))
        
        # Plot
        if current_time >= next_plot:
            plt.figure()
            imshow_grid(mg, node_type)
            next_plot += plot_interval
        
    if _TEST:
        plt.figure()
        plt.plot(tm, nsap, 'o-')
        plt.show()
        
        
        
def setup_transition_list():
    """
    Creates and returns a list of Transition() objects. This is a "custom"
    function in the sense that any particular application is determined by the
    transition rules that are created here.
    """
    xn_list = []
    
    xn_list.append( Transition(1, 3, 1., 'weathering') ) # rock-sap to sap-sap
    xn_list.append( Transition(2, 3, 1., 'weathering') ) # sap-rock to sap-sap
        
    if _DEBUG:
        print
        print 'setup_transition_list(): list has',len(xn_list),'transitions:'
        for t in xn_list:
            print '  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name
        
    return xn_list
    

xn_list = setup_transition_list()
model_grid = RasterModelGrid(3, 5, 1.0)
ns_dict = { 0 : 'rock', 1 : 'saprolite' }
node_state = model_grid.create_node_array_zeros()
node_state[6] = 1
myca = LinkCellularAutomaton(model_grid, ns_dict, xn_list, node_state)

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #main()