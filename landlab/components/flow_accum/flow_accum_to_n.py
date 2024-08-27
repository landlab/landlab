#!/usr/env/python

"""Short description.

flow_accum_to_n.py: Implementation a route-to-multiple drainage stack alorithm.


Algorithm for route to multiple (N) flow accumulation. Inspiration for data
structures and attempting O(n) efficiency taken from Braun and Willet(2013).

Algorithm constructs drainage area and (optionally) water discharge. Can
handle the case in which each node has more than one downstream receiver.

Computationally, for a grid of the same size this algorithm will take about

    1.5 x (avg number of downstream nodes per cell)
        x (duration of flow_accum_bw for same grid using route-to-one method)

So under route-to-one direction schemes, using the Braun and Willet method is
recommended.

If water discharge is calculated, the result assumes steady flow (that is,
hydrologic equilibrium).

The main public function is::

    a, q, s = flow_accumulation_to_n(r, p)

which takes the following inputs:

    r, an (np, q) array of receiver-node IDs, where np is the total number of
    nodes and q is the maximum number of receivers any node in the grid has.
    This array would be returned by the flow_routing component.

    p, an (np, q) array that identifies the proportion of flow going to each
    receiver. For each q elements along the np axis, sum(p(i, :)) must equal
    1. This array would be returned by the flow_routing component.

It returns Numpy arrays with the drainage area (a) and discharge (q) at each
node, along with an array (s) that contains the IDs of the nodes in downstream-
to-upstream order.

If you simply want the ordered list by itself, use::

    s = make_ordered_node_array_to_n(r, p, b)

Created: KRB Oct 2016 (modified from flow_accumu_bw)
"""
import numpy

from landlab.core.utils import as_id_array

from .cfuncs import _accumulate_to_n
from .cfuncs import _make_donors_to_n


class _DrainageStack_to_n:
    """Implementation of the DrainageStack_to_n class.

    The _DrainageStack_to_n() class implements a set based approach to
    constructing a stack with similar properties to the stack constructed by
    Braun & Willet (2013). It constructs an list, s, of all nodes in the grid
    such that a given node is always located earlier in the list than all
    upstream nodes that contribute to it.

    It is used by the make_ordered_node_array_to_n() function.
    """

    def __init__(self, delta, D, num_receivers):
        """Creates the stack array s and stores references to delta and D.

        Initialization of the _DrainageStack_to_n() class including
        storing delta and D.
        """

        self.num_receivers = num_receivers
        self.s = []
        self.delta = delta
        self.D = D

    def construct__stack(self, nodes):
        """Function to construct the drainage stack.

        Function to add all nodes upstream of a set of base level nodes given
        by list *nodes* in an order
        such that downstream nodes always occur before upstream nodes.

        This function contains the major algorithmic difference between the
        route to 1 method of Braun and Willet (2013) and the route to N method
        presented here.

        Rather than recursively moving up the tributary tree this method uses
        sets test that a node is downstream and add it to the stack. Both
        methods are functionally depth first searches. The method that Braun
        and Willet (2013) implement is optimized given that each node only has
        one receiver. This method is optimized to visit more than one vertex/
        node of the graph at a time.

        An important note: Since sets are un-ordered, we cannot expect the
        stack to be exactly the same each time. It will always put nodes that
        are downstream before those that are upstream, but because it will move
        up multiple branches at the same time, it may put three nodes into the
        stack at the same time that are on different branches of the flow
        network. Because these nodes are in different parts of the network,
        the relative order of them does not matter.

        For example, in the example below, the nodes 1 and 7 must be added
        after 5 but before 2 and 6.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.components.flow_accum.flow_accum_to_n import (
        ...     _DrainageStack_to_n,
        ... )
        >>> delta = np.array([0, 0, 2, 4, 4, 8, 12, 14, 17, 18, 18])
        >>> num_receivers = np.array([2, 2, 2, 2, 1, 1, 2, 2, 2, 2])
        >>> D = np.array([0, 2, 0, 3, 1, 4, 5, 7, 6, 1, 2, 7, 3, 8, 9, 6, 8, 9])
        >>> ds = _DrainageStack_to_n(delta, D, num_receivers)
        >>> ds.construct__stack(4)
        >>> ds.s[0] == 4
        True
        >>> ds.s[1] == 5
        True
        >>> ds.s[9] == 9
        True
        >>> len(set([1, 7]) - set(ds.s[2:4]))
        0
        >>> len(set([2, 6]) - set(ds.s[4:6]))
        0
        >>> len(set([0, 3, 8]) - set(ds.s[6:9]))
        0
        """
        # create base nodes set
        try:
            base = set(nodes)
        except TypeError:
            base = {nodes}

        # instantiate the time keeping variable i, and a variable to keep track
        # of the visit time. Using visit time allows us to itterate through
        # the entire graph and make sure that only put a node in the stack
        # the last time it is visited.

        i = 0
        visit_time = -1 * numpy.ones(self.delta.size - 1)
        num_visits = numpy.zeros(self.delta.size - 1)

        # deal with the first node, which goes to it
        visit_time[list(base)] = i
        num_visits[list(base)] += 1

        i = 1
        visited = set()
        for node_i in base:
            # select the nodes to visit
            visit = set(self.D[self.delta[node_i] : self.delta[node_i + 1]])
            visit = visit - base

            # record the visit time.
            visit_time[list(visit)] = i

            # record that they have been visited.
            num_visits[list(visit)] += 1

            visited.update(list(visit))

        visited = numpy.array(list(visited))
        if visited.size > 0:
            visited_enough = num_visits[visited] == self.num_receivers[visited]
            completed = set(visited[visited_enough])
        else:
            completed = {}
        # recurse through the remainder. Only look above completed nodes,
        # this prevents repeat link walking.
        while len(completed) > 0:
            # increase counter
            i += 1

            visited = set()
            new_completes = set()

            for node_i in completed:
                # select the nodes to visit
                visit = self.D[self.delta[node_i] : self.delta[node_i + 1]]
                # record the visit time.
                visit_time[visit] = i

                # record that they have been visited.
                num_visits[visit] += 1

                # add nodes that have been visited enough times to complete
                # to the upstream stack. We can ignore the rest, they will
                # be re-visited. This should reduce the number of times each
                # link is walked to the number of active links.
                visited_enough = (
                    num_visits[numpy.array(visit)]
                    == self.num_receivers[numpy.array(visit)]
                )

                visited.update(visit)
                new_completes.update(visit[visited_enough])
            completed = new_completes

        # the stack is the argsort of visit time.
        self.s = numpy.argsort(visit_time, kind="stable")


def _make_number_of_donors_array_to_n(r, p):
    """Number of donors for each node.

    Creates and returns an array containing the number of donors for each node.

    Parameters
    ----------
    r : ndarray size (np, q) where r[i,:] gives all receivers of node i. Each
        node recieves flow fom up to q donors.

    p : ndarray size (np, q) where p[i,v] give the proportion of flow going
        from node i to the receiver listed in r[i,v].

    Returns
    -------
    ndarray size (np)
        Number of donors for each node.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_to_n import (
    ...     _make_number_of_donors_array_to_n,
    ... )
    >>> r = np.array(
    ...     [
    ...         [1, 2],
    ...         [4, 5],
    ...         [1, 5],
    ...         [6, 2],
    ...         [4, -1],
    ...         [4, -1],
    ...         [5, 7],
    ...         [4, 5],
    ...         [6, 7],
    ...         [7, 8],
    ...     ]
    ... )
    >>> p = np.array(
    ...     [
    ...         [0.6, 0.4],
    ...         [0.85, 0.15],
    ...         [0.65, 0.35],
    ...         [0.9, 0.1],
    ...         [1.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.75, 0.25],
    ...         [0.55, 0.45],
    ...         [0.8, 0.2],
    ...         [0.95, 0.05],
    ...     ]
    ... )
    >>> nd = _make_number_of_donors_array_to_n(r, p)
    >>> nd
    array([0, 2, 2, 0, 4, 4, 2, 3, 1, 0])
    """
    nd = numpy.zeros(r.shape[0], dtype=int)

    # filter r based on p and flatten
    r_filter_flat = r.flatten()[p.flatten() > 0]

    max_index = numpy.amax(r_filter_flat)

    nd[: (max_index + 1)] = numpy.bincount(r_filter_flat)
    return nd


def _make_delta_array_to_n(nd):
    r"""
    Function to create the delta array.

    Creates and returns the "delta" array, which is a list containing, for each
    node, the array index where that node's donor list begins.

    Parameters
    ----------
    nd : ndarray of int
        Number of donors for each node

    Returns
    -------
    ndarray of int
        Delta array

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_to_n import _make_delta_array_to_n
    >>> nd = np.array([0, 2, 2, 0, 4, 4, 2, 3, 1, 0])
    >>> delta = _make_delta_array_to_n(nd)
    >>> delta
    array([ 0,  0,  2,  4,  4,  8,  12,  14, 17, 18, 18])
    >>> sum(nd) == max(delta)
    True
    """
    nt = sum(nd)
    np = len(nd)
    delta = numpy.zeros(np + 1, dtype=int)
    delta.fill(nt)
    delta[-2::-1] -= numpy.cumsum(nd[::-1])

    return delta


def _make_array_of_donors_to_n(r, p, delta):
    """Creates and returns an array containing the IDs of donors for each node.

    Essentially, the array is a series of lists (not in the Python list object
    sense) of IDs for each node. See Braun & Willett (2012) for details.

    The example below is from Braun & Willett (2012), and produces D_i in their
    Table 1 (except that here the ID numbers are one less, because we number
    indices from zero).

    Vectorized - inefficiently! - DEJH, 5/20/14

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_to_n import (
    ...     _make_array_of_donors_to_n,
    ... )
    >>> r = np.array(
    ...     [
    ...         [1, 2],
    ...         [4, 5],
    ...         [1, 5],
    ...         [6, 2],
    ...         [4, -1],
    ...         [4, -1],
    ...         [5, 7],
    ...         [4, 5],
    ...         [6, 7],
    ...         [7, 8],
    ...     ]
    ... )
    >>> p = np.array(
    ...     [
    ...         [0.6, 0.4],
    ...         [0.85, 0.15],
    ...         [0.65, 0.35],
    ...         [0.9, 0.1],
    ...         [1.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.75, 0.25],
    ...         [0.55, 0.45],
    ...         [0.8, 0.2],
    ...         [0.95, 0.05],
    ...     ]
    ... )
    >>> delta = np.array([0, 0, 2, 4, 4, 8, 12, 14, 17, 18, 18])
    >>> D = _make_array_of_donors_to_n(r, p, delta)
    >>> D
    array([0, 2, 0, 3, 1, 4, 5, 7, 6, 1, 2, 7, 3, 8, 9, 6, 8, 9])
    """
    np = r.shape[0]
    q = r.shape[1]
    nt = delta[-1]

    w = numpy.zeros(np, dtype=int)
    D = numpy.zeros(nt, dtype=int)

    _make_donors_to_n(np, q, w, D, delta, r, p)

    return D


def make_ordered_node_array_to_n(
    receiver_nodes, receiver_proportion, nd=None, delta=None, D=None
):
    """Create an array of node IDs.

    Creates and returns an array of node IDs that is arranged in order from
    downstream to upstream.

    The lack of a leading underscore is meant to signal that this operation
    could be useful outside of this module!

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_to_n import (
    ...     make_ordered_node_array_to_n,
    ... )
    >>> r = np.array(
    ...     [
    ...         [1, 2],
    ...         [4, 5],
    ...         [1, 5],
    ...         [6, 2],
    ...         [4, -1],
    ...         [4, -1],
    ...         [5, 7],
    ...         [4, 5],
    ...         [6, 7],
    ...         [7, 8],
    ...     ]
    ... )
    >>> p = np.array(
    ...     [
    ...         [0.6, 0.4],
    ...         [0.85, 0.15],
    ...         [0.65, 0.35],
    ...         [0.9, 0.1],
    ...         [1.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.75, 0.25],
    ...         [0.55, 0.45],
    ...         [0.8, 0.2],
    ...         [0.95, 0.05],
    ...     ]
    ... )
    >>> s = make_ordered_node_array_to_n(r, p)
    >>> s[0] == 4
    True
    >>> s[1] == 5
    True
    >>> s[9] == 9
    True
    >>> len(set([1, 7]) - set(s[2:4]))
    0
    >>> len(set([2, 6]) - set(s[4:6]))
    0
    >>> len(set([0, 3, 8]) - set(s[6:9]))
    0
    """
    node_id = numpy.arange(receiver_nodes.shape[0])
    baselevel_nodes = numpy.where(node_id == receiver_nodes[:, 0])[0]
    if nd is None:
        nd = _make_number_of_donors_array_to_n(receiver_nodes, receiver_proportion)
    if delta is None:
        delta = _make_delta_array_to_n(nd)
    if D is None:
        D = _make_array_of_donors_to_n(receiver_nodes, receiver_proportion, delta)

    num_receivers = numpy.sum(receiver_nodes >= 0, axis=1)

    dstack = _DrainageStack_to_n(delta, D, num_receivers)
    construct_it = dstack.construct__stack

    construct_it(baselevel_nodes)  # don't think this is a bottleneck, so no C++
    return dstack.s


def find_drainage_area_and_discharge_to_n(
    s, r, p, node_cell_area=1.0, runoff=1.0, boundary_nodes=None
):
    """Calculate the drainage area and water discharge at each node.

    Parameters
    ----------
    s : ndarray of int
        Ordered (downstream to upstream) array of node IDs
    r : ndarray size (np, q) where r[i, :] gives all receivers of node i. Each
        node recieves flow fom up to q donors.
    p : ndarray size (np, q) where p[i, v] give the proportion of flow going
        from node i to the receiver listed in r[i, v].
    node_cell_area : float or ndarray
        Cell surface areas for each node. If it's an array, must have same
        length as s (that is, the number of nodes).
    runoff : float or ndarray
        Local runoff rate at each cell (in water depth per time). If it's an
        array, must have same length as s (that is, the number of nodes).
        runoff *is* permitted to be negative, in which case it performs as a
        transmission loss.
    boundary_nodes: list, optional
        Array of boundary nodes to have discharge and drainage area set to
        zero. Default value is None.

    Returns
    -------
    tuple of ndarray
        drainage area and discharge

    Notes
    -----
    -  If node_cell_area not given, the output drainage area is equivalent
       to the number of nodes/cells draining through each point, including
       the local node itself.
    -  Give node_cell_area as a scalar when using a regular raster grid.
    -  If runoff is not given, the discharge returned will be the same as
       drainage area (i.e., drainage area times unit runoff rate).
    -  If using an unstructured Landlab grid, make sure that the input
       argument for node_cell_area is the cell area at each NODE rather than
       just at each CELL. This means you need to include entries for the
       perimeter nodes too. They can be zeros.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_to_n import (
    ...     find_drainage_area_and_discharge_to_n,
    ... )
    >>> r = np.array(
    ...     [
    ...         [1, 2],
    ...         [4, 5],
    ...         [1, 5],
    ...         [6, 2],
    ...         [4, -1],
    ...         [4, -1],
    ...         [5, 7],
    ...         [4, 5],
    ...         [6, 7],
    ...         [7, 8],
    ...     ]
    ... )
    >>> p = np.array(
    ...     [
    ...         [0.6, 0.4],
    ...         [0.85, 0.15],
    ...         [0.65, 0.35],
    ...         [0.9, 0.1],
    ...         [1.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.75, 0.25],
    ...         [0.55, 0.45],
    ...         [0.8, 0.2],
    ...         [0.95, 0.05],
    ...     ]
    ... )
    >>> s = np.array([4, 5, 1, 7, 2, 6, 0, 8, 3, 9])
    >>> a, q = find_drainage_area_and_discharge_to_n(s, r, p)
    >>> a.round(4)
    array([  1.    ,   2.575 ,   1.5   ,   1.    ,  10.    ,   5.2465,
             2.74  ,   2.845 ,   1.05  ,   1.    ])
    >>> q.round(4)
    array([  1.    ,   2.575 ,   1.5   ,   1.    ,  10.    ,   5.2465,
             2.74  ,   2.845 ,   1.05  ,   1.    ])
    """
    # Number of points
    np = r.shape[0]
    q = r.shape[1]

    # Initialize the drainage_area and discharge arrays. Drainage area starts
    # out as the area of the cell in question, then (unless the cell has no
    # donors) grows from there. Discharge starts out as the cell's local runoff
    # rate times the cell's surface area.
    drainage_area = numpy.zeros(np) + node_cell_area
    discharge = numpy.zeros(np) + node_cell_area * runoff

    # Optionally zero out drainage area and discharge at boundary nodes
    if boundary_nodes is not None:
        drainage_area[boundary_nodes] = 0
        discharge[boundary_nodes] = 0

    # Call the cfunc to work accumulate from upstream to downstream, permitting
    # transmission losses
    _accumulate_to_n(np, q, s, r, p, drainage_area, discharge)
    # nodes at channel heads can still be negative with this method, so...
    discharge = discharge.clip(0.0)

    return drainage_area, discharge


def find_drainage_area_and_discharge_to_n_lossy(
    s,
    r,
    link_to_receiver,
    p,
    loss_function,
    grid,
    node_cell_area=1.0,
    runoff=1.0,
    boundary_nodes=None,
):
    """Calculate the drainage area and water discharge at each node, permitting
    discharge to fall (or gain) as it moves downstream according to some
    function. Note that only transmission creates loss, so water sourced
    locally within a cell is always retained. The loss on each link is recorded
    in the 'surface_water__discharge_loss' link field on the grid; ensure this
    exists before running the function.

    Parameters
    ----------
    s : ndarray of int
        Ordered (downstream to upstream) array of node IDs
    r : ndarray size (np, q) where r[i, :] gives all receivers of node i. Each
        node receives flow fom up to q donors.
    link_to_receiver : ndarray size (np, q) where l[i, :] gives all links to receivers of
        node i.
    p : ndarray size (np, q) where p[i, v] give the proportion of flow going
        from node i to the receiver listed in r[i, v].
    loss_function : Python function(Qw, nodeID, linkID)
        Function dictating how to modify the discharge as it leaves each node.
        nodeID is the current node; linkID is the downstream link. Returns a
        float.
    grid : Landlab ModelGrid (or None)
        A grid to enable spatially variable parameters to be used in the loss
        function. If no spatially resolved parameters are needed, this can be
        a dummy variable, e.g., None.
    node_cell_area : float or ndarray
        Cell surface areas for each node. If it's an array, must have same
        length as s (that is, the number of nodes).
    runoff : float or ndarray
        Local runoff rate at each cell (in water depth per time). If it's an
        array, must have same length as s (that is, the number of nodes).
    boundary_nodes: list, optional
        Array of boundary nodes to have discharge and drainage area set to
        zero. Default value is None.

    Returns
    -------
    tuple of ndarray
        drainage area and discharge

    Notes
    -----
    -  If node_cell_area not given, the output drainage area is equivalent
       to the number of nodes/cells draining through each point, including
       the local node itself.
    -  Give node_cell_area as a scalar when using a regular raster grid.
    -  If runoff is not given, the discharge returned will be the same as
       drainage area (i.e., drainage area times unit runoff rate).
    -  If using an unstructured Landlab grid, make sure that the input
       argument for node_cell_area is the cell area at each NODE rather than
       just at each CELL. This means you need to include entries for the
       perimeter nodes too. They can be zeros.
    -  Loss cannot go negative.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_accum.flow_accum_to_n import (
    ...     find_drainage_area_and_discharge_to_n_lossy,
    ... )
    >>> r = np.array([[1, 2], [3, -1], [3, 1], [3, -1]])
    >>> p = np.array([[0.5, 0.5], [1.0, 0.0], [0.2, 0.8], [1.0, 0.0]])
    >>> s = np.array([3, 1, 2, 0])
    >>> l = np.ones_like(r, dtype=int)  # dummy

    Make here a grid that contains (too many!) links holding values for loss.
    We're only going to use the first 4 links, but illustrates the use of the
    grid for link input.

    >>> mg = RasterModelGrid((3, 3))
    >>> _ = mg.add_zeros("node", "surface_water__discharge_loss", dtype=float)
    >>> lossy = mg.add_ones("lossy", at="link", dtype=float)
    >>> lossy *= 0.5
    >>> def lossfunc(Qw, dummyn, linkID, grid):
    ...     return grid.at_link["lossy"][linkID] * Qw
    ...

    >>> a, q = find_drainage_area_and_discharge_to_n_lossy(s, r, l, p, lossfunc, mg)
    >>> a
    array([1. , 2.7, 1.5, 4. ])
    >>> q
    array([1.  , 1.75, 1.25, 2.  ])
    >>> np.allclose(mg.at_node["surface_water__discharge_loss"][:3], 0.5 * q[:3])
    True

    Note by definition no loss is occuring at the outlet node, as there are no
    nodes downstream.

    Final example of total transmission loss:

    >>> def lossfunc(Qw, dummyn, dummyl, dummygrid):
    ...     return Qw - 100.0  # huge loss
    ...
    >>> a, q = find_drainage_area_and_discharge_to_n_lossy(s, r, l, p, lossfunc, mg)
    >>> a
    array([1. , 2.7, 1.5, 4. ])
    >>> q
    array([1., 1., 1., 1.])
    """
    # Number of points
    np = r.shape[0]
    q = r.shape[1]

    # Initialize the drainage_area and discharge arrays. Drainage area starts
    # out as the area of the cell in question, then (unless the cell has no
    # donors) grows from there. Discharge starts out as the cell's local runoff
    # rate times the cell's surface area.
    drainage_area = numpy.zeros(np) + node_cell_area
    discharge = numpy.zeros(np) + node_cell_area * runoff

    # grab the field to ouput loss to

    # Optionally zero out drainage area and discharge at boundary nodes
    if boundary_nodes is not None:
        drainage_area[boundary_nodes] = 0
        discharge[boundary_nodes] = 0

    # Iterate backward through the list, which means we work from upstream to
    # downstream.
    for i in range(np - 1, -1, -1):
        donor = s[i]
        for v in range(q):
            recvr = r[donor, v]
            lrec = link_to_receiver[donor, v]
            proportion = p[donor, v]
            if proportion > 0 and donor != recvr:
                drainage_area[recvr] += proportion * drainage_area[donor]
                discharge_head = proportion * discharge[donor]
                discharge_remaining = numpy.clip(
                    loss_function(discharge_head, donor, lrec, grid),
                    0.0,
                    float("inf"),
                )
                grid.at_node["surface_water__discharge_loss"][donor] += (
                    discharge_head - discharge_remaining
                )
                discharge[recvr] += discharge_remaining

    return drainage_area, discharge


def flow_accumulation_to_n(
    receiver_nodes,
    receiver_proportions,
    node_cell_area=1.0,
    runoff_rate=1.0,
    boundary_nodes=None,
):
    """Calculate drainage area and (steady) discharge.

    Calculates and returns the drainage area and (steady) discharge at each
    node, along with a downstream-to-upstream ordered list (array) of node IDs.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_to_n import flow_accumulation_to_n
    >>> r = np.array(
    ...     [
    ...         [1, 2],
    ...         [4, 5],
    ...         [1, 5],
    ...         [6, 2],
    ...         [4, -1],
    ...         [4, -1],
    ...         [5, 7],
    ...         [4, 5],
    ...         [6, 7],
    ...         [7, 8],
    ...     ]
    ... )
    >>> p = np.array(
    ...     [
    ...         [0.6, 0.4],
    ...         [0.85, 0.15],
    ...         [0.65, 0.35],
    ...         [0.9, 0.1],
    ...         [1.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.75, 0.25],
    ...         [0.55, 0.45],
    ...         [0.8, 0.2],
    ...         [0.95, 0.05],
    ...     ]
    ... )
    >>> a, q, s = flow_accumulation_to_n(r, p)
    >>> a.round(4)
    array([  1.    ,   2.575 ,   1.5   ,   1.    ,  10.    ,   5.2465,
             2.74  ,   2.845 ,   1.05  ,   1.    ])
    >>> q.round(4)
    array([  1.    ,   2.575 ,   1.5   ,   1.    ,  10.    ,   5.2465,
             2.74  ,   2.845 ,   1.05  ,   1.    ])
    >>> s[0] == 4
    True
    >>> s[1] == 5
    True
    >>> s[9] == 9
    True
    >>> len(set([1, 7]) - set(s[2:4]))
    0
    >>> len(set([2, 6]) - set(s[4:6]))
    0
    >>> len(set([0, 3, 8]) - set(s[6:9]))
    0
    """

    assert (
        receiver_nodes.shape == receiver_proportions.shape
    ), "r and p arrays are not the same shape"

    s = as_id_array(make_ordered_node_array_to_n(receiver_nodes, receiver_proportions))
    # Note that this ordering of s DOES INCLUDE closed nodes. It really
    # shouldn't!
    # But as we don't have a copy of the grid accessible here, we'll solve this
    # problem as part of route_flow_dn.

    a, q = find_drainage_area_and_discharge_to_n(
        s,
        receiver_nodes,
        receiver_proportions,
        node_cell_area,
        runoff_rate,
        boundary_nodes,
    )

    return a, q, s


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
