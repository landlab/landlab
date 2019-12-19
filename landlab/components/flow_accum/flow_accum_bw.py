#!/usr/env/python

"""
flow_accum_bw.py: Implementation of the Braun & Willet (2012) stack alorithm.

Implementation of Braun & Willett (2012) algorithm for calculating drainage
area and (optionally) water discharge. Assumes each node has only one
downstream receiver. If water discharge is calculated, the result assumes
steady flow (that is, hydrologic equilibrium).

The main public function is::

    a, q, s = flow_accumulation(r)

which takes an array of receiver-node IDs, r (the nodes that "receive" the flow
from a each node; this array would be returned by the flow_routing component's
calc_flowdirs() method). It returns Numpy
arrays with the drainage area (a) and discharge (q) at each node, along with an
array (s) that contains the IDs of the nodes in downstream-to-upstream order.

If you simply want the ordered list by itself, use::

    s = make_ordered_node_array(r)

Created: GT Nov 2013
"""
import numpy
from six.moves import range

from landlab.core.utils import as_id_array

from .cfuncs import _accumulate_bw, _add_to_stack, _make_donors


class _DrainageStack:

    """Implements Braun & Willett's add_to_stack function.

    The _DrainageStack() class implements Braun & Willett's add_to_stack
    function (as a method) and also keeps track of the counter (j) and
    the stack (s). It is used by the make_ordered_node_array() function.
    """

    def __init__(self, delta, D):

        """Initializes the _Drainage_Stack class.

        Initializes the index counter j to zero, creates the stack array
        s, and stores references to delta and D.
        """
        self.j = 0
        self.s = numpy.zeros(len(D), dtype=int)
        self.delta = delta
        self.D = D

    def add_to_stack(self, l):

        """Adds node l to the stack and increments the current index (j).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.components.flow_accum.flow_accum_bw import(
        ... _DrainageStack)
        >>> delta = np.array([ 0,  0,  2,  2,  2,  6,  7,  9, 10, 10, 10])
        >>> D = np.array([0, 2, 1, 4, 5, 7, 6, 3, 8, 9])
        >>> ds = _DrainageStack(delta, D)
        >>> ds.add_to_stack(4)
        >>> ds.s
        array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
        """
        # we invoke cython here to attempt to suppress Python's RecursionLimit
        self.j = _add_to_stack(l, self.j, self.s, self.delta, self.D)


def _make_number_of_donors_array(r):

    """Number of donors for each node.

    Creates and returns an array containing the number of donors for each node.

    Parameters
    ----------
    r : ndarray
        ID of receiver for each node.

    Returns
    -------
    ndarray
        Number of donors for each node.

    Examples
    --------
    The example below is from Braun and Willett (2012); nd corresponds to their
    d_i in Table 1.

    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_bw import(
    ... _make_number_of_donors_array)
    >>> r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8]) - 1
    >>> nd = _make_number_of_donors_array(r)
    >>> nd
    array([0, 2, 0, 0, 4, 1, 2, 1, 0, 0])
    """
    nd = numpy.zeros(r.size, dtype=int)
    max_index = numpy.max(r)
    nd[: (max_index + 1)] = numpy.bincount(r)
    return nd


def _make_delta_array(nd):

    r"""
    Delta array.

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
    The example below is from Braun and Willett (2012), and represents
    \delta_i in their Table 1. Here, the numbers are all one less than in their
    table because here we number indices from 0 rather than 1.

    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_bw import(
    ...     _make_delta_array)
    >>> nd = np.array([0, 2, 0, 0, 4, 1, 2, 1, 0, 0])
    >>> delta = _make_delta_array(nd)
    >>> delta
    array([ 0,  0,  2,  2,  2,  6,  7,  9, 10, 10, 10])
    """
    np = len(nd)
    delta = numpy.zeros(np + 1, dtype=int)
    delta.fill(np)
    delta[-2::-1] -= numpy.cumsum(nd[::-1])
    return delta


def _make_array_of_donors(r, delta):

    """Creates and returns an array containing the IDs of donors for each node.

    Essentially, the array is a series of lists (not in the Python list object
    sense) of IDs for each node. See Braun & Willett (2012) for details.

    The example below is from Braun & Willett (2012), and produces D_i in their
    Table 1 (except that here the ID numbers are one less, because we number
    indices from zero).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum.flow_accum_bw import(
    ... _make_array_of_donors)
    >>> r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
    >>> delta = np.array([ 0,  0,  2,  2,  2,  6,  7,  9, 10, 10, 10])
    >>> D = _make_array_of_donors(r, delta)
    >>> D
    array([0, 2, 1, 4, 5, 7, 6, 3, 8, 9])
    """
    np = len(r)
    w = numpy.zeros(np, dtype=int)
    D = numpy.zeros(np, dtype=int)

    _make_donors(np, w, D, delta, r)

    return D


def make_ordered_node_array(receiver_nodes):

    """Create an array of node IDs that is arranged in order from.

    Creates and returns an array of node IDs that is arranged in order from
    downstream to upstream.

    The lack of a leading underscore is meant to signal that this operation
    could be useful outside of this module!

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum import make_ordered_node_array
    >>> r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
    >>> s = make_ordered_node_array(r)
    >>> s
    array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    """
    node_id = numpy.arange(receiver_nodes.size)
    baselevel_nodes = numpy.where(node_id == receiver_nodes)[0]
    nd = _make_number_of_donors_array(receiver_nodes)
    delta = _make_delta_array(nd)
    D = _make_array_of_donors(receiver_nodes, delta)
    dstack = _DrainageStack(delta, D)
    add_it = dstack.add_to_stack
    for k in baselevel_nodes:
        add_it(k)  # don't think this is a bottleneck, so no C++

    return dstack.s


def find_drainage_area_and_discharge(
    s, r, node_cell_area=1.0, runoff=1.0, boundary_nodes=None
):

    """Calculate the drainage area and water discharge at each node.

    Parameters
    ----------
    s : ndarray of int
        Ordered (downstream to upstream) array of node IDs
    r : ndarray of int
        Receiver IDs for each node
    node_cell_area : float or ndarray
        Cell surface areas for each node. If it's an array, must have same
        length as s (that is, the number of nodes).
    runoff : float or ndarray
        Local runoff rate at each cell (in water depth per time). If it's an
        array, must have same length as s (that is, the number of nodes).
    boundary_nodes: list, optional
        Array of boundary nodes to have discharge and drainage area set to zero.
        Default value is None.
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
    >>> from landlab.components.flow_accum import (
    ...     find_drainage_area_and_discharge)
    >>> r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
    >>> s = np.array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    >>> a, q = find_drainage_area_and_discharge(s, r)
    >>> a
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    >>> q
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    """
    # Number of points
    np = len(s)

    # Initialize the drainage_area and discharge arrays. Drainage area starts
    # out as the area of the cell in question, then (unless the cell has no
    # donors) grows from there. Discharge starts out as the cell's local runoff
    # rate times the cell's surface area.
    drainage_area = numpy.zeros(np, dtype=int) + node_cell_area
    discharge = numpy.zeros(np, dtype=int) + node_cell_area * runoff

    # Optionally zero out drainage area and discharge at boundary nodes
    if boundary_nodes is not None:
        drainage_area[boundary_nodes] = 0
        discharge[boundary_nodes] = 0

    # Call the cfunc to work accumulate from upstream to downstream, permitting
    # transmission losses
    _accumulate_bw(np, s, r, drainage_area, discharge)
    # nodes at channel heads can still be negative with this method, so...
    discharge = discharge.clip(0.0)

    return drainage_area, discharge


def find_drainage_area_and_discharge_lossy(
    s, r, l, loss_function, grid, node_cell_area=1.0, runoff=1.0, boundary_nodes=None
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
    r : ndarray of int
        Receiver node IDs for each node
    l : ndarray of int
        Link to receiver node IDs for each node
    loss_function : Python function(Qw, nodeID, linkID, grid)
        Function dictating how to modify the discharge as it leaves each node.
        nodeID is the current node; linkID is the downstream link, grid is a
        ModelGrid. Returns a float.
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
        Array of boundary nodes to have discharge and drainage area set to zero.
        Default value is None.

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
    >>> from landlab.components.flow_accum import (
    ...     find_drainage_area_and_discharge)
    >>> r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
    >>> s = np.array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    >>> l = np.ones(10, dtype=int)  # dummy
    >>> nodes_wo_outlet = np.array([0, 1, 2, 3, 5, 6, 7, 8, 9])

    >>> def lossfunc(Qw, dummyn, dummyl, dummygrid):
    ...     return 0.5 * Qw
    >>> mg = RasterModelGrid((3, 4))  # some grid big enough to make go
    >>> _ = mg.add_zeros('node', 'surface_water__discharge_loss', dtype=float)
    >>> a, q = find_drainage_area_and_discharge_lossy(s, r, l, lossfunc, mg)
    >>> a
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    >>> q
    array([  1.  ,   2.  ,   1.  ,   1.  ,  3.75,   2.  ,   2.  ,   1.5 ,   1.  ,   1.  ])
    >>> np.allclose(
    ...     mg.at_node['surface_water__discharge_loss'][nodes_wo_outlet],
    ...     0.5*q[nodes_wo_outlet])
    True
    >>> np.isclose(mg.at_node['surface_water__discharge_loss'][4], 0.)
    True

    >>> lossfield = mg.add_ones('node', 'loss_field', dtype=float)
    >>> lossfield *= 0.5
    >>> def lossfunc2(Qw, nodeID, dummyl, grid):
    ...     return grid.at_node['loss_field'][nodeID] * Qw
    >>> a, q = find_drainage_area_and_discharge_lossy(s, r, l, lossfunc2, mg)
    >>> a
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    >>> q
    array([  1.  ,   2.  ,   1.  ,   1.  ,  3.75,   2.  ,   2.  ,   1.5 ,   1.  ,   1.  ])
    >>> np.allclose(
    ...     mg.at_node['surface_water__discharge_loss'][nodes_wo_outlet],
    ...     lossfield[nodes_wo_outlet] * q[nodes_wo_outlet])
    True

    >>> def lossfunc3(Qw, nodeID, dummyl, dummygrid):
    ...     return Qw - 100.  # a huge loss
    >>> a, q = find_drainage_area_and_discharge_lossy(s, r, l, lossfunc3, mg)
    >>> a
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    >>> q
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    """
    # Number of points
    np = len(s)

    # Initialize the drainage_area and discharge arrays. Drainage area starts
    # out as the area of the cell in question, then (unless the cell has no
    # donors) grows from there. Discharge starts out as the cell's local runoff
    # rate times the cell's surface area.
    drainage_area = numpy.zeros(np, dtype=int) + node_cell_area
    discharge = numpy.zeros(np, dtype=int) + node_cell_area * runoff
    # note no loss occurs at a node until the water actually moves along a link

    # Optionally zero out drainage area and discharge at boundary nodes
    if boundary_nodes is not None:
        drainage_area[boundary_nodes] = 0
        discharge[boundary_nodes] = 0

    # Iterate backward through the list, which means we work from upstream to
    # downstream.
    for i in range(np - 1, -1, -1):
        donor = s[i]
        recvr = r[donor]
        lrec = l[donor]
        if donor != recvr:
            drainage_area[recvr] += drainage_area[donor]
            discharge_remaining = numpy.clip(
                loss_function(discharge[donor], donor, lrec, grid), 0.0, float("inf")
            )
            grid.at_node["surface_water__discharge_loss"][donor] = (
                discharge[donor] - discharge_remaining
            )
            discharge[recvr] += discharge_remaining

    return drainage_area, discharge


def flow_accumulation(
    receiver_nodes, node_cell_area=1.0, runoff_rate=1.0, boundary_nodes=None
):

    """Calculate drainage area and (steady) discharge.

    Calculates and returns the drainage area and (steady) discharge at each
    node, along with a downstream-to-upstream ordered list (array) of node IDs.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.flow_accum import flow_accumulation
    >>> r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
    >>> a, q, s = flow_accumulation(r)
    >>> a
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    >>> q
    array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
    >>> s
    array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    """

    s = as_id_array(make_ordered_node_array(receiver_nodes))
    # Note that this ordering of s DOES INCLUDE closed nodes. It really shouldn't!
    # But as we don't have a copy of the grid accessible here, we'll solve this
    # problem as part of route_flow_dn.

    a, q = find_drainage_area_and_discharge(
        s, receiver_nodes, node_cell_area, runoff_rate, boundary_nodes
    )

    return a, q, s


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
