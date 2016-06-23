# -*- coding: utf-8 -*-
"""
hex_lattice_tectonicizer.py.

Models discrete normal-fault offset on a 2D hex lattice with a rectangular
shape and with one orientation of the nodes being vertical.

The intention here is to use a particle (LCA) model to represent the evolution
of a 2D hillslope, with the hex_lattice_tectonicizer serving to shift the nodes
either upward (simple vertical uplift relative to baselevel), or up and
sideways (representing motion on a fault plane).

Created on Mon Nov 17 08:01:49 2014

@author: gtucker
"""

from landlab import HexModelGrid
from numpy import amax, zeros, arange, array, sqrt

_DEFAULT_NUM_ROWS = 5
_DEFAULT_NUM_COLS = 5
_TAN60 = 1.732


class HexLatticeTectonicizer(object):
    """Handles tectonics and baselevel for CellLab-CTS models.

    This is the base class from which classes to represent particular
    baselevel/fault geometries are derived.

    Examples
    --------
    >>> hlt = HexLatticeTectonicizer()
    >>> hlt.grid.number_of_nodes
    25
    >>> hlt.nr
    5
    >>> hlt.nc
    5
    """

    def __init__(self, grid=None, node_state=None, propid=None, prop_data=None,
                 prop_reset_value=None):
        """
        Create and initialize a HexLatticeTectonicizer.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid(6, 6, shape='rect')
        >>> hlt = HexLatticeTectonicizer()
        >>> hlt.grid.number_of_nodes
        25
        >>> hlt.nr
        5
        >>> hlt.nc
        5
        """
        # If needed, create grid
        if grid is None:
            num_rows = _DEFAULT_NUM_ROWS
            num_cols = _DEFAULT_NUM_COLS
            self.grid = HexModelGrid(num_rows, num_cols, dx=1.0,
                                     orientation='vertical',
                                     shape='rect', reorient_links=True)
        else:
            # Make sure caller passed the right type of grid
            assert (grid.orientation == 'vertical'), \
                   'Grid must have vertical orientation'

            # Keep a reference to the grid
            self.grid = grid

        # If needed, create node-state grid
        if node_state is None:
            self.node_state = self.grid.add_zeros('node', 'node_state_map',
                                                  dtype=int)
        else:
            self.node_state = node_state

        # Remember the # of rows and cols
        self.nr = self.grid.number_of_node_rows
        self.nc = self.grid.number_of_node_columns

        # propid should be either a reference to a CA model's "property id"
        # array, or None
        self.propid = propid
        self.prop_data = prop_data
        self.prop_reset_value = prop_reset_value


class LatticeNormalFault(HexLatticeTectonicizer):
    """Handles normal-fault displacement in CellLab-CTS models.

    Represents a 60 degree, left-dipping normal fault, and handles discrete
    offsets for a hex grid that has vertical columns and a rectangular shape.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import HexModelGrid
    >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
    >>> pid = np.arange(25, dtype=int)
    >>> pdata = np.arange(25)
    >>> ns = np.arange(25, dtype=int)
    >>> grid = HexModelGrid(5, 5, 1.0, orientation='vertical', shape='rect', reorient_links=True)
    >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
    >>> lnf.first_fw_col
    1
    >>> lnf.num_fw_rows
    array([0, 1, 3, 4, 5])
    >>> lnf.incoming_node
    array([1, 3, 4, 6])
    >>> lnf.outgoing_node
    array([12, 17, 19, 22])

    >>> pid = arange(16, dtype=int)
    >>> ns = arange(16, dtype=int)
    >>> pdata = arange(16)
    >>> grid = HexModelGrid(4, 4, 1.0, orientation='vertical', shape='rect', reorient_links=True)
    >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
    >>> lnf.num_fw_rows
    array([0, 1, 3, 4])
    >>> lnf.incoming_node
    array([1, 2, 5])
    >>> lnf.outgoing_node
    array([ 7, 11, 15])
    >>> lnf.do_offset(rock_state=16)
    >>> ns
    array([ 0, 16, 16, 16,  4, 16,  6,  1,  8,  2, 10,  5, 12, 13, 14,  9])
    >>> lnf.propid
    array([ 0,  7, 11,  3,  4, 15,  6,  1,  8,  2, 10,  5, 12, 13, 14,  9])

    >>> pid = arange(20, dtype=int)
    >>> ns = arange(20, dtype=int)
    >>> pdata = arange(20)
    >>> grid = HexModelGrid(4, 5, 1.0, orientation='vertical', shape='rect', reorient_links=True)
    >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
    >>> lnf.incoming_node
    array([1, 3, 4, 6])
    >>> lnf.outgoing_node
    array([12, 14, 17, 19])
    >>> lnf.do_offset(rock_state=20)
    >>> ns
    array([ 0, 20, 20, 20, 20,  5, 20, 20,  8,  1, 10,  3,  4, 13,  6, 15, 16,
            9, 18, 11])
    >>> lnf.propid
    array([ 0, 12,  2, 14, 17,  5, 19,  7,  8,  1, 10,  3,  4, 13,  6, 15, 16,
            9, 18, 11])
    """

    def __init__(self, fault_x_intercept=0.0, grid=None, node_state=None,
                 propid=None, prop_data=None, prop_reset_value=None):
        """
        Create and initialize a LatticeNormalFault object.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ns = np.arange(25, dtype=int)
        >>> grid = HexModelGrid(5, 5, 1.0, orientation='vertical', shape='rect', reorient_links=True)
        >>> lnf = LatticeNormalFault(-0.01, grid, ns, pid, pdata, 0.0)
        >>> lnf.first_fw_col
        0
        >>> lnf.num_fw_rows
        array([1, 2, 4, 5, 5])
        >>> lnf.incoming_node
        array([0, 1, 3, 4, 6])
        >>> lnf.outgoing_node
        array([12, 17, 19, 22, 24])

        >>> pid = np.arange(16, dtype=int)
        >>> pdata = np.arange(16)
        >>> ns = np.arange(16, dtype=int)
        >>> grid = HexModelGrid(4, 4, 1.0, orientation='vertical', shape='rect', reorient_links=True)
        >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
        >>> lnf.first_fw_col
        1
        >>> lnf.num_fw_rows
        array([0, 1, 3, 4])
        >>> lnf.incoming_node
        array([1, 2, 5])
        >>> lnf.outgoing_node
        array([ 7, 11, 15])

        >>> pid = np.arange(45, dtype=int)
        >>> pdata = np.arange(45)
        >>> ns = np.arange(45, dtype=int)
        >>> grid = HexModelGrid(5, 9, 1.0, orientation='vertical', shape='rect', reorient_links=True)
        >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
        >>> lnf.first_fw_col
        1
        >>> lnf.num_fw_rows
        array([0, 1, 3, 4, 5, 5, 5, 5, 5])
        >>> lnf.incoming_node
        array([ 1,  2,  3,  5,  6,  7,  8, 10, 11, 12])
        >>> lnf.outgoing_node
        array([22, 31, 33, 34, 35, 38, 39, 40, 43, 44])
        """
        # Do the base class init
        super(LatticeNormalFault, self).__init__(grid, node_state, propid,
                                                 prop_data, prop_reset_value)
        # Set up data structures:
        #   Make sure the footwall location is such that the fault actually
        #   cuts across the grid. This means the x intercept has to be, at
        #   the very least, no smaller than the biggest x-coordinate, and if
        #   there is an even number of columns, it must be smaller than that
        #   number minus 1/tangent 60 degrees (0.57735)
        assert fault_x_intercept < (amax(self.grid.node_x) - 0.57735), 'err'

        #   Figure out which nodes are and are not within the footwall
        in_footwall = (self.grid.node_y < _TAN60 * (self.grid.node_x -
                       fault_x_intercept))

        # Helpful to have an array of node IDs for the bottom full row. Because
        # the nodes in the bottom row are staggered in a vertical, rectangular
        # hex grid, the IDs go: 0, M, 1, M+1, 2, M+2, ... etc., where M is half
        # the number of columns, rounded up (so, for example, 3 for a 5- or 6-
        # column grid, etc.)
        half_num_cols = (self.nc + 1) // 2
        bottom_row_node_id = (arange(self.nc) // 2 +
                              (arange(self.nc) % 2) * half_num_cols)

        # Also useful to remember the number of even-numbered columns
        self.n_even_cols = (self.nc + 1) // 2

        #   Find the first of the bottom-row nodes that lies in the footwall.
        # This loop exploits the fact that nodes are numbered in an order
        # sorted by x then y, and that the bottom row is staggered, with node
        # zero being "low", like: ,',', etc.
        self.first_fw_col = 0
        n = 0
        while not in_footwall[bottom_row_node_id[n]]:
            n += 1
            self.first_fw_col += 1
            assert n < self.nc, 'overflow in loop'

        #   Remember the number of footwall rows in each column
        self.num_fw_rows = zeros(self.nc, dtype=int)
        for c in range(self.nc):
            current_row = 0
            while (current_row < self.nr and
                   in_footwall[bottom_row_node_id[c] + self.nc*current_row]):
                self.num_fw_rows[c] += 1
                current_row += 1

        # If we're handling properties and property IDs, we need to do some
        # setup
        if self.propid is not None:

            # We want to remember the node IDs of two sets of nodes: those
            # whose contents will vanish off the right-hand side (and possibly
            # the top) of the grid with each offset step, and those that gain
            # new ("rock") contents with each such step.
            #
            # First, let's find the latter set, which we'll call
            # "incoming_node" because material is flowing into these nodes from
            # below. To do this, we'll go column-by-column, starting with
            # the first column that has nodes in the footwall, and going to
            # the next-to-last column. Even-numbered columns have *two*
            # incoming nodes at their base (if the footwall reaches that high),
            # whereas odd-numbered columns have only one. Note that we'll start
            # with a list (so we can append), and then convert to an array
            # that belongs to this class.
            incoming_node_list = []
            lower_right_node = (self.nc % 2) * (self.nc // 2) + \
                               ((self.nc + 1) % 2) * (self.nc - 1)
            for n in range(0, self.nc + (self.nc + 1) // 2):
                if in_footwall[n] and ((n - lower_right_node) % self.nc) != 0:
                    incoming_node_list.append(n)

            # Convert to a numpy array that belongs to the class (so we can
            # use it in the do_offset() method).
            self.incoming_node = array(incoming_node_list, dtype=int)

            # Next, find the IDs the "outgoing" nodes. There will always be
            # some of these along the right side of the grid. Depending on the
            # height of the grid and the position of the fault, there may or
            # may not also be some along the top.
            #
            # To find out which nodes will be exiting the grid, we use
            # geometry: simply find out which nodes are (1) within the footwall
            # and (2) the tectonic offset would take them beyond the grid
            # perimeter. We already have information about the first condition.
            # For the second, let's start by defining the grid edges. The
            # y coordinates of the top full row of nodes will either be NR - 1
            # (even-numbered columns) or NR - 1/2. So we'll take as a
            # reasonable boundary line NR - 1/4: any node that would move above
            # this point is definitely out of the grid. Note that the vertical
            # offset of nodes during each earthquake will be 1.5, so if we were
            # to offset a top-row node, it would move to y = NR + 1/2 or
            # y = NR + 1. Odd-numbered columns in the next-from-top row will
            # move to y = NR, which is out of bounds, so we want to flag these
            # too, and therefore need our cutoff below this. Even-numbered
            # columns in the next-from-top row will end up at y = NR - 1/2,
            # which is within the grid. So our cutoff must be between NR - 1/2
            # and NR. Hence the choise of NR - 1/4 as the effective "top" of
            # the grid.
            top_grid_edge = self.nr - 0.25

            # The right-hand edge of the grid is a simpler case, because our
            # grid is vertically oriented and there is no stagger on the right
            # and left sides. So we simply make it half a column beyond the
            # right-most column. Column width is sqrt(3), the last column is
            # NC - 1, so the right edge y coordinate is (NC - 1/2) x sqrt(3)/2
            right_grid_edge = (self.nc - 0.5) * (sqrt(3.0) / 2.0)

            # To save a repeated calculation in a loop, we'll find a threshold
            # x and y coordinate beyond which any node offset would take them
            # off the grid.
            x_threshold = right_grid_edge - (sqrt(3.0) / 2.0)
            y_threshold = top_grid_edge - 1.5

            # Actually it turns out there is a third criterion. One or two
            # nodes in the lower-right corner could be counted as both
            # "incoming" (they're at the bottom of the grid) AND outgoing
            # (they're on the right-hand side). We ignored these in setting up
            # incoming, so we should ignore them for outgoing too. This is
            # easy: any nodes on the right side (x > x_threshold) that are also
            # near the bottom (y < 1.25) should be ignored.

            # Now march through all nodes, placing those on the list that meet
            # our criteria. Yes, it's slow to do this as a Python loop, but we
            # only do it once.
            outgoing_node_list = []
            for n in range(self.grid.number_of_nodes):
                if (((self.grid.node_x[n] > x_threshold and
                      self.grid.node_y[n] > 1.25) or
                     self.grid.node_y[n] > y_threshold) and in_footwall[n]):
                    outgoing_node_list.append(n)

            # Finally, convert the outgoing node list to an array stored in
            # this object
            self.outgoing_node = array(outgoing_node_list, dtype=int)

    def do_offset(self, rock_state=1):
        """Apply 60-degree normal-fault offset.

        Offset is applied to a hexagonal grid with vertical node orientation
        and rectangular arrangement of nodes.

        Parameters
        ----------
        rock_state : int
            State code to apply to new cells introduced along bottom row.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
        >>> from landlab import HexModelGrid
        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ns = np.arange(25, dtype=int)
        >>> grid = HexModelGrid(5, 5, 1.0, orientation='vertical', shape='rect', reorient_links=True)
        >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
        >>> lnf.do_offset(rock_state=25)
        >>> ns
        array([ 0, 25, 25, 25, 25,  5, 25, 25,  8,  1, 10,  3,  4, 13,  6, 15, 16,
                9, 18, 11, 20, 21, 14, 23, 24])
        >>> lnf.propid
        array([ 0, 12,  2, 17, 19,  5, 22,  7,  8,  1, 10,  3,  4, 13,  6, 15, 16,
                9, 18, 11, 20, 21, 14, 23, 24])
        """

        # If we need to shift the property ID numbers, we'll first need to
        # record the property IDs in those nodes that are about to "shift off
        # the grid" (or rather, their contents will shift) due to tectonic
        # motion. We'll call these "nodes to replace".
        if self.propid is not None:

            # Now, remember the property IDs of the nodes along the right side
            # and possibly top that are about to shift off the grid
            propids_for_incoming_nodes = self.propid[self.outgoing_node]
            #print 'pid for new base:',propids_for_incoming_nodes

        # We go column-by-column, starting from the right side
        for c in range(self.grid.number_of_node_columns - 1,
                       self.first_fw_col - 1, -1):

            # Odd-numbered rows are shifted up in the hexagonal, vertically
            # oriented lattice
            row_offset = 2 - (c % 2)

            # Number of base nodes in the footwall in this column (1 or 2).
            n_base_nodes = min(self.num_fw_rows[c], row_offset)

            # ID of the bottom footwall node in this column
            bottom_node = (c // 2) + ((c % 2) * self.n_even_cols)

            # The bottom 1 or 2 nodes in this column are set to rock
            self.node_state[bottom_node] = rock_state
            if n_base_nodes == 2:
                self.node_state[bottom_node+self.nc] = rock_state

            # "indices" here contains the array indices of those nodes in this
            # column that are to be replaced by the ones in the column to the
            # left and down one or two nodes. We do this replacement if indices
            # contains any data.
            first_repl = bottom_node + n_base_nodes * self.nc
            last_repl = first_repl + (self.num_fw_rows[c] - \
                        (n_base_nodes + 1)) * self.nc
            indices = arange(first_repl, last_repl + 1, self.nc)

            if len(indices) > 0:
                offset = (self.nc + ((self.nc + 1) // 2) +
                          ((c + 1) % 2) * ((self.nc + 1) % 2))
                self.node_state[indices] = self.node_state[indices-offset]
                if self.propid is not None:
                    self.propid[indices] = self.propid[indices-offset]

        if self.propid is not None:
            self.propid[self.incoming_node] = propids_for_incoming_nodes
            self.prop_data[self.propid[self.incoming_node]] = self.prop_reset_value

        if self.first_fw_col==0:
            self.node_state[:self.n_footwall_rows[0]] = rock_state


class LatticeUplifter(HexLatticeTectonicizer):
    """Handles vertical uplift of interior (not edges) for a hexagonal lattice
    with vertical node orientation and rectangular node arrangement.
    """
    def __init__(self, grid=None, node_state=None, propid=None, prop_data=None, prop_reset_value=None):
        """
        Create and initialize a LatticeUplifter

        Examples
        --------
        >>> lu = LatticeUplifter()
        >>> lu.inner_base_row_nodes
        array([1, 3, 4])

        >>> hg = HexModelGrid(5, 6, 1.0, orientation='vertical', shape='rect', reorient_links=True)
        >>> lu = LatticeUplifter(grid=hg)
        >>> lu.inner_base_row_nodes
        array([1, 2, 3, 4])
        """
        # Do the base class init
        super(LatticeUplifter, self).__init__(grid, node_state, propid,
                                              prop_data, prop_reset_value)

        # Remember the IDs of nodes on the bottom row
        self.inner_base_row_nodes = zeros(self.nc - 2, dtype=int)
        n_in_lower = (self.nc // 2) - 1  # num inner nodes bottom row
        upper_start = (self.nc + 1) // 2
        n_in_upper = upper_start - 1
        self.inner_base_row_nodes[:n_in_lower] = arange(1, n_in_lower + 1)
        self.inner_base_row_nodes[n_in_lower:] = arange(upper_start,
                                                        upper_start + \
                                                        n_in_upper)

        if self.propid is not None:
            self.inner_top_row_nodes = self.inner_base_row_nodes + \
                                       ((self.nr - 1) * self.nc)


    def uplift_interior_nodes(self, rock_state=1):
        """
        Simulate 'vertical' displacement by shifting contents of node_state

        Examples
        --------
        >>> lu = LatticeUplifter()
        >>> lu.node_state[:] = arange(len(lu.node_state))
        >>> lu.uplift_interior_nodes(rock_state=25)
        >>> lu.node_state # doctest: +NORMALIZE_WHITESPACE
        array([ 0, 25,  2, 25, 25,
                5,  1,  7,  3,  4,
               10,  6, 12,  8,  9,
               15, 11, 17, 13, 14,
               20, 16, 22, 18, 19])
        """

        # Shift the node states up by a full row. A "full row" includes two
        # staggered rows.
        for r in range(self.nr - 1, 0, -1):
            # This row gets the contents of the nodes 1 row down
            self.node_state[self.inner_base_row_nodes+self.nc*r] = \
                    self.node_state[self.inner_base_row_nodes+self.nc*(r-1)]

        # Fill the bottom rows with "fresh material" (code = rock_state)
        self.node_state[self.inner_base_row_nodes] = rock_state

        # Shift the node states up by two rows: two because the grid is
        # staggered, and we don't want any horizontal offset.

        # STILL TO DO: MAKE SURE THIS HANDLES WRAP PROPERLY (I DON'T THINK
        # IT DOES NOW)
        # If propid (property ID or index) is defined, shift that too.
        if self.propid is not None:
            top_row_propid = self.propid[self.inner_top_row_nodes]
            for r in range(self.nr-1, 1, -1):
                self.propid[self.inner_base_row_nodes+self.nc*r] =  \
                            self.propid[self.inner_base_row_nodes+self.nc*(r-2)]
            self.propid[self.inner_base_row_nodes] = top_row_propid
            self.prop_data[self.propid[self.inner_base_row_nodes]] = self.prop_reset_value


if __name__=='__main__':
    import doctest
    doctest.testmod()
