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
from numpy import amax, zeros, arange, array
from pylab import figure, show

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
        Create and initialize a HexLatticeTectonicizer

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
            self.node_state = self.grid.add_zeros('node', 'node_state_map', \
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
    >>> lnf = LatticeNormalFault(0.01, grid, ns, pid, pdata, 0.0)
    >>> lnf.first_fw_col
    1
    >>> lnf.num_fw_rows
    array([0, 1, 3, 4, 5])
    >>> lnf.incoming_node
    array([ 5, 10, 11, 15])
    >>> lnf.outgoing_node
    array([22, 23, 24, 18])
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
        >>> lnf = LatticeNormalFault(0.01, grid, ns, pid, pdata, 0.0)
        >>> lnf.first_fw_col
        1
        >>> lnf.num_fw_rows
        array([0, 1, 3, 4, 5])
        >>> lnf.incoming_node
        array([ 5, 10, 11, 15])
        >>> lnf.outgoing_node
        array([22, 23, 24, 18])
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
        assert fault_x_intercept< amax(self.grid.node_x) - 0.57735, 'err'

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
            while (current_row<self.nr and
                   in_footwall[bottom_row_node_id[c] + self.nc*current_row]):
                self.num_fw_rows[c] += 1
                current_row += 1

        # If we're handling properties and property IDs, we need to do some setup
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
            for c in range(self.first_fw_col, self.nc - 1):

                # This little loop appends to the incoming node list:
                # either just the bottom node in this column (ID=self.nr*c)
                # or that plus the node above it (ID=self.nr*c+1). Note that
                # 2-c%s evaluates to 2 for even-numbered columns and 1 for
                # odd-numbered columns. So the loop is either 1 or 2 iterations.
                for node_id in range(self.nr * c, self.nr * c + (2 - c % 2)):
                    incoming_node_list.append(node_id)

            # Convert to a numpy array that belongs to the class (so we can
            # use it in the do_offset() method).
            self.incoming_node = array(incoming_node_list, dtype=int)

            # Next, find the IDs the "outgoing" nodes. There will always be
            # some of these along the right side of the grid. Depending on the
            # height of the grid and the position of the fault, there may or
            # may not also be some along the top.
            #
            # Which of the nodes on the right side outgoing? The lowermost one
            # won't be. If the right-most column is even-numbered, the next one
            # up won't be either. So call the ID of the first potential
            # outgoing node onthe right side the "base_id". Then we also have
            # the "top_id", which is the ID of *either* the top-most node of
            # the column (ID=# of nodes in grid-1), *or* the top footwall node.
            #
            # The next line translates as: Take the ID of the bottom node in
            # the right-most column (nr x (nc-1)) and add to it either 1
            # (if the last column is odd-numbered) or 2 (if even-numbered).
            # Example: in a 5x5 grid, it is node ID 22 (two up from the base)
            base_id = self.nr*(self.nc-1)+(1+self.nc%2)
            #
            # The next line finds the ID of the top outgoing node on the right
            # side. Here's how it works. Find the ID of the topmost node in
            # the column: that's the # of nodes in the grid minus one.
            # Example: in a 5x5 grid, that's node number 24.
            #
            # Then find the ID of the topmost node IN THE FOOTWALL. To get
            # this, we start with the ID at the top of the column to the left,
            # and we start we start number of footwall nodes in the right-most
            # column. Example: in a 5x5 grid, the top of the next-to-rightmost
            # column is node 19. If the fault position is x=0.0, there will be
            # 5 footwall nodes here. 19+5 = 24 (meaning the top of the column
            # and the top of the footwall happen to the same)
            top_id = min(self.nr * self.nc - 1,
                         (self.nr * (self.nc - 1) - 1) + \
                         self.num_fw_rows[self.nc - 1])

            # Having found the top and the base, we now append all these to
            # the list of outgoing nodes.
            outgoing_node_list = []
            for node_id in range(base_id, top_id+1):
                outgoing_node_list.append(node_id)

            # Now, at this point, we might have found all the outgoing nodes.
            # Or there might be some more along the top edge of the grid. The
            # latter will be true if the number of outgoing nodes found so far
            # is less than the number there should be (same as # of incoming).
            #
            # start with the next-to-rightmost column; if we need to, we'll
            # work right-to-left from column to column
            col = self.nc - 2
            # If the following is true, we still have more nodes to add along
            # the top
            while len(outgoing_node_list) < len(self.incoming_node):

                # Add the top-most footwall node in this column
                outgoing_node_list.append(self.nr*col+self.num_fw_rows[col]-1)

                # If we're on an odd-numbered column, and not already full,
                # add the next one down too
                if col%2==1 and len(outgoing_node_list)<len(self.incoming_node):
                    outgoing_node_list.append(self.nr*col+self.num_fw_rows[col]-2)

            # Finally, convert the outgoing node list to an array stored in this
            # object
            self.outgoing_node = array(outgoing_node_list, dtype=int)

    def do_offset(self, rock_state=1):
        """Applies 60-degree normal-fault offset to a hexagonal grid with
        vertical node orientation and rectangular arrangement of nodes.

        Parameters
        ----------
        rock_state : int
            State code to apply to new cells introduced along bottom row.

        Examples
        --------
#        >>> import numpy as np
#        >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
#        >>> from landlab import HexModelGrid
#        >>> pid = np.arange(25, dtype=int)
#        >>> pdata = np.arange(25)
#        >>> ns = np.arange(25, dtype=int)
#        >>> grid = HexModelGrid(5, 5, 1.0, orientation='vertical', shape='rect', reorient_links=True)
#        >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
#        >>> lnf.do_offset()
#        >>> lnf.propid
#        array([ 0,  1,  2,  3,  4, 22,  6,  7,  8,  9, 23, 24,  5, 13, 14, 18, 10,
#               11, 12, 19, 20, 21, 15, 16, 17])
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
        for c in range(self.grid.number_of_node_columns-1, self.first_fw_col-1, -1):

            # Odd-numbered rows are shifted up in the hexagonal, vertically
            # oriented lattice
            row_offset = 2 - (c % 2)

            # Number of base nodes in the footwall in this column. Either 1 or 2.
            n_base_nodes = min(self.num_fw_rows[c], row_offset)

            # ID of the bottom footwall node
            bottom_node = c*self.grid.number_of_node_rows

            # The bottom 1 or 2 nodes in this column are set to rock
            self.node_state[bottom_node:(bottom_node+n_base_nodes)] = rock_state

            # "indices" here contains the array indices of those nodes in this
            # column that are to be replaced by the ones in the column to the
            # left and down one or two nodes. We do this replacement if indices
            # contains any data.
            indices = arange(n_base_nodes, self.num_fw_rows[c], dtype=int)+bottom_node
            if len(indices)>0:
                self.node_state[indices] = self.node_state[indices-(self.nr+row_offset)]
                if self.propid is not None:
                    #print 'in col',c,'replacing nodes',indices,'with',indices-(self.nr+row_offset)
                    #if c==2:
                        #print 'node 18 propid changing from',self.propid[18],'(',self.prop_data[self.propid[18]],')'
                    self.propid[indices] = self.propid[indices-(self.nr+row_offset)]
                    #if c==2:
                        #print '    to',self.propid[18],'(',self.prop_data[self.propid[18]],')'

        if self.propid is not None:
            self.propid[self.incoming_node] = propids_for_incoming_nodes
            self.prop_data[self.propid[self.incoming_node]] = self.prop_reset_value
            #print 'pid after:',self.propid
            #print 'propdata after:',self.prop_data

        if self.first_fw_col==0:
            self.node_state[:self.n_footwall_rows[0]] = rock_state


class LatticeUplifter(HexLatticeTectonicizer):
    """Handles vertical uplift of interior (not edges) for a hexagonal lattice
    with vertical node orientation and rectangular node arrangement.

    Examples
    --------
#    >>> lu = LatticeUplifter()
#    >>> lu.base_row_nodes
#    array([0, 1, 2, 3, 4])
    """
    def __init__(self, grid=None, node_state=None, propid=None, prop_data=None, prop_reset_value=None):
        """
        Create and initialize a LatticeUplifter

        Examples
        --------
        >>> lu = LatticeUplifter()
        >>> lu.inner_base_row_nodes
        array([1, 3, 4])
        """
        # Do the base class init
        super(LatticeUplifter, self).__init__(grid, node_state, propid, prop_data, prop_reset_value)

        # Remember the IDs of nodes on the bottom row
        if self.nc % 2 == 0: # if even num cols
            self.inner_base_row_nodes = arange(1, self.nc - 1, dtype=int)
        else: # if odd num cols
            self.inner_base_row_nodes = zeros(self.nc - 2, dtype=int)
            n_in_row1 = self.nc // 2    # num inner nodes 2nd row from bottom
            n_in_row0 = n_in_row1 - 1  # num inner nodes bottom row
            self.inner_base_row_nodes[:n_in_row0] = arange(1, n_in_row0 + 1)
            self.inner_base_row_nodes[n_in_row0:] = arange(n_in_row0 + 2, 
                                                           self.nc)
        #print 'LU INIT HERE******************'
        if self.propid is not None:
            self.inner_top_row_nodes = self.inner_base_row_nodes+(self.nr-1)*self.nc
            #print 'top:',self.top_row_nodes
            #print 'base:',self.base_row_nodes

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
        #print 'in uin, ns is'
        #print self.node_state
        #print 'and propid before is'
        #print self.propid
        #print 'and prop before is'
        #print self.prop_data[self.propid]

        # Shift the node states up by a full row. A "full row" includes two
        # staggered rows.
        for r in range(self.nr-1, 0, -1):
            # This row gets the contents of the nodes 1 row down
            self.node_state[self.inner_base_row_nodes+self.nc*r] = \
                    self.node_state[self.inner_base_row_nodes+self.nc*(r-1)]

        # Fill the bottom rows with "fresh material" (code = rock_state)
        self.node_state[self.inner_base_row_nodes] = rock_state

        # Shift the node states up by two rows: two because the grid is
        # staggered, and we don't want any horizontal offset.
#        for r in range(self.nr-1, 1, -1):
#            # This row gets the contents of the nodes 2 rows down
#            self.node_state[self.base_row_nodes+self.nc*r] = \
#                    self.node_state[self.base_row_nodes+self.nc*(r-2)]
#
#        # Fill the bottom two rows with "fresh material" (code = rock_state)
#        self.node_state[self.base_row_nodes] = rock_state
#        self.node_state[self.base_row_nodes+self.nc] = rock_state

        # STILL TO DO: MAKE SURE THIS HANDLES WRAP PROPERLY (I DON'T THINK
        # IT DOES NOW)
        # If propid (property ID or index) is defined, shift that too.
        #print 'uid1: propid 11-12=',self.propid[11:13],self.prop_data[self.propid[11:13]]
        if self.propid is not None:
            top_row_propid = self.propid[self.inner_top_row_nodes]
            for r in range(self.nr-1, 1, -1):
                self.propid[self.inner_base_row_nodes+self.nc*r] =  \
                            self.propid[self.inner_base_row_nodes+self.nc*(r-2)]
            #print 'uid2: propid 11-12=',self.propid[11:13],self.prop_data[self.propid[11:13]]
            self.propid[self.inner_base_row_nodes] = top_row_propid
            #print 'uid3: propid 11-12=',self.propid[11:13],self.prop_data[self.propid[11:13]]
            self.prop_data[self.propid[self.inner_base_row_nodes]] = self.prop_reset_value
            #print 'uid4: propid 11-12=',self.propid[11:13],self.prop_data[self.propid[11:13]]
            #print 'in UIN, pid is'
            #print self.propid
            #print self.prop_data[self.propid]


#def test_create_lnf(nr, nc):
#
#    pid = arange(nr*nc, dtype=int)
#    ns = arange(nr*nc, dtype=int)
#    pdata = arange(nr*nc)
#    grid = HexModelGrid(nr, nc, 1.0, orientation='vertical', shape='rect', reorient_links=True)
#    lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
#    #for i in range(grid.number_of_nodes):
#    #    print i, grid.node_x[i], grid.node_y[i]
#    return lnf


def main():
    """The main function is used just to do some unit tests.

    Examples
    --------
#    >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import test_create_lnf
#    >>> lnf = test_create_lnf(4, 4)
#    >>> lnf.incoming_node
#    array([4, 8, 9])
#    >>> lnf.outgoing_node
#    array([13, 14, 15])
#    >>> lnf.do_offset()
#    >>> lnf.propid
#    array([ 0,  1,  2,  3, 13,  5,  6,  7, 14, 15,  4, 11, 12,  8,  9, 10])
#    >>> lnf = test_create_lnf(4, 5)
#    >>> lnf.incoming_node
#    array([ 4,  8,  9, 12])
#    >>> lnf.outgoing_node
#    array([18, 19, 15, 14])
#    >>> lnf.do_offset()
#    >>> lnf.propid
#    array([ 0,  1,  2,  3, 18,  5,  6,  7, 19, 15,  4, 11, 14,  8,  9, 10, 16,
#           17, 12, 13])
    """
    pid = arange(_DEFAULT_NUM_ROWS*_DEFAULT_NUM_COLS, dtype=int)
    pdata = arange(_DEFAULT_NUM_ROWS*_DEFAULT_NUM_COLS)
    lnf = LatticeNormalFault(propid=pid, prop_data=pdata, prop_reset_value=0.0)

    for i in range(3):
        lnf.do_offset()
    lnf.grid.hexplot(lnf.node_state)
    show()

    lu = LatticeUplifter()
    lu.uplift_interior_nodes()
    figure(2)
    for i in range(2):
        lu.uplift_interior_nodes()
        lu.grid.hexplot(lu.node_state)

if __name__=='__main__':
    #main()
    import doctest
    doctest.testmod()
