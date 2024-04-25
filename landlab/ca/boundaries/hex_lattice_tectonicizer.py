"""hex_lattice_tectonicizer.py.

Models discrete normal-fault offset on a 2D hex lattice with a rectangular
node layout and with one orientation of the nodes being vertical.

The intention here is to use a particle (LCA) model to represent the evolution
of a 2D hillslope, with the hex_lattice_tectonicizer serving to shift the nodes
either upward (simple vertical uplift relative to baselevel), or up and
sideways (representing motion on a fault plane).

Created on Mon Nov 17 08:01:49 2014

@author: gtucker
"""

from numpy import amax
from numpy import arange
from numpy import array
from numpy import cos
from numpy import logical_and
from numpy import logical_or
from numpy import logical_xor
from numpy import pi
from numpy import sqrt
from numpy import tan
from numpy import where
from numpy import zeros

from landlab import HexModelGrid
from landlab import LinkStatus
from landlab.core.utils import as_id_array

from ..cfuncs import get_next_event_new  # , update_link_state_new

_DEFAULT_NUM_ROWS = 5
_DEFAULT_NUM_COLS = 5
_TAN60 = 1.732
_NEVER = 1.0e50  # this arbitrarily large val is also defined in ..cfuncs.pyx


def is_interior_link(link, grid):
    """Return True if both nodes are core; False otherwise."""
    return (
        grid.status_at_node[grid.node_at_link_tail[link]] == grid.BC_NODE_IS_CORE
        and grid.status_at_node[grid.node_at_link_head[link]] == grid.BC_NODE_IS_CORE
    )


def is_perim_link(link, grid):
    """Return True if both nodes are boundaries; False otherwise.

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> import numpy as np
    >>> mg = HexModelGrid(
    ...     (3, 4), spacing=1.0, orientation="vertical", node_layout="rect"
    ... )
    >>> is_perim_link(1, mg)
    True
    >>> is_perim_link(11, mg)
    False
    """
    return (
        grid.status_at_node[grid.node_at_link_tail[link]] != grid.BC_NODE_IS_CORE
        and grid.status_at_node[grid.node_at_link_head[link]] != grid.BC_NODE_IS_CORE
    )


class HexLatticeTectonicizer:
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

    def __init__(
        self,
        grid=None,
        node_state=None,
        propid=None,
        prop_data=None,
        prop_reset_value=None,
    ):
        """Create and initialize a HexLatticeTectonicizer.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid((6, 6), node_layout="rect")
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
            self.grid = HexModelGrid(
                (num_rows, num_cols),
                spacing=1.0,
                orientation="vertical",
                node_layout="rect",
                reorient_links=True,
            )
        else:
            # Make sure caller passed the right type of grid
            assert grid.orientation == "vertical", "Grid must have vertical orientation"

            # Keep a reference to the grid
            self.grid = grid

        # If needed, create node-state grid
        if node_state is None:
            self.node_state = self.grid.add_zeros("node", "node_state_map", dtype=int)
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
    offsets for a hex grid with vertical columns and rectangular node_layout.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import HexModelGrid
    >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
    >>> pid = np.arange(25, dtype=int)
    >>> pdata = np.arange(25)
    >>> ns = np.arange(25, dtype=int)
    >>> grid = HexModelGrid(
    ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
    ... )
    >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
    >>> lnf.first_fw_col
    1
    >>> lnf.num_fw_rows
    array([0, 1, 3, 4, 5])
    >>> lnf.incoming_node
    array([1, 3, 4, 6])
    >>> lnf.outgoing_node
    array([12, 17, 19, 22])

    >>> pid = np.arange(16, dtype=int)
    >>> ns = np.arange(16, dtype=int)
    >>> pdata = np.arange(16)
    >>> grid = HexModelGrid(
    ...     (4, 4), spacing=1.0, orientation="vertical", node_layout="rect"
    ... )
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

    >>> pid = np.arange(20, dtype=int)
    >>> ns = np.arange(20, dtype=int)
    >>> pdata = np.arange(20)
    >>> grid = HexModelGrid(
    ...     (4, 5), spacing=1.0, orientation="vertical", node_layout="rect"
    ... )
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

    def __init__(
        self,
        fault_x_intercept=0.0,
        grid=None,
        node_state=None,
        propid=None,
        prop_data=None,
        prop_reset_value=None,
    ):
        """Create and initialize a LatticeNormalFault object.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import (
        ...     LatticeNormalFault,
        ... )

        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ns = np.arange(25, dtype=int)
        >>> grid = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
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
        >>> grid = HexModelGrid(
        ...     (4, 4), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
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
        >>> grid = HexModelGrid(
        ...     (5, 9), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
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
        super().__init__(grid, node_state, propid, prop_data, prop_reset_value)
        # Set up data structures:
        #   Make sure the footwall location is such that the fault actually
        #   cuts across the grid. This means the x intercept has to be, at
        #   the very least, no smaller than the biggest x-coordinate, and if
        #   there is an even number of columns, it must be smaller than that
        #   number minus 1/tangent 60 degrees (0.57735)
        assert fault_x_intercept < (amax(self.grid.node_x) - 0.57735), "err"

        #   Figure out which nodes are and are not within the footwall
        in_footwall = self.grid.node_y < _TAN60 * (self.grid.node_x - fault_x_intercept)

        # Set up array of link offsets: when slip occurs, what link's data get
        # copied into the present link? Which links get cut by fault plane and
        # need to have their states reset?
        self._setup_link_offsets(in_footwall)
        self._setup_links_to_update_after_offset(in_footwall)

        # Helpful to have an array of node IDs for the bottom full row. Because
        # the nodes in the bottom row are staggered in a vertical, rectangular
        # hex grid, the IDs go: 0, M, 1, M+1, 2, M+2, ... etc., where M is half
        # the number of columns, rounded up (so, for example, 3 for a 5- or 6-
        # column grid, etc.)
        half_num_cols = (self.nc + 1) // 2
        bottom_row_node_id = (
            arange(self.nc) // 2 + (arange(self.nc) % 2) * half_num_cols
        )

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
            assert n < self.nc, "overflow in loop"

        #   Remember the number of footwall rows in each column
        self.num_fw_rows = zeros(self.nc, dtype=int)
        for c in range(self.nc):
            current_row = 0
            while (
                current_row < self.nr
                and in_footwall[bottom_row_node_id[c] + self.nc * current_row]
            ):
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
            lower_right_node = (self.nc % 2) * (self.nc // 2) + ((self.nc + 1) % 2) * (
                self.nc - 1
            )
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
                if (
                    (self.grid.node_x[n] > x_threshold and self.grid.node_y[n] > 1.25)
                    or self.grid.node_y[n] > y_threshold
                ) and in_footwall[n]:
                    outgoing_node_list.append(n)

            # Finally, convert the outgoing node list to an array stored in
            # this object
            self.outgoing_node = array(outgoing_node_list, dtype=int)

    def _link_in_footwall(self, link, node_in_footwall):
        """Return True of both nodes are in footwall, False otherwise."""
        return (
            node_in_footwall[self.grid.node_at_link_tail[link]]
            and node_in_footwall[self.grid.node_at_link_head[link]]
        )

    def _get_link_orientation(self, link):
        """Return link orientation code for given link."""
        assert self.grid.orientation[0] == "v", "assumes vertical orientation"
        head = self.grid.node_at_link_head[link]
        tail = self.grid.node_at_link_tail[link]
        dx = self.grid.x_of_node[head] - self.grid.x_of_node[tail]
        dy = self.grid.y_of_node[head] - self.grid.y_of_node[tail]
        if dy > dx:
            return 0  # vertical
        elif dy > 0:
            return 1  # right and up
        else:
            return 2  # right and down

    def _setup_link_offsets(self, node_in_footwall):
        """Set up array with link IDs for shifting link data up and right.

        Notes
        -----
        The array contains the ID of the link TO WHICH the contents get
        shifted upon fault slip.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import (
        ...     LatticeNormalFault,
        ... )
        >>> from landlab import HexModelGrid
        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ns = np.arange(25, dtype=int)
        >>> grid = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> lnf = LatticeNormalFault(-0.01, grid, ns, pid, pdata, 0.0)
        >>> lnf.link_offset_id[14:22]
        array([35, 15, 16, 17, 38, 19, 20, 41])
        >>> lnf.first_link_shifted_from
        14
        >>> lnf.first_link_shifted_to
        35

        >>> pid = np.arange(36, dtype=int)
        >>> pdata = np.arange(36)
        >>> ns = np.arange(36, dtype=int)
        >>> grid = HexModelGrid(
        ...     (6, 6), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> lnf = LatticeNormalFault(-0.1, grid, ns, pid, pdata, 0.0)
        >>> lnf.first_link_shifted_from
        17
        >>> lnf.first_link_shifted_to
        42
        >>> lnf.link_offset_id[17:39]
        array([42, 43, 19, 20, 21, 46, 23, 24, 50, 51, 27, 28, 29, 55, 31, 32, 33,
               59, 35, 36, 37, 62])
        """
        self.link_offset_id = arange(self.grid.number_of_links, dtype=int)
        nc = self.grid.number_of_node_columns
        default_offset = 2 * nc + 2 * (nc - 1) + nc // 2
        self.first_link_shifted_from = 0
        self.first_link_shifted_to = 0

        for ln in range(self.grid.number_of_links - (default_offset + 1)):
            if self._link_in_footwall(ln, node_in_footwall) and is_interior_link(
                ln, self.grid
            ):
                tail_node = self.grid.node_at_link_tail[ln]
                (_, c) = self.grid.node_row_and_column(tail_node)
                link_orientation = self._get_link_orientation(ln)
                offset = default_offset
                if nc % 2 == 1:  # odd number of columns
                    if (link_orientation + ((c - 1) % 2)) == 2:
                        offset += 1
                else:  # even number of columns
                    if (c % 2) == 0 and link_orientation == 0:
                        offset -= 1
                if is_interior_link(ln + offset, self.grid):
                    self.link_offset_id[ln] = ln + offset
                    if self.first_link_shifted_from == 0:
                        self.first_link_shifted_from = ln
                        self.first_link_shifted_to = ln + offset

    def _setup_links_to_update_after_offset(self, in_footwall):
        """Create and store array with IDs of links for which to update
        transitions after fault offset.

        These are: all active boundary links with at least one node in the
        footwall, plus the lowest non-boundary links, including the
        next-to-lowest vertical links and those angling that are below them,
        plus the fault-crossing links.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid((5, 5), orientation="vertical", node_layout="rect")
        >>> lu = LatticeNormalFault(fault_x_intercept=-0.01, grid=hg)
        >>> lu.first_link_shifted_to
        35
        >>> lu.links_to_update
        array([ 5,  8,  9, 11, 12, 13, 14, 15, 16, 18, 20, 21, 22, 23, 24, 25, 27,
               28, 29, 31, 34, 36, 40, 42, 44, 48, 49, 51])
        """
        g = self.grid
        lower_active = logical_and(
            arange(g.number_of_links) < self.first_link_shifted_to,
            g.status_at_link == LinkStatus.ACTIVE,
        )
        link_in_fw = logical_or(
            in_footwall[g.node_at_link_tail], in_footwall[g.node_at_link_head]
        )
        lower_active_fw = logical_and(lower_active, link_in_fw)
        active_bnd = logical_and(
            g.status_at_link == LinkStatus.ACTIVE,
            logical_or(
                g.status_at_node[g.node_at_link_tail] != 0,
                g.status_at_node[g.node_at_link_head] != 0,
            ),
        )
        active_bnd_fw = logical_and(active_bnd, link_in_fw)
        crosses_fw = logical_and(
            g.status_at_link == LinkStatus.ACTIVE,
            logical_xor(
                in_footwall[g.node_at_link_tail], in_footwall[g.node_at_link_head]
            ),
        )
        update = logical_or(logical_or(lower_active_fw, active_bnd_fw), crosses_fw)
        self.links_to_update = as_id_array(where(update)[0])

    def assign_new_link_state_and_transition(self, link, ca, current_time):
        """Update state and schedule new transition for given link."""
        tail_state = ca.node_state[self.grid.node_at_link_tail[link]]
        head_state = ca.node_state[self.grid.node_at_link_head[link]]
        orientation = ca.link_orientation[link]
        new_link_state = int(
            orientation * ca.num_node_states_sq
            + tail_state * ca.num_node_states
            + head_state
        )
        ca.update_link_state_new(link, new_link_state, current_time)

    def shift_scheduled_transitions(self, ca, current_time):
        """Update link IDs in scheduled events at offset links.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.oriented_hex_cts import OrientedHexCTS
        >>> from landlab.ca.celllab_cts import Transition
        >>> import numpy as np

        >>> mg = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> nsd = {0: "yes", 1: "no"}
        >>> xnlist = []
        >>> xnlist.append(Transition((0, 0, 0), (1, 1, 0), 1.0, "test"))
        >>> xnlist.append(Transition((0, 0, 1), (1, 1, 1), 1.0, "test"))
        >>> xnlist.append(Transition((0, 0, 2), (1, 1, 2), 1.0, "test"))
        >>> nsg = mg.add_zeros("node", "node_state_grid")
        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ohcts = OrientedHexCTS(mg, nsd, xnlist, nsg)
        >>> lnf = LatticeNormalFault(-0.1, grid=mg)
        >>> pq = ohcts.priority_queue._queue
        >>> (int(1000 * pq[11][0]), pq[11][1:])
        (752, (11, 21))
        >>> (int(1000 * pq[12][0]), pq[12][1:])
        (483, (9, 18))
        >>> (int(1000 * pq[30][0]), pq[30][1:])
        (575, (6, 14))
        >>> lnf.do_offset(ca=ohcts)
        >>> (int(1000 * pq[48][0]), pq[48][1:])
        (752, (11, 41))
        >>> (int(1000 * pq[54][0]), pq[54][1:])
        (483, (9, 38))
        >>> (int(1000 * pq[61][0]), pq[61][1:])
        (575, (6, 35))
        """
        for i in range(len(ca.priority_queue._queue)):
            link = ca.priority_queue._queue[i][2]
            if self.link_offset_id[link] != link:
                ca.priority_queue._queue[i] = (
                    ca.priority_queue._queue[i][0],
                    ca.priority_queue._queue[i][1],
                    self.link_offset_id[link],
                )

    def shift_link_states(self, ca, current_time):
        """Shift link data up and right.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.oriented_hex_cts import OrientedHexCTS
        >>> from landlab.ca.celllab_cts import Transition
        >>> import numpy as np

        >>> mg = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> nsd = {0: "yes", 1: "no"}
        >>> xnlist = []
        >>> xnlist.append(Transition((1, 0, 0), (1, 1, 0), 1.0, "frogging"))
        >>> xnlist.append(Transition((1, 0, 1), (1, 1, 1), 1.0, "frogging"))
        >>> xnlist.append(Transition((1, 0, 2), (1, 1, 2), 1.0, "frogging"))
        >>> nsg = mg.add_zeros("node", "node_state_grid")
        >>> nsg[:10] = 1
        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ohcts = OrientedHexCTS(mg, nsd, xnlist, nsg)
        >>> ohcts.link_state
        array([ 0,  0,  0,  0,  0,  3,  0,  0,  7, 11,  0,  3,  3,  7, 11,  7, 11,
                0,  2,  0,  9,  6,  9,  6,  2,  2,  4,  8,  4,  8,  0,  0,  0,  8,
                4,  8,  4,  0,  0,  4,  8,  4,  8,  0,  0,  0,  8,  4,  8,  4,  0,
                0,  0,  0,  0,  0])
        >>> lnf = LatticeNormalFault(-0.01, mg, nsg, pid, pdata, 0.0)
        >>> lnf.first_link_shifted_from
        14
        >>> lnf.link_offset_id
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 35, 15, 16,
               17, 38, 19, 20, 41, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
               34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
               51, 52, 53, 54, 55])
        >>> lnf.shift_link_states(ohcts, 0.0)
        >>> ohcts.link_state
        array([ 0,  0,  0,  0,  0,  3,  0,  0,  7, 11,  0,  3,  3,  7, 11,  7, 11,
                0,  2,  0,  9,  6,  9,  6,  2,  2,  4,  8,  4,  8,  0,  0,  0,  8,
                4, 11,  4,  0,  2,  4,  8,  6,  8,  0,  0,  0,  8,  4,  8,  4,  0,
                0,  0,  0,  0,  0])
        """
        num_links = self.grid.number_of_links
        for lnk in range(num_links - 1, self.first_link_shifted_from - 1, -1):
            link_offset = self.link_offset_id[lnk]
            if link_offset != lnk:
                ca.link_state[link_offset] = ca.link_state[lnk]
                ca.next_trn_id[link_offset] = ca.next_trn_id[lnk]
                ca.next_update[link_offset] = ca.next_update[lnk]

        self.shift_scheduled_transitions(ca, current_time)

        for lnk in self.links_to_update:
            self.assign_new_link_state_and_transition(lnk, ca, current_time)

    def do_offset(self, ca=None, current_time=0.0, rock_state=1):
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
        >>> from landlab.ca.boundaries.hex_lattice_tectonicizer import (
        ...     LatticeNormalFault,
        ... )
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.oriented_hex_cts import OrientedHexCTS
        >>> from landlab.ca.celllab_cts import Transition

        >>> pid = np.arange(25, dtype=int)
        >>> pdata = np.arange(25)
        >>> ns = np.arange(25, dtype=int)
        >>> grid = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
        >>> lnf.do_offset(rock_state=25)
        >>> ns
        array([ 0, 25, 25, 25, 25,  5, 25, 25,  8,  1, 10,  3,  4, 13,  6, 15, 16,
                9, 18, 11, 20, 21, 14, 23, 24])
        >>> lnf.propid
        array([ 0, 12,  2, 17, 19,  5, 22,  7,  8,  1, 10,  3,  4, 13,  6, 15, 16,
                9, 18, 11, 20, 21, 14, 23, 24])

        >>> ns[5:] = 0
        >>> ns[:5] = 1
        >>> nsd = {0: "yes", 1: "no"}
        >>> xnlist = []
        >>> xnlist.append(Transition((0, 0, 0), (1, 1, 0), 1.0, "test"))
        >>> ohcts = OrientedHexCTS(grid, nsd, xnlist, ns)
        >>> lnf = LatticeNormalFault(0.0, grid, ns, pid, pdata, 0.0)
        >>> ohcts.link_state[5:33]
        array([2, 0, 0, 6, 9, 0, 2, 2, 4, 8, 4, 8, 0, 0, 0, 8, 4, 8, 4, 0, 0, 4,
               8, 4, 8, 0, 0, 0])
        >>> lnf.do_offset(ca=ohcts, current_time=0.0, rock_state=1)
        >>> ohcts.link_state[5:33]
        array([ 3,  0,  0,  7, 11,  0,  2,  3,  4,  9,  7, 11,  0,  3,  0,  8,  5,
               11,  7,  0,  2,  4,  9,  6,  9,  0,  2,  0])
        """

        # If we need to shift the property ID numbers, we'll first need to
        # record the property IDs in those nodes that are about to "shift off
        # the grid" (or rather, their contents will shift) due to tectonic
        # motion. We'll call these "nodes to replace".
        if self.propid is not None:
            # Now, remember the property IDs of the nodes along the right side
            # and possibly top that are about to shift off the grid
            propids_for_incoming_nodes = self.propid[self.outgoing_node]

        # We go column-by-column, starting from the right side
        for c in range(self.grid.number_of_node_columns - 1, self.first_fw_col - 1, -1):
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
                self.node_state[bottom_node + self.nc] = rock_state

            # "indices" here contains the array indices of those nodes in this
            # column that are to be replaced by the ones in the column to the
            # left and down one or two nodes. We do this replacement if indices
            # contains any data.
            first_repl = bottom_node + n_base_nodes * self.nc
            last_repl = (
                first_repl + (self.num_fw_rows[c] - (n_base_nodes + 1)) * self.nc
            )
            indices = arange(first_repl, last_repl + 1, self.nc)

            if len(indices) > 0:
                offset = (
                    self.nc + ((self.nc + 1) // 2) + ((c + 1) % 2) * ((self.nc + 1) % 2)
                )
                self.node_state[indices] = self.node_state[indices - offset]
                if self.propid is not None:
                    self.propid[indices] = self.propid[indices - offset]

        if self.propid is not None:
            self.propid[self.incoming_node] = propids_for_incoming_nodes
            self.prop_data[self.propid[self.incoming_node]] = self.prop_reset_value

        if ca is not None:
            self.shift_link_states(ca, current_time)


class LatticeUplifter(HexLatticeTectonicizer):
    """Handles vertical uplift of interior (not edges) for a hexagonal lattice
    with vertical node orientation and rectangular node arrangement."""

    def __init__(
        self,
        grid=None,
        node_state=None,
        propid=None,
        prop_data=None,
        prop_reset_value=None,
        opt_block_layer=False,
        block_ID=9,
        block_layer_dip_angle=0.0,
        block_layer_thickness=1.0,
        layer_left_x=0.0,
        y0_top=0.0,
    ):
        """Create and initialize a LatticeUplifter.

        Examples
        --------
        >>> lu = LatticeUplifter()
        >>> lu.inner_base_row_nodes
        array([1, 3, 4])

        >>> hg = HexModelGrid(
        ...     (5, 6), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> lu = LatticeUplifter(grid=hg)
        >>> lu.inner_base_row_nodes
        array([1, 2, 3, 4])
        """
        # Do the base class init
        super().__init__(grid, node_state, propid, prop_data, prop_reset_value)

        # Remember the IDs of nodes on the bottom row
        self.inner_base_row_nodes = zeros(self.nc - 2, dtype=int)
        n_in_lower = (self.nc // 2) - 1  # num inner nodes bottom row
        upper_start = (self.nc + 1) // 2
        n_in_upper = upper_start - 1
        self.inner_base_row_nodes[:n_in_lower] = arange(1, n_in_lower + 1)
        self.inner_base_row_nodes[n_in_lower:] = arange(
            upper_start, upper_start + n_in_upper
        )

        if self.propid is not None:
            self.inner_top_row_nodes = self.inner_base_row_nodes + (
                (self.nr - 1) * self.nc
            )

        self._setup_links_to_update_after_uplift()

        # Handle option for a layer of "blocks"
        self.opt_block_layer = opt_block_layer
        if opt_block_layer:
            self.cum_uplift = 0.0
            self.block_ID = block_ID
            self.block_layer_thickness = block_layer_thickness
            self.block_layer_dip_angle = block_layer_dip_angle
            self.layer_left_x = layer_left_x
            self.y0_top = y0_top

    def _setup_links_to_update_after_uplift(self):
        """Create and store array with IDs of links for which to update
        transitions after uplift.

        These are: all active boundary links, plus the lowest non-boundary
        links, including the next-to-lowest vertical links and those angling
        that are below them.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid((6, 6), orientation="vertical", node_layout="rect")
        >>> lu = LatticeUplifter(grid=hg)
        >>> lu.links_to_update
        array([ 6,  7,  9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 22, 23, 24, 28,
               32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 70, 71, 72, 73, 74, 75, 77,
               78])
        >>> hg = HexModelGrid((5, 5), orientation="vertical", node_layout="rect")
        >>> lu = LatticeUplifter(grid=hg)
        >>> lu.links_to_update
        array([ 5,  8,  9, 11, 12, 13, 14, 15, 16, 18, 20, 23, 26, 29, 33, 36, 39,
               42, 44, 46, 47, 48, 49, 50, 51])
        """
        g = self.grid
        nc = g.number_of_node_columns
        max_link_id = 3 * (nc - 1) + 2 * ((nc + 1) // 2) + nc // 2
        lower_active = logical_and(
            arange(g.number_of_links) < max_link_id, g.status_at_link == 0
        )
        boundary = logical_or(
            g.status_at_node[g.node_at_link_tail] != 0,
            g.status_at_node[g.node_at_link_head] != 0,
        )
        active_bnd = logical_and(boundary, g.status_at_link == 0)
        self.links_to_update = as_id_array(
            where(logical_or(lower_active, active_bnd))[0]
        )

    def _get_new_base_nodes(self, rock_state):
        """Return an array (or scalar) of states for the newly uplifted bottom
        inner row.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.hex_cts import HexCTS
        >>> from landlab.ca.celllab_cts import Transition
        >>> mg = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> nsd = {}  # node state dict
        >>> for i in range(10):
        ...     nsd[i] = i
        ...
        >>> xnlist = []
        >>> xnlist.append(Transition((0, 0, 0), (1, 1, 0), 1.0, "frogging"))
        >>> nsg = mg.add_zeros("node", "node_state_grid")
        >>> ca = HexCTS(mg, nsd, xnlist, nsg)

        >>> lu = LatticeUplifter(opt_block_layer=True)
        >>> lu._get_new_base_nodes(rock_state=7)
        array([9, 9, 9])
        >>> lu.uplift_interior_nodes(ca, current_time=0.0, rock_state=7)
        >>> lu.node_state[:5]
        array([0, 9, 0, 9, 9])
        >>> lu = LatticeUplifter(
        ...     opt_block_layer=True,
        ...     block_layer_thickness=2,
        ...     block_layer_dip_angle=90.0,
        ...     layer_left_x=1.0,
        ... )
        >>> lu._get_new_base_nodes(rock_state=7)
        array([9, 7, 9])
        >>> lu.uplift_interior_nodes(ca, current_time=0.0, rock_state=7)
        >>> lu.node_state[:5]
        array([0, 9, 0, 7, 9])
        >>> lu = LatticeUplifter(
        ...     opt_block_layer=True,
        ...     block_layer_thickness=1,
        ...     block_layer_dip_angle=45.0,
        ...     y0_top=-1.0,
        ... )
        >>> lu._get_new_base_nodes(rock_state=7)
        array([9, 7, 9])
        >>> lu.uplift_interior_nodes(ca, current_time=0.0, rock_state=7)
        >>> lu.node_state[:5]
        array([0, 9, 0, 7, 9])
        """

        new_base_nodes = zeros(len(self.inner_base_row_nodes), dtype=int)

        if self.block_layer_dip_angle == 0.0:  # horizontal
            if self.cum_uplift < self.block_layer_thickness:
                new_base_nodes[:] = self.block_ID
            else:
                new_base_nodes[:] = rock_state

        elif self.block_layer_dip_angle == 90.0:  # vertical
            layer_right_x = self.layer_left_x + self.block_layer_thickness
            inside_layer = where(
                logical_and(
                    self.grid.x_of_node[self.inner_base_row_nodes] >= self.layer_left_x,
                    self.grid.x_of_node[self.inner_base_row_nodes] <= layer_right_x,
                )
            )[0]
            new_base_nodes[:] = rock_state
            new_base_nodes[inside_layer] = self.block_ID

        else:
            x = self.grid.x_of_node[self.inner_base_row_nodes]
            y = self.grid.y_of_node[self.inner_base_row_nodes]
            m = tan(pi * self.block_layer_dip_angle / 180.0)
            y_top = m * x + self.y0_top
            y_bottom = y_top - (
                self.block_layer_thickness
                / cos(pi * self.block_layer_dip_angle / 180.0)
            )
            inside_layer = where(logical_and(y >= y_bottom, y <= y_top))
            new_base_nodes[:] = rock_state
            new_base_nodes[inside_layer] = self.block_ID

        return new_base_nodes

    def shift_link_and_transition_data_upward(self, ca, current_time):
        """Applies uplift to links and transitions.

        For each link that lies above the y = 1.5 cells line, assign the
        properties of the link one row down.

        (For an example, see unit test:
            test_shift_link_and_transition_data_upward)
        """

        # Find the ID of the first link above the y = 1.5 line
        nc = self.grid.number_of_node_columns
        first_link = (
            ((nc - 1) // 2)  # skip bottom horizontals
            + (3 * (nc - 1))  # skip 3 sets of diagonals
            + nc  # skip a full row of verticals
            + ((nc + 1) // 2)
        )  # skip a an even row of verticals

        # Define the offset in ID between a link and its neighbor one row up
        # (or down)
        shift = nc + 2 * (nc - 1)

        # Loop from top to bottom of grid, shifting the following link data
        # upward: state of link, ID of its next transition, and time of its
        # next transition.
        for lnk in range(self.grid.number_of_links - 1, first_link - 1, -1):
            ca.link_state[lnk] = ca.link_state[lnk - shift]
            ca.next_trn_id[lnk] = ca.next_trn_id[lnk - shift]
            ca.next_update[lnk] = ca.next_update[lnk - shift]

        # Sweep through event queue, shifting links upward. Do NOT shift links
        # with IDs greater than NL - [SHIFT + (NC - 1)], because these are so
        # close to the top of the grid that either the events would refer to
        # non-existent links (>= NL) or would involve shifting an event onto
        # an upper-boundary link. Note that because the event data are stored
        # in a tuple, we have to replace the entire tuple (can't simply change
        # the one item, because tuples are immutable)
        first_no_shift_id = self.grid.number_of_links - (shift + (nc - 1))
        for i in range(len(ca.priority_queue._queue)):
            if ca.priority_queue._queue[i][2] < first_no_shift_id:
                ca.priority_queue._queue[i] = (
                    ca.priority_queue._queue[i][0],
                    ca.priority_queue._queue[i][1],
                    (ca.priority_queue._queue[i][2] + shift),
                )

        # Update state of links along the boundaries.
        for lk in self.links_to_update:
            # Update link state
            fns = self.node_state[self.grid.node_at_link_tail[lk]]
            tns = self.node_state[self.grid.node_at_link_head[lk]]
            orientation = ca.link_orientation[lk]
            new_link_state = (
                orientation * ca.num_node_states_sq + fns * ca.num_node_states + tns
            )

            # Schedule a new transition, if applicable
            ca.link_state[lk] = new_link_state
            if ca.n_trn[new_link_state] > 0:
                (event_time, this_trn_id) = get_next_event_new(
                    lk, new_link_state, current_time, ca.n_trn, ca.trn_id, ca.trn_rate
                )
                ca.priority_queue.push(lk, event_time)
                ca.next_update[lk] = event_time
                ca.next_trn_id[lk] = this_trn_id
            else:
                ca.next_update[lk] = _NEVER
                ca.next_trn_id[lk] = -1

    def uplift_property_ids(self):
        """Shift property IDs upward by one row."""
        top_row_propid = self.propid[self.inner_top_row_nodes]
        for r in range(self.nr - 1, 0, -1):
            self.propid[self.inner_base_row_nodes + self.nc * r] = self.propid[
                self.inner_base_row_nodes + self.nc * (r - 1)
            ]
        self.propid[self.inner_base_row_nodes] = top_row_propid
        self.prop_data[self.propid[self.inner_base_row_nodes]] = self.prop_reset_value

    def uplift_interior_nodes(self, ca, current_time, rock_state=1):
        """Simulate 'vertical' displacement by shifting contents of node_state.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> from landlab.ca.hex_cts import HexCTS
        >>> from landlab.ca.celllab_cts import Transition
        >>> mg = HexModelGrid(
        ...     (5, 5), spacing=1.0, orientation="vertical", node_layout="rect"
        ... )
        >>> nsd = {}
        >>> for i in range(26):
        ...     nsd[i] = i
        ...
        >>> xnlist = []
        >>> xnlist.append(Transition((0, 0, 0), (1, 1, 0), 1.0, "frogging", True))
        >>> nsg = mg.add_zeros("node", "node_state_grid")
        >>> ca = HexCTS(mg, nsd, xnlist, nsg)
        >>> pd = mg.add_zeros("node", "propdata")
        >>> lu = LatticeUplifter(propid=ca.propid, prop_data=pd)
        >>> lu.node_state[:] = np.arange(len(lu.node_state))
        >>> lu.uplift_interior_nodes(ca, rock_state=25, current_time=0.0)
        >>> lu.node_state
        array([ 0, 25,  2, 25, 25,
                5,  1,  7,  3,  4,
               10,  6, 12,  8,  9,
               15, 11, 17, 13, 14,
               20, 16, 22, 18, 19])
        >>> lu.propid
        array([ 0, 21,  2, 23, 24,
                5,  1,  7,  3,  4,
               10,  6, 12,  8,  9,
               15, 11, 17, 13, 14,
               20, 16, 22, 18, 19])
        """

        # Shift the node states up by a full row. A "full row" includes two
        # staggered rows.
        for r in range(self.nr - 1, 0, -1):
            # This row gets the contents of the nodes 1 row down
            self.node_state[self.inner_base_row_nodes + self.nc * r] = self.node_state[
                self.inner_base_row_nodes + self.nc * (r - 1)
            ]

        # Fill the bottom rows with "fresh material" (code = rock_state), or
        # if using a block layer, with the right pattern of states.
        if self.opt_block_layer:
            new_base_nodes = self._get_new_base_nodes(rock_state)
            self.cum_uplift += 1.0
            self.y0_top += 1.0
        else:
            new_base_nodes = rock_state
        self.node_state[self.inner_base_row_nodes] = new_base_nodes

        # If propid (property ID) is defined, shift that too.
        if self.propid is not None:
            self.uplift_property_ids()

        self.shift_link_and_transition_data_upward(ca, current_time)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
