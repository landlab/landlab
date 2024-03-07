#!/usr/bin/env python2
"""
Created on Fri Jun 15 09:31:10 2018

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal

from landlab import HexModelGrid
from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeUplifter
from landlab.ca.celllab_cts import Transition
from landlab.ca.oriented_hex_cts import OrientedHexCTS


def test_links_to_update():
    """Test that update list includes lower 2 rows and fault-crossing links"""

    # Create a 6x6 test grid
    hg = HexModelGrid((6, 6), node_layout="rect", orientation="vertical")

    lnf = LatticeNormalFault(grid=hg, fault_x_intercept=-0.1)

    assert_array_equal(
        lnf.links_to_update,
        [
            6,
            7,
            9,
            10,
            11,
            12,
            13,
            14,
            16,
            17,
            18,
            19,
            20,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            33,
            34,
            35,
            36,
            38,
            39,
            41,
            44,
            49,
            52,
            54,
            58,
            60,
            66,
            68,
            71,
            74,
            75,
            78,
        ],
    )


def test_shift_link_and_transition_data_upward():
    """Test the LatticeUplifter method that uplifts link data and tr'ns."""

    mg = HexModelGrid((4, 3), spacing=1.0, orientation="vertical", node_layout="rect")
    nsd = {0: "yes", 1: "no"}
    xnlist = []
    xnlist.append(Transition((0, 0, 0), (1, 1, 0), 1.0, "frogging"))
    xnlist.append(Transition((0, 0, 1), (1, 1, 1), 1.0, "frogging"))
    xnlist.append(Transition((0, 0, 2), (1, 1, 2), 1.0, "frogging"))
    nsg = mg.add_zeros("node_state_grid", at="node")
    ohcts = OrientedHexCTS(mg, nsd, xnlist, nsg)

    assert_array_equal(
        ohcts.link_state[mg.active_links], [0, 4, 8, 8, 4, 0, 4, 8, 8, 4, 0]
    )

    assert_array_equal(
        ohcts.next_trn_id[mg.active_links], [0, 1, 2, 2, 1, 0, 1, 2, 2, 1, 0]
    )

    assert_array_equal(
        np.round(ohcts.next_update[mg.active_links], 2),
        [0.8, 1.26, 0.92, 0.79, 0.55, 1.04, 0.58, 2.22, 3.31, 0.48, 1.57],
    )

    pq = ohcts.priority_queue

    assert_equal(pq._queue[0][2], 19)  # link for first event = 19, not shifted
    assert_equal(round(pq._queue[0][0], 2), 0.48)  # trn scheduled for t = 0.48
    assert_equal(pq._queue[2][2], 14)  # this event scheduled for link 15...
    assert_equal(round(pq._queue[2][0], 2), 0.58)  # ...trn sched for t = 0.58

    lu = LatticeUplifter(grid=mg)
    lu.shift_link_and_transition_data_upward(ohcts, 0.0)

    # note new events lowest 5 links
    assert_array_equal(
        np.round(ohcts.next_update[mg.active_links], 2),
        [0.75, 0.84, 2.6, 0.07, 0.09, 0.8, 0.02, 1.79, 1.51, 2.04, 3.85],
    )
    assert_equal(pq._queue[0][2], 14)  # new soonest event
    assert_equal(pq._queue[9][2], 13)  # was previously 7, now shifted up...
    assert_equal(round(pq._queue[9][0], 2), 0.8)  # ...still sched for t = 0.80
