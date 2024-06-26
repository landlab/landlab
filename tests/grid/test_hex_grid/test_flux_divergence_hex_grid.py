"""
Created on Wed Nov  4 08:10:22 2015

@author: gtucker

List of needs:

grid.link_at_face: IDs of links at each face
gradient functions:
    - of a node value at faces
    - of a node value at links
    - of a node value at junctions
a consistent way to sort
    - nodes: by y then x
    - links and faces: by junction coords (y then x)


Notes:

    grads_at_links:
        lg = (nv[node_at_link_head]-nv[node_at_link_tail])/link_length
    grads_at_faces:
        fg = lg[link_at_face]
        OR
        fg = (nv[node_at_face_right]-nv[node_at_face_left])/link_length[link_at_face]
        OR
        fg = (nv[node_at_face_right]-nv[node_at_face_left])/link_length_face

CURRENT STICKING POINT, 11/6/15:
WANT TO IMPLEMENT THE 3 VERSIONS ABOVE. BUT NO PROPER FACE-LINK DATA STRUCTURE,
AND CELL AREAS SEEMS TO BE BROKEN...
RELATES TO THE NEED TO IDENTIFY PERIMETER NODES IN A VORONOI GRID.
IN A RECTANGULAR RASTER, NODES WITH CELLS ARE EASY TO IDENTIFY.
IN A HEX-SHAPED OR RECTANGULAR HEX, SLIGHTLY MORE COMPLICATED BUT NOT TOO BAD.
IN A GENERAL VORONOI...??
FIRST THING PROBABLY IS TO FIX LINK_AT_FACE FOR RASTER, THEN IMPLEMENT FOR HEX.
"""

import numpy as np
from numpy.testing import assert_array_equal

from landlab import HexModelGrid

MAX_NUM_LINKS = 6


def gt_grads_at_faces1(grid, nv):
    # lg = (nv[grid.node_at_link_head] - nv[grid.node_at_link_tail]) / grid.length_of_link
    # return lg[grid.link_at_face]
    return None  # temporary


def gt_link_flux_divergence_at_cells_with_2darray(grid, f, out=None):
    if out is not None:
        net_flux = out
    else:
        net_flux = np.zeros(grid.number_of_nodes)

    for r in range(MAX_NUM_LINKS):
        links = grid.gt_links_at_node[r, grid.node_at_cell]
        net_flux += f[links] * grid.gt_link_dirs_at_node[r, :]


def gt_calc_gradients_at_faces(grid, vn):
    # vn = values at nodes
    print(grid.link_at_face)


def make_links_at_node_array(grid):
    """Make array with links at each node"""

    # Create arrays for link-at-node information
    grid.gt_links_at_node = -np.ones(
        (MAX_NUM_LINKS, grid.number_of_nodes), dtype=np.int32
    )
    grid.gt_link_dirs_at_node = np.zeros(
        (MAX_NUM_LINKS, grid.number_of_nodes), dtype=np.int8
    )
    grid.gt_active_link_dirs_at_node = np.zeros(
        (MAX_NUM_LINKS, grid.number_of_nodes), dtype=np.int8
    )
    grid.gt_num_links_at_node = np.zeros(
        grid.number_of_nodes, dtype=np.uint8
    )  # assume <256 links at any node
    grid.gt_num__active_links_at_node = np.zeros(
        grid.number_of_nodes, dtype=np.uint8
    )  # assume <256 links at any node

    # Sweep over all links
    for lk in range(grid.number_of_links):
        # Find the ID of the tail node
        t = grid.node_at_link_tail[lk]

        # Its row in the 2D array is equal to the number of links we've found
        # so far for this node
        index = grid.gt_num_links_at_node[t]

        # Add this link to the list for this node, set the direction (outgoing,
        # indicated by -1), and increment the number found so far
        grid.gt_links_at_node[index][t] = lk
        grid.gt_link_dirs_at_node[index][t] = -1
        grid.gt_num_links_at_node[t] += 1

        # Find the ID of the tail node
        h = grid.node_at_link_head[lk]

        # Its row in the 2D array is equal to the number of links we've found
        # so far for this node
        index = grid.gt_num_links_at_node[h]

        # Add this link to the list for this node, set the direction
        # (incoming, indicated by +1), and increment the number found so far
        grid.gt_links_at_node[index][h] = lk
        grid.gt_link_dirs_at_node[index][h] = 1
        grid.gt_num_links_at_node[h] += 1


def testing_flux_divergence_with_hex():
    r"""Test flux divergence function(s).

    Notes
    -----
    Test grid looks like this:

        (7)-17-.(8)-18-.(9)
        . .     . .     . .
      11  12  13  14  15  16
      /     \ /     \ /     \
    (3)--8-.(4)--9-.(5)-10-.(6)
      .     . .     . .     .
       2   3   4   5   6   7
        \ /     \ /     \ /
        (0)--0-.(1)--1-.(2)

    Node numbers in parentheses; others are link numbers; period indicates
    link head.
    """
    hmg = HexModelGrid((3, 3), reorient_links=True)

    f = hmg.add_zeros("test_flux", at="link")
    f[:] = np.arange(hmg.number_of_links)

    make_links_at_node_array(hmg)

    assert_array_equal(hmg.gt_num_links_at_node, [3, 4, 3, 3, 6, 6, 3, 3, 4, 3])

    assert_array_equal(
        hmg.gt_links_at_node,
        [
            [0, 0, 1, 2, 3, 5, 7, 11, 13, 15],
            [2, 1, 6, 8, 4, 6, 10, 12, 14, 16],
            [3, 4, 7, 11, 8, 9, 16, 17, 17, 18],
            [-1, 5, -1, -1, 9, 10, -1, -1, 18, -1],
            [-1, -1, -1, -1, 12, 14, -1, -1, -1, -1],
            [-1, -1, -1, -1, 13, 15, -1, -1, -1, -1],
        ],
    )

    assert_array_equal(
        hmg.gt_link_dirs_at_node,
        [
            [-1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [-1, -1, -1, -1, 1, 1, 1, 1, 1, 1],
            [-1, -1, -1, -1, 1, 1, -1, -1, 1, 1],
            [0, -1, 0, 0, -1, -1, 0, 0, -1, 0],
            [0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, -1, -1, 0, 0, 0, 0],
        ],
    )

    raw_net_flux = np.zeros(hmg.number_of_nodes)
    for r in range(MAX_NUM_LINKS):
        raw_net_flux += f[hmg.gt_links_at_node[r, :]] * hmg.gt_link_dirs_at_node[r, :]
    assert_array_equal(
        raw_net_flux, [-5.0, -10.0, -12.0, -17.0, -19.0, -19.0, 1.0, 6.0, 26.0, 49.0]
    )

    nv = np.arange(hmg.number_of_nodes)
    # gt_calc_gradients_at_faces(hmg, nv)
    # gt_link_flux_divergence_at_cells_with_2darray(hmg, f)

    for _ in range(1000):
        gt_grads_at_faces1(hmg, nv)


if __name__ == "__main__":
    testing_flux_divergence_with_hex()
