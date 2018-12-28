#! /usr/bin/env python

import numpy as np
from .ext.perimeternodes import (
    fill_perimeter_nodes_rect_horizontal,
    fill_perimeter_nodes_rect_vertical,
    fill_perimeter_nodes_hex_horizontal,
    fill_perimeter_nodes_hex_vertical,
)


def number_of_perimeter_nodes(shape, orientation="horizontal", node_layout="rect"):
    """Calculate the number of perimeter nodes on a hex grid.
    """
    if orientation == "horizontal":
        n_perimeter_nodes = 2 * shape[0] + 2 * (shape[1] - 2)
    elif orientation == "vertical":
        n_perimeter_nodes = 2 * shape[1] + 2 * (shape[0] - 2)
    else:
        raise ValueError(
            "orientation not understood. Must be one of 'horizontal' or 'vertical'"
        )

    if node_layout in ("hex", "rect1"):
        if orientation == "horizontal" and shape[0] % 2 == 0:
            n_perimeter_nodes += 1
            # n_perimeter_nodes += (shape[0] + 1) % 2
        elif orientation == "vertical" and shape[1] % 2 == 0:
            n_perimeter_nodes += 1
            # n_perimeter_nodes += (shape[1] + 1) % 2

    return n_perimeter_nodes


def perimeter_nodes(shape, orientation="horizontal", node_layout="rect", out=None):
    """Get the nodes along the perimeter of a hex grid.

    Parameters
    ----------
    shape : iterable of int
        Shape of the hex grid as (rows, columns).
    orientation : {"horizontal", "vertical"}
        Orientation of the hex grid: either 'horizontal' or 'vertical'.
    node_layout : {"rect", "hex"}
        Overall layout of the nodes.
    Returns
    -------
    ndarray of int
        Nodes along the perimeter of the grid.
    """

    if len(shape) != 2:
        raise ValueError("shape must be (rows, cols)")

    if node_layout not in ("hex", "rect"):
        raise ValueError( "node_layout not understood. Must be one of 'hex' or 'rect'")

    n_perimeter_nodes = number_of_perimeter_nodes(
        shape, orientation=orientation, node_layout=node_layout
    )

    if out is None:
        perimeter_nodes = np.empty(n_perimeter_nodes, dtype=int)
    else:
        perimeter_nodes = out

    if node_layout == "hex":
        if orientation == "horizontal":
            fill_perimeter_nodes_hex_horizontal(shape, perimeter_nodes)
        else:
            fill_perimeter_nodes_hex_vertical(shape, perimeter_nodes)
    else:
        if orientation == "horizontal":
            fill_perimeter_nodes_rect_horizontal(shape, perimeter_nodes)
        else:
            fill_perimeter_nodes_rect_vertical(shape, perimeter_nodes)

    return perimeter_nodes
