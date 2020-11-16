#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Grid element mappers that are specific to hex grids.

Mapping functions unique to hex grids
+++++++++++++++++++++++++++++++++++++

.. autosummary::

    ~landlab.grid.hex_mappers.map_link_vector_components_to_node_hex
"""
import enum

import numpy as np


class LinkAtNode(enum.IntEnum):
    EAST = 0  # array column of link to east/right of node, horizontal grid
    ENE = 0  # array column of link to ene/upper right of node, vertical grid
    NNE = 1  # array column of link to nne/upper right of node, horizontal grid
    NORTH = 1  # array column of link to north/top of node, vertical grid
    NNW = 2  # array column of link to nnw/upper left of node, horizontal grid
    WNW = 2  # array column of link to wnw/upper left of node, vertical grid
    WEST = 3  # array column of link to west/left of node, horizontal grid
    WSW = 3  # array column of link to wsw/lower left of node, vertical grid
    SSW = 4  # array column of link to ssw/lower left of node, horizontal grid
    SOUTH = 4  # array column of link to south/bottom of node, vertical grid
    SSE = 5  # array column of link to sse/lower right of node, horizontal grid
    ESE = 5  # array column of link to ese/lower right of node, vertical grid


SIN60 = np.sin(np.deg2rad(60.0))


def map_link_vector_components_to_node_hex(grid, data_at_link):
    """Map (x,y) components of link data data_at_link onto nodes of hex grid.

    Parameters
    ----------
    grid : HexModelGrid
        Landlab HexModelGrid object
    data_at_links : ndarray of float x number of links
        Data to be mapped

    Returns
    -------
    (x_component, y_component) : tuple of ndarray
        The *x* and *y* components of the field at each node. Both *x*
        and *y* components for non-*core* nodes are set to zero.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import HexModelGrid
    >>> from landlab.grid.mappers import map_link_vector_components_to_node

    >>> grid = HexModelGrid((3, 3))
    >>> link_data = np.full(grid.number_of_links, 0.5 * 3.0 ** 0.5)
    >>> link_data[np.isclose(grid.angle_of_link, 0.0)] = 0.0

    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> vx
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
    >>> vy
    array([ 0.,  0.,  0.,  0.,  1.,  1.,  0.,  0.,  0.,  0.])

    >>> link_data = np.arange(grid.number_of_links)
    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> vx
    array([ 0. ,  0. ,  0. ,  0. ,  8.5,  9.5,  0. ,  0. ,  0. ,  0. ])

    >>> link_data = np.full(grid.number_of_links, 0.5 * 3.0 ** 0.5)
    >>> link_data[np.isclose(grid.angle_of_link, 2.0 / 3.0 * np.pi)] = 0.0

    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> np.round(vx, 3)
    array([ 0. ,  0. ,  0. ,  0. ,  0.866,  0.866,  0. ,  0. ,  0. ,  0. ])
    >>> vy
    array([ 0. ,  0. ,  0. ,  0. ,  0.5,  0.5,  0. ,  0. ,  0. ,  0. ])

    >>> grid = HexModelGrid((3, 3), orientation='vertical')
    >>> link_data = np.arange(grid.number_of_links)

    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> vy
    array([ 0. , 0. ,  0. ,  5.5,  0. ,  0. , 12.5,  0. ,  0. ,  0. ])

    >>> link_data = np.full(grid.number_of_links, 0.5 * 3.0 ** 0.5)
    >>> link_data[np.isclose(grid.angle_of_link, np.pi / 2.0)] = 0.0

    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> vx
    array([ 0.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  0.])
    >>> vy
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

    Notes
    -----
    Calculation is only made for core nodes; boundary nodes receive zeros.
        For a grid with orientation='horizontal', one of the 3 link orientations
    is horizontal. The x component is therefore taken as the average value of
    the links to the east (right) and west (left). For a grid with
    orientation='vertical', the same principle applies for the y component:
    it is the average of the link values to the north/top and south/bottom.
        In general, the theory behind the approach is that there exists a
    "true" vector field that has been projected onto the links.
    Let $V = (v_x, v_y)$ be the true vector. The projection of $V$ onto two links
    with unit vectors $l_1 = (l_{1x}, l_{1y})$ and $l_2 = (l_{2x}, l_{2y})$
    yields the following formula for the scalar magnitude, or component, of
    the two link vectors:

        L_1 = V dot l_1

        L_2 = V dot l_2

    We know $L_1$ and $L_2$: they are the values associated with two adjacent
    links, and here we're interested either in the nne and nnw oriented links
    (in a horizontally oriented grid) or the ene and ese oriented links (on a
    vertically oriented grid). We also know the unit vectors $l_1$ and $l_2$,
    which derive from the link orientations (they involve sin 60 and
    cos 60). So we have two equations with two unknowns: the vector
    components $v_x$ and $v_y$.
        In practice, we use this math to obtain the $y$ component for a
    horizontal grid, and the $x$ component for a vertical grid. The opposite
    component is found directly from the horizontal or vertical links,
    respectively.
        Note that in the above doc tests, we take advantage of the fact that
    sin 60 deg = half the square root of 3 (no need to import math or numpy).
    """
    cores = grid.core_nodes

    x_component = np.zeros(grid.number_of_nodes)
    y_component = np.zeros(grid.number_of_nodes)

    if grid.orientation[0] == "h":
        east = grid.links_at_node[cores, LinkAtNode.EAST]
        nne = grid.links_at_node[cores, LinkAtNode.NNE]
        nnw = grid.links_at_node[cores, LinkAtNode.NNW]
        west = grid.links_at_node[cores, LinkAtNode.WEST]
        sse = grid.links_at_node[cores, LinkAtNode.SSE]
        ssw = grid.links_at_node[cores, LinkAtNode.SSW]

        x_component[cores] = (data_at_link[west] + data_at_link[east]) / 2
        vyn = (3.0 * data_at_link[nnw] - data_at_link[nne]) / (2 * SIN60)
        vys = (3.0 * data_at_link[ssw] - data_at_link[sse]) / (2 * SIN60)
        y_component[cores] = (vyn + vys) / 2
    else:
        ene = grid.links_at_node[cores, LinkAtNode.ENE]
        north = grid.links_at_node[cores, LinkAtNode.NORTH]
        wnw = grid.links_at_node[cores, LinkAtNode.WNW]
        wsw = grid.links_at_node[cores, LinkAtNode.WSW]
        south = grid.links_at_node[cores, LinkAtNode.SOUTH]
        ese = grid.links_at_node[cores, LinkAtNode.ESE]
        y_component[cores] = (data_at_link[north] + data_at_link[south]) / 2

        vxe = (3.0 * data_at_link[ese] - data_at_link[ene]) / (2 * SIN60)
        vxw = (3.0 * data_at_link[wsw] - data_at_link[wnw]) / (2 * SIN60)
        x_component[cores] = (vxe + vxw) / 2

    return x_component, y_component
