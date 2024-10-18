"""
IcosphereGlobalGrid class

Greg Tucker, CU Boulder, November 2023
"""

import numpy as np

from landlab.graph.quasi_spherical.dual_icosphere import DualIcosphereGraph
from landlab.grid.base import ModelGrid
from landlab.grid.nodestatus import NodeStatus


class IcosphereGlobalGrid(DualIcosphereGraph, ModelGrid):
    """
    Icosphere-based quasi-spherical ("global") Landlab grid.

    The default configuration is a spherical grid of unit radius that
    forms the spherical version an icosahedron (20 triangular patches,
    12 nodes), with the dual complement representing a dodecahedron
    (12 hexagonal cells, 20 corners). The mesh_densification_level
    parameter allows you to densify this initial shape by subdividing
    each triangular patch into four triangles, with corresponding
    addition of nodes (as the triangle vertices), together with
    corresponding cells and corners.

    If you specify mesh_densification_level=1, you get a soccer ball:
    a combination of pentagons and hexagons as cells. Further
    densification produces more hexagons (as cells; the patches are
    always triangles).

    Because there are no boundaries, there is a 1:1 relationship between
    nodes and cells (every node has a cell, and the ID of every cell is
    the same as the ID of its node), and similarly for corners and patches.

    Link length is calculated as the arc-length of the sphere between
    the link's two nodes. Patch area is calculated as the spherical
    (not flat) triangle area. Cell area is calculated by summing
    spherical triangles (pentagonal cell area is the sum of 5 triangles
    within the pentagon, while hexagonal cell area is the sum of 6
    triangles).

    Topography relative to the sphere can be created by adding a field
    that represents height relative to the sphere's surface. Topography
    could, if desired, be configured to represent an ellipsoid or
    geoid surface, but the patch/cell areas and link/face lengths are
    computed as if on a sphere (i.e., the present version of the
    component does not include algorithms to calculate distances or
    areas on an ellipsoid).

    The grid-generation and refinement algorithms are implemented in the
    superclass DualIcosphereGraph.

    Parameters
    ----------
    radius : float, optional
        Radius of the icosphere (default 1)
    mesh_densification_level : int, optional
        Number of times to densify the initial icosahedron (default 0)

    Examples
    --------
    >>> ico = IcosphereGlobalGrid()
    >>> ico.number_of_nodes
    12
    >>> ico.number_of_patches
    20
    >>> ico.number_of_corners
    20
    >>> ico.number_of_cells
    12
    """

    def __init__(self, radius=1.0, mesh_densification_level=0):
        """Initialize the IcosphereGlobalGrid"""
        DualIcosphereGraph.__init__(self, radius, mesh_densification_level)
        ModelGrid.__init__(self)

        self._node_status = np.full(
            self.number_of_nodes, NodeStatus.CORE, dtype=np.uint8
        )
