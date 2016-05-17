#! /usr/env/python
"""
Python implementation of VoronoiDelaunayGrid, a class used to create and manage
unstructured, irregular grids for 2D numerical models.

Getting Information about a Grid
--------------------------------
The following attributes, properties, and methods provide data about the grid,
its geometry, and the connectivity among the various elements. Each grid
element has an ID number, which is also its position in an array that
contains information about that type of element. For example, the *x*
coordinate of node 5 would be found at `grid.node_x[5]`.

The naming of grid-element arrays is *attribute*`_at_`*element*, where
*attribute* is the name of the data in question, and *element* is the element
to which the attribute applies. For example, the property `node_at_cell`
contains the ID of the node associated with each cell. For example,
`node_at_cell[3]` contains the *node ID* of the node associated with cell 3.
The *attribute* is singular if there is only one value per element; for
example, there is only one node associated with each cell. It is plural when
there are multiple values per element; for example, the `faces_at_cell` array
contains multiple faces for each cell. Exceptions to these general rules are
functions that return indices of a subset of all elements of a particular type.
For example, you can obtain an array with IDs of only the core nodes using
`core_nodes`, while `active_links` provides an array of IDs of active links
(only). Finally, attributes that represent a measurement of something, such as
the length of a link or the surface area of a cell, are described using `_of_`,
as in the example `area_of_cell`.

Information about the grid as a whole
+++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.axis_name
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.axis_units
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.move_origin
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_axis_coordinates
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_elements
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.save

Information about nodes
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_link_dirs_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.all_node_azimuths_map
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.all_node_distances_map
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.cell_area_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.cell_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.closed_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.downwind_links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_value_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_node_is_downwind
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_node_is_upwind
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_dirs_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.neighbors_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_core_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_link_head
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_link_tail
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_axis_coordinates
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_is_boundary
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_x
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_y
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.nodes_at_patch
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.open_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.patches_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_closed
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.unit_vector_sum_xcomponent_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.unit_vector_sum_ycomponent_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.upwind_links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.x_of_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.y_of_node

Information about links
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_link_dirs_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.angle_of_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.downwind_links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.face_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.length_of_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_face
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_node_is_downwind
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_node_is_upwind
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_dirs_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_link_head
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_link_tail
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_links_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.resolve_values_on_active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.resolve_values_on_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.unit_vector_xcomponent_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.unit_vector_ycomponent_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.upwind_links_at_node

Information about cells
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.area_of_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.cell_area_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.cell_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.faces_at_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_core_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_faces_at_cell

Information about faces
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.face_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.faces_at_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_face
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_faces_at_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.width_of_face

Information about patches
+++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.nodes_at_patch
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.patches_at_node

Data Fields in ModelGrid
------------------------
:class:`~.ModelGrid` inherits from the :class:`~.ModelDataFields` class. This
provides `~.ModelGrid`, and its subclasses, with the ability to, optionally,
store data values that are associated with the different types grid elements
(nodes, cells, etc.). In particular, as part of ``ModelGrid.__init__()``,
data field *groups* are added to the `ModelGrid` that provide containers to
put data fields into. There is one group for each of the eight grid elements
(node, cell, link, face, core_node, core_cell, active_link, and active_face).

To access these groups, use the same methods as accessing groups with
`~.ModelDataFields`. ``ModelGrid.__init__()`` adds the following attributes to
itself that provide access to the values groups:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.at_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.at_face

Each of these attributes returns a ``dict``-like object whose keys are value
names as strings and values are numpy arrays that gives quantities at
grid elements.


Create Field Arrays
+++++++++++++++++++
:class:`~.ModelGrid` inherits several useful methods for creating new data
fields and adding new data fields to a ModelGrid instance. Methods to add or
create a new data array follow the ``numpy`` syntax for creating arrays. The
folowing methods create and, optionally, initialize new arrays. These arrays
are of the correct size but a new field will not be added to the field:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.field.grouped.ModelDataFields.empty
    ~landlab.field.grouped.ModelDataFields.ones
    ~landlab.field.grouped.ModelDataFields.zeros

Add Fields to a ModelGrid
+++++++++++++++++++++++++
Unlike with the equivalent numpy functions, these do not take a size argument
as the size of the returned arrays is determined from the size of the
ModelGrid. However, the keyword arguments are the same as those of the numpy
equivalents.

The following methods will create a new array and add a reference to that
array to the ModelGrid:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.add_empty
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.add_field
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.add_ones
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.add_zeros
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.delete_field
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_units

These methods operate in the same way as the previous set except that, in
addition to creating a new array, the newly-created array is added to the
ModelGrid. The calling signature is the same but with the addition of an
argument that gives the name of the new field as a string. The additional
method, :meth:`~.ModelDataFields.add_field`, adds a previously allocation
array to the ModelGrid. If the array is of the incorrect size it will raise
``ValueError``.

Query Fields
++++++++++++
Use the following methods/attributes get information about the stored data
fields:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.field.grouped.ModelDataFields.size
    ~landlab.field.grouped.ModelDataFields.keys
    ~landlab.field.grouped.ModelDataFields.has_group
    ~landlab.field.grouped.ModelDataFields.has_field
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.field_units
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.field_values
    ~landlab.field.grouped.ModelDataFields.groups

    i.e., call, e.g. mg.has_field('node', 'my_field_name')

    # START HERE check that all functions listed below are included above,
    # ignore ones that start with underscores(_)

Gradients, fluxes, and divergences on the grid
----------------------------------------------

Landlab is designed to easily calculate gradients in quantities across the
grid, and to construct fluxes and flux divergences from them. Because these
calculations tend to be a little more involved than property lookups, the
methods tend to start with `calc_`.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_diff_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_flux_div_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_grad_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_grad_at_patch
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_net_flux_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_slope_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_slope_at_patch
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_unit_normal_at_patch

Mappers
-------

These methods allow mapping of values defined on one grid element type onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_downwind_node_link_max_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_downwind_node_link_mean_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_head_node_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_tail_node_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_vector_to_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_max_of_link_nodes_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_max_of_node_links_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_mean_of_link_nodes_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_min_of_link_nodes_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_min_of_node_links_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_node_to_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_upwind_node_link_max_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_upwind_node_link_mean_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_downwind_node_link_max_to_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_max_node_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_min_node_to_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_upwind_node_link_max_to_node


Boundary condition control
--------------------------

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to :meth:`~.ModelGrid.status_at_node` and
:meth:`~.ModelGrid.status_at_node` will automatically update the conditions
defined at other grid elements automatically.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.closed_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_value_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_core_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_is_boundary
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.open_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_closed
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_node

Identifying node subsets
------------------------

These methods are useful in identifying subsets of nodes, e.g., closest node
to a point; nodes at edges.

(None are available for this grid type)

Surface analysis
----------------

These methods permit the kinds of surface analysis that you might expect to
find in GIS software.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_aspect_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_hillshade_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_slope_at_node

Notes
-----
It is important that when creating a new grid class that inherits from
``ModelGrid``, to call ``ModelGrid.__init__()`` in the new grid's
``__init__()``. For example, the new class's __init__ should contain the
following code,

.. code-block:: python

    class NewGrid(ModelGrid):
        def __init__(self, *args, **kwds):
            ModelGrid.__init__(self, **kwds)
            # Code that initializes the NewGrid

Without this, the new grid class will not have the ``at_*`` attributes.
"""

import numpy
from six.moves import range

from landlab.grid.base import (ModelGrid, CORE_NODE, BAD_INDEX_VALUE,
                               INACTIVE_LINK)
from landlab.core.utils import (as_id_array, sort_points_by_x_then_y,
                                argsort_points_by_x_then_y,
                                anticlockwise_argsort_points)
from .decorators import return_readonly_id_array

from scipy.spatial import Voronoi


def simple_poly_area(x, y):
    """Calculates and returns the area of a 2-D simple polygon.

    Input vertices must be in sequence (clockwise or counterclockwise). *x*
    and *y* are arrays that give the x- and y-axis coordinates of the
    polygon's vertices.

    Parameters
    ----------
    x : ndarray
        x-coordinates of of polygon vertices.
    y : ndarray
        y-coordinates of of polygon vertices.

    Returns
    -------
    out : float
        Area of the polygon

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.voronoi import simple_poly_area
    >>> x = np.array([3., 1., 1., 3.])
    >>> y = np.array([1.5, 1.5, 0.5, 0.5])
    >>> simple_poly_area(x, y)
    2.0

    If the input coordinate arrays are 2D, calculate the area of each polygon.
    Note that when used in this mode, all polygons must have the same
    number of vertices, and polygon vertices are listed column-by-column.

    >>> x = np.array([[ 3.,  1.,  1.,  3.],
    ...               [-2., -2., -1., -1.]]).T
    >>> y = np.array([[1.5, 1.5, 0.5, 0.5],
    ...               [ 0.,  1.,  2.,  0.]]).T
    >>> simple_poly_area(x, y)
    array([ 2. ,  1.5])
    """
    # For short arrays (less than about 100 elements) it seems that the
    # Python sum is faster than the numpy sum. Likewise for the Python
    # built-in abs.
    return .5 * abs(sum(x[:-1] * y[1:] - x[1:] * y[:-1]) +
                    x[-1] * y[0] - x[0] * y[-1])


def calculate_link_lengths(pts, link_from, link_to):
    """Calculates and returns length of links between nodes.

    Parameters
    ----------
    pts : Nx2 numpy array containing (x,y) values
    link_from : 1D numpy array containing index numbers of nodes at starting
                point ("from") of links
    link_to : 1D numpy array containing index numbers of nodes at ending point
              ("to") of links

    Returns
    -------
    out : ndarray
        1D numpy array containing horizontal length of each link

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.voronoi import calculate_link_lengths
    >>> pts = np.array([[0.,0.],[3.,0.],[3.,4.]]) # 3:4:5 triangle
    >>> lfrom = np.array([0,1,2])
    >>> lto = np.array([1,2,0])
    >>> calculate_link_lengths(pts, lfrom, lto)
    array([ 3.,  4.,  5.])
    """
    dx = pts[link_to, 0] - pts[link_from, 0]
    dy = pts[link_to, 1] - pts[link_from, 1]
    link_length = numpy.sqrt(dx * dx + dy * dy)
    return link_length


class VoronoiDelaunayGrid(ModelGrid):
    """
    This inherited class implements an unstructured grid in which cells are
    Voronoi polygons and nodes are connected by a Delaunay triangulation. Uses
    scipy.spatial module to build the triangulation.

    Create an unstructured grid from points whose coordinates are given
    by the arrays *x*, *y*.

    Parameters
    ----------
    x : array_like
        x-coordinate of points
    y : array_like
        y-coordinate of points
    reorient_links (optional) : bool
        whether to point all links to the upper-right quadrant

    Returns
    -------
    VoronoiDelaunayGrid
        A newly-created grid.

    Examples
    --------
    >>> from numpy.random import rand
    >>> from landlab.grid import VoronoiDelaunayGrid
    >>> x, y = rand(25), rand(25)
    >>> vmg = VoronoiDelaunayGrid(x, y)  # node_x_coords, node_y_coords
    >>> vmg.number_of_nodes
    25

    >>> import numpy as np
    >>> x = [0, 0, 0, 0,
    ...      1, 1, 1, 1,
    ...      2, 2, 2, 2,]
    >>> y = [0, 1, 2, 3,
    ...      0, 1, 2, 3,
    ...      0, 1, 2, 3]
    >>> vmg = VoronoiDelaunayGrid(x, y)
    >>> vmg.node_x # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  1.,  2.,
            0.,  1.,  2.,
            0.,  1.,  2.,
            0.,  1.,  2.])
    >>> vmg.node_y # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  0.,
            1.,  1.,  1.,
            2.,  2.,  2.,
            3.,  3.,  3.])
    """

    def __init__(self, x=None, y=None, reorient_links=True, **kwds):
        """
        Create a Voronoi Delaunay grid from a set of points.

        Create an unstructured grid from points whose coordinates are given
        by the arrays *x*, *y*.

        Parameters
        ----------
        x : array_like
            x-coordinate of points
        y : array_like
            y-coordinate of points
        reorient_links (optional) : bool
            whether to point all links to the upper-right quadrant

        Returns
        -------
        VoronoiDelaunayGrid
            A newly-created grid.

        Examples
        --------
        >>> from numpy.random import rand
        >>> from landlab.grid import VoronoiDelaunayGrid
        >>> x, y = rand(25), rand(25)
        >>> vmg = VoronoiDelaunayGrid(x, y)  # node_x_coords, node_y_coords
        >>> vmg.number_of_nodes
        25
        """
        if (x is not None) and (y is not None):
            self._initialize(x, y, reorient_links)
        super(VoronoiDelaunayGrid, self).__init__(**kwds)

    def _initialize(self, x, y, reorient_links=True):
        """
        Creates an unstructured grid around the given (x,y) points.
        """
        x = numpy.asarray(x, dtype=float).reshape((-1, ))
        y = numpy.asarray(y, dtype=float).reshape((-1, ))

        if x.size != y.size:
            raise ValueError('x and y arrays must have the same size')

        # Make a copy of the points in a 2D array (useful for calls to geometry
        # routines, but takes extra memory space).
        pts = numpy.zeros((len(x), 2))
        pts[:, 0] = x
        pts[:, 1] = y
        self.pts = sort_points_by_x_then_y(pts)
        x = self.pts[:, 0]
        y = self.pts[:, 1]

        # NODES AND CELLS: Set up information pertaining to nodes and cells:
        #   - number of nodes
        #   - node x, y coordinates
        #   - default boundary status
        #   - interior and boundary nodes
        #   - nodes associated with each cell and active cell
        #   - cells and active cells associated with each node
        #     (or BAD_VALUE_INDEX if none)
        #
        # Assumptions we make here:
        #   - all interior (non-perimeter) nodes have cells (this should be
        #       guaranteed in a Delaunay triangulation, but there may be
        #       special cases)
        #   - all cells are active (later we'll build a mechanism for the user
        #       specify a subset of cells as active)
        #
        self._node_x = x
        self._node_y = y
        [self._node_status, self._core_nodes, self._boundary_nodes] = \
            self._find_perimeter_nodes_and_BC_set(pts)
        [self._cell_at_node, self._node_at_cell] = \
            self._node_to_cell_connectivity(self._node_status,
                                            self.number_of_cells)
        active_cell_at_node = self.cell_at_node[self.core_nodes]

        # ACTIVE CELLS: Construct Voronoi diagram and calculate surface area of
        # each active cell.
        vor = Voronoi(self.pts)
        self.vor = vor
        self._area_of_cell = numpy.zeros(self.number_of_cells)
        for node in self._node_at_cell:
            xv = vor.vertices[vor.regions[vor.point_region[node]], 0]
            yv = vor.vertices[vor.regions[vor.point_region[node]], 1]
            self._area_of_cell[self.cell_at_node[node]] = (
                simple_poly_area(xv, yv))

        # LINKS: Construct Delaunay triangulation and construct lists of link
        # "from" and "to" nodes.
        (self._node_at_link_tail,
         self._node_at_link_head,
         _,
         self._face_width) = \
            self._create_links_and_faces_from_voronoi_diagram(vor)
        self._status_at_link = numpy.full(len(self._node_at_link_tail),
                                          INACTIVE_LINK, dtype=int)

        # Sort them by midpoint coordinates
        self._sort_links_by_midpoint()

        # Optionally re-orient links so that they all point within upper-right
        # semicircle
        if reorient_links:
            self._reorient_links_upper_right()

        # LINKS: Calculate link lengths
        self._link_length = calculate_link_lengths(self.pts,
                                                   self.node_at_link_tail,
                                                   self.node_at_link_head)

        # LINKS: inlink and outlink matrices
        # SOON TO BE DEPRECATED
        self._setup_inlink_and_outlink_matrices()

        # ACTIVE LINKS: Create list of active links, as well as "from" and "to"
        # nodes of active links.
        self._reset_link_status_list()

        # NODES & LINKS: IDs and directions of links at each node
        self._create_links_and_link_dirs_at_node()

        # LINKS: set up link unit vectors and node unit-vector sums
        self._create_link_unit_vectors()

    @property
    def number_of_patches(self):
        """Number of patches.

        Returns the number of patches over the grid.
        """
        try:
            return self._number_of_patches
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self.pts, self.vor)
            return self._number_of_patches

    @property
    def nodes_at_patch(self):
        """Get the four nodes at the corners of each patch in a regular grid.
        """
        try:
            return self._nodes_at_patch
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self.pts, self.vor)
            return self._nodes_at_patch

    @property
    @return_readonly_id_array
    def patches_at_node(self):
        """
        Return a (nnodes, max_voronoi_polygon_sides) array of patches at nodes.

        The patches are returned in LL standard order (ccw from E), with any
        nonexistent patches recorded after the ids of existing faces.
        Nonexistent patches are ID'ed as -1.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> mg = HexModelGrid(3, 3)
        >>> mg.patches_at_node # doctest: +SKIP
        array([[ 0,  2, -1, -1, -1, -1],
               [ 1,  3,  0, -1, -1, -1],
               [ 4,  1, -1, -1, -1, -1],
               [ 5,  2, -1, -1, -1, -1],
               [ 6,  8,  5,  2,  0,  3],
               [ 7,  9,  6,  3,  1,  4],
               [ 7,  4, -1, -1, -1, -1],
               [ 5,  8, -1, -1, -1, -1],
               [ 8,  6,  9, -1, -1, -1],
               [ 9,  7, -1, -1, -1, -1]])
        """
        try:
            return self._patches_at_node
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self.pts, self.vor)
            return self._patches_at_node

    def _find_perimeter_nodes_and_BC_set(self, pts):
        """
        Uses a convex hull to locate the perimeter nodes of the Voronoi grid,
        then sets them as fixed value boundary nodes.
        It then sets/updates the various relevant node lists held by the grid,
        and returns *node_status*, *core_nodes*, *boundary_nodes*.
        """

        # Calculate the convex hull for the set of points
        from scipy.spatial import ConvexHull
        hull = ConvexHull(pts, qhull_options='Qc')  # see below why we use 'Qt'

        # The ConvexHull object lists the edges that form the hull. We need to
        # get from this list of edges the unique set of nodes. To do this, we
        # first flatten the list of vertices that make up all the hull edges
        # ("simplices"), so it becomes a 1D array. With that, we can use the
        # set() function to turn the array into a set, which removes duplicate
        # vertices. Then we turn it back into an array, which now contains the
        # set of IDs for the nodes that make up the convex hull.
        #   The next thing to worry about is the fact that the mesh perimeter
        # might contain nodes that are co-planar (that is, co-linear in our 2D
        # world). For example, if you make a set of staggered points for a
        # hexagonal lattice using make_hex_points(), there will be some
        # co-linear points along the perimeter. The ones of these that don't
        # form convex corners won't be included in convex_hull_nodes, but they
        # are nonetheless part of the perimeter and need to be included in
        # the list of boundary_nodes. To deal with this, we pass the 'Qt'
        # option to ConvexHull, which makes it generate a list of coplanar
        # points. We include these in our set of boundary nodes.
        convex_hull_nodes = numpy.array(list(set(hull.simplices.flatten())))
        coplanar_nodes = hull.coplanar[:, 0]
        boundary_nodes = as_id_array(numpy.concatenate(
            (convex_hull_nodes, coplanar_nodes)))

        # Now we'll create the "node_status" array, which contains the code
        # indicating whether the node is interior and active (=0) or a
        # boundary (=1). This means that all perimeter (convex hull) nodes are
        # initially flagged as boundary code 1. An application might wish to
        # change this so that, for example, some boundaries are inactive.
        node_status = numpy.zeros(len(pts[:, 0]), dtype=numpy.int8)
        node_status[boundary_nodes] = 1

        # It's also useful to have a list of interior nodes
        core_nodes = as_id_array(numpy.where(node_status == 0)[0])

        # save the arrays and update the properties
        self._node_status = node_status
        self._core_cells = numpy.arange(len(core_nodes), dtype=numpy.int)
        self._node_at_cell = core_nodes
        self._boundary_nodes = boundary_nodes

        # Return the results
        return node_status, core_nodes, boundary_nodes

    def _create_cell_areas_array(self):
        """Set up an array of cell areas."""
        self._cell_areas = self.active_cell_areas
        return self._cell_areas

    @staticmethod
    def _node_to_cell_connectivity(node_status, ncells):
        """Set up node connectivity.

        Creates and returns the following arrays:

        *  For each node, the ID of the corresponding cell, or
           BAD_INDEX_VALUE if the node has no cell.
        *  For each cell, the ID of the corresponding node.

        Parameters
        ----------
        node_status : ndarray of ints
            1D array containing the boundary status code for each node.
        ncells : ndarray of ints
            Number of cells (must equal the number of occurrences of CORE_NODE
            in node_status).

        Examples
        --------
        >>> from landlab import VoronoiDelaunayGrid as vdg
        >>> import numpy as np
        >>> from landlab.grid import BAD_INDEX_VALUE
        >>> ns = np.array([1, 0, 0, 1, 0])  # 3 interior, 2 boundary nodes
        >>> [node_cell, cell_node] = vdg._node_to_cell_connectivity(ns, 3)
        >>> node_cell[1:3]
        array([0, 1])
        >>> node_cell[0] == BAD_INDEX_VALUE
        True
        >>> cell_node
        array([1, 2, 4])
        """
        assert ncells == numpy.count_nonzero(node_status == CORE_NODE), \
            'ncells must equal number of CORE_NODE values in node_status'

        cell = 0
        node_cell = numpy.ones(len(node_status), dtype=int) * BAD_INDEX_VALUE
        cell_node = numpy.zeros(ncells, dtype=int)
        for node in range(len(node_cell)):
            if node_status[node] == CORE_NODE:
                node_cell[node] = cell
                cell_node[cell] = node
                cell += 1

        return node_cell, cell_node

    @staticmethod
    def _create_links_from_triangulation(tri):
        """Create links from a Delaunay triangulation.

        From a Delaunay Triangulation of a set of points, contained in a
        scipy.spatial.Delaunay object "tri", creates and returns:

        *  a numpy array containing the ID of the "from" node for each link
        *  a numpy array containing the ID of the "to" node for each link
        *  the number of links in the triangulation

        Examples
        --------
        >>> from scipy.spatial import Delaunay
        >>> import numpy as np
        >>> from landlab.grid import VoronoiDelaunayGrid as vdg
        >>> pts = np.array([[ 0., 0.], [ 1., 0.], [ 1., 0.87],
        ...                 [-0.5, 0.87], [ 0.5, 0.87], [ 0., 1.73],
        ...                 [ 1., 1.73]])
        >>> dt = Delaunay(pts)
        >>> [myfrom,myto,nl] = vdg._create_links_from_triangulation(dt)
        >>> print myfrom, myto, nl # doctest: +SKIP
        [5 3 4 6 4 3 0 4 1 1 2 6] [3 4 5 5 6 0 4 1 0 2 4 2] 12
        """

        # Calculate how many links there will be and create the arrays.
        #
        # The number of links equals 3 times the number of triangles minus
        # half the number of shared links. Finding out the number of shared
        # links is easy: for every shared link, there is an entry in the
        # tri.neighbors array that is > -1 (indicating that the triangle has a
        # neighbor opposite a given vertex; in other words, two triangles are
        # sharing an edge).
        num_shared_links = numpy.count_nonzero(tri.neighbors > -1)
        num_links = 3 * tri.nsimplex - num_shared_links // 2
        link_fromnode = numpy.zeros(num_links, dtype=int)
        link_tonode = numpy.zeros(num_links, dtype=int)

        # Sweep through the list of triangles, assigning "from" and "to" nodes
        # to the list of links.
        #
        # The basic algorithm works as follows. For each triangle, we will add
        # its 3 edges as links. However, we have to make sure that each shared
        # edge is added only once. To do this, we keep track of whether or not
        # each triangle has been processed yet using a boolean array called
        # "tridone". When we look at a given triangle, we check each vertex in
        # turn. If there is no neighboring triangle opposite that vertex, then
        # we need to add the corresponding edge. If there is a neighboring
        # triangle but we haven't processed it yet, we also need to add the
        # edge. If neither condition is true, then this edge has already been
        # added, so we skip it.
        link_id = 0
        tridone = numpy.zeros(tri.nsimplex, dtype=bool)
        for t in range(tri.nsimplex):  # loop over triangles
            for i in range(0, 3):       # loop over vertices & neighbors
                if tri.neighbors[t, i] == -1 or not tridone[
                        tri.neighbors[t, i]]:
                    link_fromnode[link_id] = tri.simplices[
                        t, numpy.mod(i + 1, 3)]
                    link_tonode[link_id] = tri.simplices[
                        t, numpy.mod(i + 2, 3)]
                    link_id += 1
            tridone[t] = True

        # save the results
        # self.node_at_link_tail = link_fromnode
        # self.node_at_link_head = link_tonode

        # Return the results
        return link_fromnode, link_tonode, num_links

    @staticmethod
    def _is_valid_voronoi_ridge(vor, n):

        SUSPICIOUSLY_BIG = 40000000.0
        return (vor.ridge_vertices[n][0] != -1 and
                vor.ridge_vertices[n][1] != -1 and
                numpy.amax(numpy.abs(vor.vertices[
                    vor.ridge_vertices[n]])) < SUSPICIOUSLY_BIG)

    @staticmethod
    def _create_links_and_faces_from_voronoi_diagram(vor):
        """
        From a Voronoi diagram object created by scipy.spatial.Voronoi(),
        builds and returns:
        1. Arrays of link tail and head nodes
        2. Array of link IDs for each active link
        3. Array containing with of each face

        Parameters
        ----------
        vor : scipy.spatial.Voronoi
            Voronoi object initialized with the grid nodes.

        Returns
        -------
        out : tuple of ndarrays
            - link_fromnode = "from" node for each link (len=num_links)
            - link_tonode   = "to" node for each link (len=num_links)
            - active_links  = link ID for each active link
                                                    (len=num_active_links)
            - face_width    = width of each face (len=num_active_links

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.grid import VoronoiDelaunayGrid as vdg
        >>> pts = np.array([[0., 0.], [1., 0.], [-0.5, 0.87], [0.5, 0.87],
        ...                 [1.5, 0.87], [0., 1.73], [1., 1.73]])
        >>> from scipy.spatial import Voronoi
        >>> vor = Voronoi(pts)
        >>> [tn,hn,al,fw] = vdg._create_links_and_faces_from_voronoi_diagram(vor)
        >>> tn
        array([0, 0, 0, 1, 1, 2, 3, 2, 3, 6, 6, 6])
        >>> hn
        array([1, 2, 3, 3, 4, 3, 4, 5, 5, 3, 4, 5])
        >>> al
        array([2, 3, 5, 6, 8, 9])
        >>> fw
        array([ 0.57669199,  0.57669199,  0.575973  ,  0.575973  ,  0.57836419,
                0.57836419])
        """
        # Each Voronoi "ridge" corresponds to a link. The Voronoi object has an
        # attribute ridge_points that contains the IDs of the nodes on either
        # side (including ridges that have one of their endpoints undefined).
        # So, we set the number of links equal to the number of ridges.
        num_links = len(vor.ridge_points)

        # Create the arrays for link from and to nodes
        link_fromnode = -numpy.ones(num_links, dtype=int)
        link_tonode = -numpy.ones(num_links, dtype=int)

        # Ridges along the perimeter of the grid will have one of their
        # endpoints undefined. The endpoints of each ridge are contained in
        # vor.ridge_vertices, and an undefined vertex is flagged with -1.
        # Ridges with both vertices defined correspond to faces and active
        # links, while ridges with an undefined vertex correspond to inactive
        # links. So, to find the number of active links, we subtract from the
        # total number of links the number of occurrences of an undefined
        # vertex.
        num_active_links = num_links \
            - numpy.count_nonzero(numpy.array(vor.ridge_vertices) == -1)

        # Create arrays for active links and width of faces (which are Voronoi
        # ridges).
        active_links = -numpy.ones(num_active_links, dtype=int)
        face_width = -numpy.ones(num_active_links)

        # Find the order to sort by link midpoints
        link_midpoints = numpy.zeros((num_links, 2))
        for i in range(num_links):
            link_midpoints[i][:] = (vor.points[vor.ridge_points[i, 0]] +
                                    vor.points[vor.ridge_points[i, 1]])/2.
        ind = argsort_points_by_x_then_y(link_midpoints)

        # Loop through the list of ridges. For each ridge, there is a link, and
        # its "from" and "to" nodes are the associated "points". In addition,
        # if the ridge endpoints are defined, we have a face and an active
        # link, so we add them to our arrays as well.
        j = 0
        for i in range(num_links):
            link_fromnode[i] = vor.ridge_points[ind[i], 0]
            link_tonode[i] = vor.ridge_points[ind[i], 1]
            face_corner1 = vor.ridge_vertices[ind[i]][0]
            face_corner2 = vor.ridge_vertices[ind[i]][1]
            # means it's a valid face
            if VoronoiDelaunayGrid._is_valid_voronoi_ridge(vor, ind[i]):
                dx = vor.vertices[face_corner2, 0] - \
                    vor.vertices[face_corner1, 0]
                dy = vor.vertices[face_corner2, 1] - \
                    vor.vertices[face_corner1, 1]
                face_width[j] = numpy.sqrt(dx * dx + dy * dy)
                active_links[j] = i
                j += 1

        return link_fromnode, link_tonode, active_links, face_width

    def _reorient_links_upper_right(self):
        """Reorient links to all point within the upper-right semi-circle.

        Notes
        -----
        "Upper right semi-circle" means that the angle of the link with respect
        to the vertical (measured clockwise) falls between -45 and +135. More
        precisely, if :math:`\theta' is the angle,
        :math:`-45 \ge \theta < 135`.
        For example, the link could point up and left as much as -45, but not
        -46. It could point down and right as much as 134.9999, but not 135. It
        will never point down and left, or up-but-mostly-left, or
        right-but-mostly-down.

        Examples
        --------
        >>> from landlab.grid import HexModelGrid
        >>> hg = HexModelGrid(3, 2, 1., reorient_links=True)
        >>> hg.node_at_link_tail
        array([0, 0, 0, 1, 1, 2, 3, 2, 3, 3, 4, 5])
        >>> hg.node_at_link_head
        array([1, 2, 3, 3, 4, 3, 4, 5, 5, 6, 6, 6])
        """

        # Calculate the horizontal (dx) and vertical (dy) link offsets
        link_dx = self.node_x[self.node_at_link_head] - \
            self.node_x[self.node_at_link_tail]
        link_dy = self.node_y[self.node_at_link_head] - \
            self.node_y[self.node_at_link_tail]

        # Calculate the angle, clockwise, with respect to vertical, then rotate
        # by 45 degrees counter-clockwise (by adding pi/4)
        link_angle = numpy.arctan2(link_dx, link_dy) + numpy.pi / 4

        # The range of values should be -180 to +180 degrees (but in radians).
        # It won't be after the above operation, because angles that were
        # > 135 degrees will now have values > 180. To correct this, we
        # subtract 360 (i.e., 2 pi radians) from those that are > 180 (i.e.,
        # > pi radians).
        link_angle -= 2 * numpy.pi * (link_angle >= numpy.pi)

        # Find locations where the angle is negative; these are the ones we
        # want to flip
        (flip_locs, ) = numpy.where(link_angle < 0.)

        # If there are any flip locations, proceed to switch their fromnodes
        # and tonodes; otherwise, we're done
        if len(flip_locs) > 0:

            # Temporarily story the fromnode for these
            fromnode_temp = self.node_at_link_tail[flip_locs]

            # The fromnodes now become the tonodes, and vice versa
            self._node_at_link_tail[
                flip_locs] = self.node_at_link_head[flip_locs]
            self._node_at_link_head[flip_locs] = fromnode_temp

    def _create_patches_from_delaunay_diagram(self, pts, vor):
        """
        Uses a delaunay diagram drawn from the provided points to
        generate an array of patches and patch-node-link connectivity.
        Returns ...
        DEJH, 10/3/14, modified May 16.
        """
        from scipy.spatial import Delaunay
        tri = Delaunay(pts)
        assert numpy.array_equal(tri.points, vor.points)
        nodata = -1
        self._nodes_at_patch = tri.simplices
        # self._nodes_at_patch = numpy.empty_like(_nodes_at_patch)
        self._number_of_patches = tri.simplices.shape[0]
        # get the patches in order:
        patches_xy = numpy.empty((self._number_of_patches, 2), dtype=float)
        patches_xy[:, 0] = numpy.mean(self.node_x[self._nodes_at_patch],
                                      axis=1)
        patches_xy[:, 1] = numpy.mean(self.node_y[self._nodes_at_patch],
                                      axis=1)
        orderforsort = argsort_points_by_x_then_y(patches_xy)
        self._nodes_at_patch = self._nodes_at_patch[orderforsort, :]
        patches_xy = patches_xy[orderforsort, :]
        # get the nodes around the patch in order:
        nodes_xy = numpy.empty((3, 2), dtype=float)
        for i in range(self._number_of_patches):
            these_nodes = self._nodes_at_patch[i]
            nodes_xy[:, 0] = self.node_x[these_nodes]
            nodes_xy[:, 1] = self.node_y[these_nodes]
            sortorder = anticlockwise_argsort_points(nodes_xy)
            try:
                self._nodes_at_patch[i, :] = these_nodes[sortorder]
            except TypeError:  # sortorder was an int
                pass
        max_dimension = 0
        # need to build a squared off, masked array of the patches_at_node
        # the max number of patches for a node in the grid is the max sides of
        # the side-iest voronoi region.
        for i in range(len(vor.regions)):
            if len(vor.regions[i]) > max_dimension:
                max_dimension = len(vor.regions[i])
        _patches_at_node = numpy.empty(
            (self.number_of_nodes, max_dimension), dtype=int)
        _patches_at_node.fill(nodata)
        for i in range(self.number_of_nodes):
            patches_with_node = numpy.argwhere(
                numpy.equal(self._nodes_at_patch, i))[:, 0]
            _patches_at_node[
                i, :patches_with_node.size] = patches_with_node[:]
        self._patches_at_node = _patches_at_node
        self._patches_created = True

    def save(self, path, clobber=False):
        """Save a grid and fields.

        This method uses cPickle to save a Voronoi grid as a cPickle file.
        At the time of coding, this is the only convenient output format
        for Voronoi grids, but support for netCDF is likely coming.

        All fields will be saved, along with the grid.

        The recommended suffix for the save file is '.grid'. This will
        be added to your save if you don't include it.

        This method is equivalent to
        :py:func:`~landlab.io.native_landlab.save_grid`, and
        :py:func:`~landlab.io.native_landlab.load_grid` can be used to
        load these files.

        Caution: Pickling can be slow, and can produce very large files.
        Caution 2: Future updates to Landlab could potentially render old
        saves unloadable.

        Parameters
        ----------
        path : str
            Path to output file.
        clobber : bool (defaults to false)
            Set to true to allow overwriting

        Examples
        --------
        >>> from landlab import VoronoiDelaunayGrid
        >>> import numpy as np
        >>> import os
        >>> x = np.random.rand(20)
        >>> y = np.random.rand(20)
        >>> vmg = VoronoiDelaunayGrid(x,y)
        >>> vmg.save('./mytestsave.grid')
        >>> os.remove('mytestsave.grid') #to remove traces of this test
        """
        import os
        from six.moves import cPickle

        if os.path.exists(path) and not clobber:
            raise ValueError('file exists')

        (base, ext) = os.path.splitext(path)
        if ext != '.grid':
            ext = ext + '.grid'
        path = base + ext

        with open(path, 'wb') as fp:
            cPickle.dump(self, fp)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
