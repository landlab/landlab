#! /usr/env/python
"""
Python implementation of ModelGrid, a base class used to create and manage
grids for 2D numerical models.

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

    ~landlab.grid.radial.RadialModelGrid.axis_name
    ~landlab.grid.radial.RadialModelGrid.axis_units
    ~landlab.grid.radial.RadialModelGrid.from_dict
    ~landlab.grid.radial.RadialModelGrid.move_origin
    ~landlab.grid.radial.RadialModelGrid.node_axis_coordinates
    ~landlab.grid.radial.RadialModelGrid.number_of_active_faces
    ~landlab.grid.radial.RadialModelGrid.number_of_active_links
    ~landlab.grid.radial.RadialModelGrid.number_of_cells
    ~landlab.grid.radial.RadialModelGrid.number_of_core_cells
    ~landlab.grid.radial.RadialModelGrid.number_of_core_nodes
    ~landlab.grid.radial.RadialModelGrid.number_of_elements
    ~landlab.grid.radial.RadialModelGrid.number_of_faces
    ~landlab.grid.radial.RadialModelGrid.number_of_fixed_links
    ~landlab.grid.radial.RadialModelGrid.number_of_links
    ~landlab.grid.radial.RadialModelGrid.number_of_nodes
    ~landlab.grid.radial.RadialModelGrid.number_of_nodes_in_shell
    ~landlab.grid.radial.RadialModelGrid.number_of_patches
    ~landlab.grid.radial.RadialModelGrid.number_of_shells
    ~landlab.grid.radial.RadialModelGrid.save
    ~landlab.grid.radial.RadialModelGrid.spacing_of_shells


Information about nodes
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.radial.RadialModelGrid.active_link_dirs_at_node
    ~landlab.grid.radial.RadialModelGrid.all_node_azimuths_map
    ~landlab.grid.radial.RadialModelGrid.all_node_distances_map
    ~landlab.grid.radial.RadialModelGrid.boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.cell_area_at_node
    ~landlab.grid.radial.RadialModelGrid.cell_at_node
    ~landlab.grid.radial.RadialModelGrid.closed_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.core_nodes
    ~landlab.grid.radial.RadialModelGrid.downwind_links_at_node
    ~landlab.grid.radial.RadialModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.link_at_node_is_downwind
    ~landlab.grid.radial.RadialModelGrid.link_at_node_is_upwind
    ~landlab.grid.radial.RadialModelGrid.link_dirs_at_node
    ~landlab.grid.radial.RadialModelGrid.links_at_node
    ~landlab.grid.radial.RadialModelGrid.neighbors_at_node
    ~landlab.grid.radial.RadialModelGrid.node_at_cell
    ~landlab.grid.radial.RadialModelGrid.node_at_core_cell
    ~landlab.grid.radial.RadialModelGrid.node_at_link_head
    ~landlab.grid.radial.RadialModelGrid.node_at_link_tail
    ~landlab.grid.radial.RadialModelGrid.node_axis_coordinates
    ~landlab.grid.radial.RadialModelGrid.node_is_boundary
    ~landlab.grid.radial.RadialModelGrid.node_x
    ~landlab.grid.radial.RadialModelGrid.node_y
    ~landlab.grid.radial.RadialModelGrid.nodes
    ~landlab.grid.radial.RadialModelGrid.nodes_at_patch
    ~landlab.grid.radial.RadialModelGrid.number_of_core_nodes
    ~landlab.grid.radial.RadialModelGrid.number_of_links_at_node
    ~landlab.grid.radial.RadialModelGrid.number_of_nodes
    ~landlab.grid.radial.RadialModelGrid.number_of_nodes_in_shell
    ~landlab.grid.radial.RadialModelGrid.open_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.patches_at_node
    ~landlab.grid.radial.RadialModelGrid.radius_at_node
    ~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.radial.RadialModelGrid.status_at_node
    ~landlab.grid.radial.RadialModelGrid.unit_vector_sum_xcomponent_at_node
    ~landlab.grid.radial.RadialModelGrid.unit_vector_sum_ycomponent_at_node
    ~landlab.grid.radial.RadialModelGrid.upwind_links_at_node
    ~landlab.grid.radial.RadialModelGrid.x_of_node
    ~landlab.grid.radial.RadialModelGrid.y_of_node


Information about links
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.radial.RadialModelGrid.active_link_dirs_at_node
    ~landlab.grid.radial.RadialModelGrid.active_links
    ~landlab.grid.radial.RadialModelGrid.angle_of_link
    ~landlab.grid.radial.RadialModelGrid.downwind_links_at_node
    ~landlab.grid.radial.RadialModelGrid.face_at_link
    ~landlab.grid.radial.RadialModelGrid.fixed_links
    ~landlab.grid.radial.RadialModelGrid.length_of_link
    ~landlab.grid.radial.RadialModelGrid.link_at_face
    ~landlab.grid.radial.RadialModelGrid.link_at_node_is_downwind
    ~landlab.grid.radial.RadialModelGrid.link_at_node_is_upwind
    ~landlab.grid.radial.RadialModelGrid.link_dirs_at_node
    ~landlab.grid.radial.RadialModelGrid.links_at_node
    ~landlab.grid.radial.RadialModelGrid.node_at_link_head
    ~landlab.grid.radial.RadialModelGrid.node_at_link_tail
    ~landlab.grid.radial.RadialModelGrid.number_of_active_links
    ~landlab.grid.radial.RadialModelGrid.number_of_fixed_links
    ~landlab.grid.radial.RadialModelGrid.number_of_links
    ~landlab.grid.radial.RadialModelGrid.number_of_links_at_node
    ~landlab.grid.radial.RadialModelGrid.resolve_values_on_active_links
    ~landlab.grid.radial.RadialModelGrid.resolve_values_on_links
    ~landlab.grid.radial.RadialModelGrid.status_at_link
    ~landlab.grid.radial.RadialModelGrid.unit_vector_xcomponent_at_link
    ~landlab.grid.radial.RadialModelGrid.unit_vector_ycomponent_at_link
    ~landlab.grid.radial.RadialModelGrid.upwind_links_at_node

Information about cells
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.radial.RadialModelGrid.area_of_cell
    ~landlab.grid.radial.RadialModelGrid.cell_area_at_node
    ~landlab.grid.radial.RadialModelGrid.cell_at_node
    ~landlab.grid.radial.RadialModelGrid.core_cells
    ~landlab.grid.radial.RadialModelGrid.faces_at_cell
    ~landlab.grid.radial.RadialModelGrid.node_at_cell
    ~landlab.grid.radial.RadialModelGrid.node_at_core_cell
    ~landlab.grid.radial.RadialModelGrid.number_of_cells
    ~landlab.grid.radial.RadialModelGrid.number_of_core_cells
    ~landlab.grid.radial.RadialModelGrid.number_of_faces_at_cell

Information about faces
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.radial.RadialModelGrid.active_faces
    ~landlab.grid.radial.RadialModelGrid.face_at_link
    ~landlab.grid.radial.RadialModelGrid.faces_at_cell
    ~landlab.grid.radial.RadialModelGrid.link_at_face
    ~landlab.grid.radial.RadialModelGrid.number_of_active_faces
    ~landlab.grid.radial.RadialModelGrid.number_of_faces
    ~landlab.grid.radial.RadialModelGrid.number_of_faces_at_cell
    ~landlab.grid.radial.RadialModelGrid.width_of_face

Information about patches
+++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.radial.RadialModelGrid.nodes_at_patch
    ~landlab.grid.radial.RadialModelGrid.number_of_patches
    ~landlab.grid.radial.RadialModelGrid.patches_at_node

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

    ~landlab.grid.radial.RadialModelGrid.at_node
    ~landlab.grid.radial.RadialModelGrid.at_cell
    ~landlab.grid.radial.RadialModelGrid.at_link
    ~landlab.grid.radial.RadialModelGrid.at_face

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

    ~landlab.grid.radial.RadialModelGrid.add_empty
    ~landlab.grid.radial.RadialModelGrid.add_field
    ~landlab.grid.radial.RadialModelGrid.add_ones
    ~landlab.grid.radial.RadialModelGrid.add_zeros
    ~landlab.grid.radial.RadialModelGrid.delete_field
    ~landlab.grid.radial.RadialModelGrid.set_units

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
    ~landlab.grid.radial.RadialModelGrid.field_units
    ~landlab.grid.radial.RadialModelGrid.field_values
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

    ~landlab.grid.radial.RadialModelGrid.calc_diff_at_link
    ~landlab.grid.radial.RadialModelGrid.calc_flux_div_at_node
    ~landlab.grid.radial.RadialModelGrid.calc_grad_at_link
    ~landlab.grid.radial.RadialModelGrid.calc_grad_at_patch
    ~landlab.grid.radial.RadialModelGrid.calc_net_flux_at_node
    ~landlab.grid.radial.RadialModelGrid.calc_slope_at_node
    ~landlab.grid.radial.RadialModelGrid.calc_slope_at_patch
    ~landlab.grid.radial.RadialModelGrid.calc_unit_normal_at_patch

Mappers
-------

These methods allow mapping of values defined on one grid element type onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

    ~landlab.grid.radial.RadialModelGrid.map_downwind_node_link_max_to_node
    ~landlab.grid.radial.RadialModelGrid.map_downwind_node_link_mean_to_node
    ~landlab.grid.radial.RadialModelGrid.map_link_head_node_to_link
    ~landlab.grid.radial.RadialModelGrid.map_link_tail_node_to_link
    ~landlab.grid.radial.RadialModelGrid.map_link_vector_to_nodes
    ~landlab.grid.radial.RadialModelGrid.map_max_of_link_nodes_to_link
    ~landlab.grid.radial.RadialModelGrid.map_max_of_node_links_to_node
    ~landlab.grid.radial.RadialModelGrid.map_mean_of_link_nodes_to_link
    ~landlab.grid.radial.RadialModelGrid.map_min_of_link_nodes_to_link
    ~landlab.grid.radial.RadialModelGrid.map_min_of_node_links_to_node
    ~landlab.grid.radial.RadialModelGrid.map_node_to_cell
    ~landlab.grid.radial.RadialModelGrid.map_upwind_node_link_max_to_node
    ~landlab.grid.radial.RadialModelGrid.map_upwind_node_link_mean_to_node
    ~landlab.grid.radial.RadialModelGrid.map_value_at_downwind_node_link_max_to_node
    ~landlab.grid.radial.RadialModelGrid.map_value_at_max_node_to_link
    ~landlab.grid.radial.RadialModelGrid.map_value_at_min_node_to_link
    ~landlab.grid.radial.RadialModelGrid.map_value_at_upwind_node_link_max_to_node

Boundary condition control
--------------------------

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to :meth:`~.ModelGrid.status_at_node` and
:meth:`~.ModelGrid.status_at_node` will automatically update the conditions
defined at other grid elements automatically.

.. autosummary::
    :toctree: generated/

   .....:     
    ~landlab.grid.radial.RadialModelGrid.active_faces
    ~landlab.grid.radial.RadialModelGrid.active_links
    ~landlab.grid.radial.RadialModelGrid.boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.closed_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.core_cells
    ~landlab.grid.radial.RadialModelGrid.core_nodes
    ~landlab.grid.radial.RadialModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.fixed_links
    ~landlab.grid.radial.RadialModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.node_at_core_cell
    ~landlab.grid.radial.RadialModelGrid.node_is_boundary
    ~landlab.grid.radial.RadialModelGrid.number_of_active_faces
    ~landlab.grid.radial.RadialModelGrid.number_of_active_links
    ~landlab.grid.radial.RadialModelGrid.number_of_core_cells
    ~landlab.grid.radial.RadialModelGrid.number_of_core_nodes
    ~landlab.grid.radial.RadialModelGrid.number_of_fixed_links
    ~landlab.grid.radial.RadialModelGrid.open_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.radial.RadialModelGrid.status_at_link
    ~landlab.grid.radial.RadialModelGrid.status_at_node

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

    ~landlab.grid.radial.RadialModelGrid.calc_aspect_at_node
    ~landlab.grid.radial.RadialModelGrid.calc_hillshade_at_node
    ~landlab.grid.radial.RadialModelGrid.calc_slope_at_node

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

from .voronoi import VoronoiDelaunayGrid
from landlab.utils.decorators import deprecated


class RadialModelGrid(VoronoiDelaunayGrid):

    """Grid of concentric circles.

    This inherited class implements a circular grid in which grid nodes are
    placed at regular radial and semi-regular arc-wise intervals. That is,
    if the radial spacing between *shells* is *dr*, the nodes are placed around
    the circular shell at regular intervals that get as close as possible to
    *dr*. The points are then arranged in a Delaunay triangulation with Voronoi
    cells. Within each ring, nodes are numbered according to Landlab
    convention, from the first node counterclockwise of east. Numbering
    begins at the centermost node and works outwards through the rings.

    Parameters
    ----------
    num_shells : int
        Number of rings in the grid.
    dr : float, optional
        Radial interval for rings.
    origin_x : float, optional
        x-coordinate of origin node.
    origin_y : float, optional
        y-coordinate of origin node.

    Returns
    -------
    RadialModelGrid
        A newly-created grid.

    Examples
    --------
    A grid with just one ring will have a node at the origin surrounded
    by six other nodes.

    >>> from landlab import RadialModelGrid
    >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0., origin_y=0.)
    >>> omg.number_of_nodes
    7
    >>> omg.number_of_cells
    1

    A second rings will have 13 nodes.

    >>> omg = RadialModelGrid(2)
    >>> omg.number_of_nodes
    20
    """

    def __init__(self, num_shells=0, dr=1.0, origin_x=0.0, origin_y=0.0,
                 **kwds):
        """Create a circular grid.

        Create a circular grid in which grid nodes are placed at regular
        radial and semi-regular arc-wise intervals. That is, if the radial
        spacing between *shells* is *dr*, the nodes are placed around the
        circular shell at regular intervals that get as close as possible to
        *dr*.  The points are then arranged in a Delaunay triangulation with
        Voronoi cells.

        Parameters
        ----------
        num_shells : int
            Number of rings in the grid.
        dr : float, optional
            Radial interval for rings.
        origin_x : float, optional
            x-coordinate of origin node.
        origin_y : float, optional
            y-coordinate of origin node.

        Returns
        -------
        RadialModelGrid
            A newly-created grid.

        Examples
        --------
        A grid with just one ring will have a node at the origin surrounded
        by six other nodes.

        >>> from landlab import RadialModelGrid
        >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0.,
        ...                       origin_y=0.)
        >>> omg.number_of_nodes
        7
        >>> omg.number_of_cells
        1

        A second rings will have 13 nodes.

        >>> omg = RadialModelGrid(2)
        >>> omg.number_of_nodes
        20
        """
        # Set number of nodes, and initialize if caller has given dimensions
        self._origin_x = origin_x
        self._origin_y = origin_y
        if num_shells > 0:
            self._initialize(num_shells, dr, origin_x, origin_y)
        super(RadialModelGrid, self).__init__(**kwds)

    @classmethod
    def from_dict(cls, params):
        num_shells = params['num_shells']
        dr = params.get('dr', 1.)
        origin = params.get('origin', (0., 0.))

        return cls(num_shells=num_shells, dr=dr, origin_x=origin[0],
                   origin_y=origin[1])

    def _initialize(self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        [pts, npts] = self._create_radial_points(num_shells, dr)
        self._n_shells = int(num_shells)
        self._dr = dr
        super(RadialModelGrid, self)._initialize(pts[:, 0], pts[:, 1])

    def _create_radial_points(self, num_shells, dr, origin_x=0.0,
                              origin_y=0.0):
        """Create a set of points on concentric circles.

        Creates and returns a set of (x,y) points placed in a series of
        concentric circles around the origin.
        """
        shells = numpy.arange(0, num_shells) + 1
        twopi = 2 * numpy.pi
        # number of points in each shell
        n_pts_in_shell = numpy.round(twopi * shells)
        dtheta = twopi / n_pts_in_shell
        npts = int(sum(n_pts_in_shell) + 1)
        pts = numpy.zeros((npts, 2))
        r = shells * dr
        startpt = 1
        for i in numpy.arange(0, num_shells):
            theta = (dtheta[i] * numpy.arange(0, n_pts_in_shell[i]) +
                     dtheta[i] / (i + 1))
            ycoord = r[i] * numpy.sin(theta)
            if numpy.isclose(ycoord[-1], 0.):
                # this modification necessary to force the first ring to
                # follow our new CCW from E numbering convention (DEJH, Nov15)
                ycoord[-1] = 0.
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    0] = numpy.roll(r[i] * numpy.cos(theta), 1)
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    1] = numpy.roll(ycoord, 1)
            else:
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    0] = r[i] * numpy.cos(theta)
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    1] = ycoord
            startpt += int(n_pts_in_shell[i])
        pts[:, 0] += origin_x
        pts[:, 1] += origin_y

        return pts, npts

    @property
    def number_of_shells(self):
        """Number of node shells in grid.

        Returns
        -------
        int
            The number of node shells in the radial grid (not counting the
            center node).
        """
        return self._n_shells

    @property
    @deprecated(use='spacing_of_shells', version=1.0)
    def shell_spacing(self):
        """Fixed distance between shells."""
        return self._dr

    @property
    def spacing_of_shells(self):
        """Fixed distance between shells."""
        return self._dr

    @property
    def number_of_nodes_in_shell(self):
        """Number of nodes in each shell.

        Returns
        -------
        int
            Number of nodes in each shell, excluding the center node.
        """
        try:
            return self._nnodes_inshell
        except AttributeError:
            n_pts_in_shell = numpy.round(2. * numpy.pi * (
                numpy.arange(self.number_of_shells, dtype=float) + 1.))
            self._nnodes_inshell = n_pts_in_shell.astype(int)
            return self._nnodes_inshell

    @property
    def radius_at_node(self):
        """Distance for center node to each node.

        Returns
        -------
        ndarray of float
            The distance from the center node of each node.

        >>> mg = RadialModelGrid(num_shells=2)
        >>> mg.radius_at_node
        array([ 2.,  2.,  2.,  2.,  2.,  1.,  1.,  2.,  0.,  1.,  1.,  2.,  2.,
                1.,  1.,  2.,  2.,  2.,  2.,  2.])
        """
        try:
            return self._node_radii
        except AttributeError:
            self._node_radii = numpy.sqrt(numpy.square(self.node_x -
                                                       self._origin_x) +
                                          numpy.square(self.node_y -
                                                       self._origin_y))
            return self._node_radii
