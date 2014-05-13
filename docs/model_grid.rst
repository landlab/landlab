========================================================================
Building Simple Models with Landlab's Gridding Library: A Tutorial Guide
========================================================================

(NOTE: currently, this web version of the Model Grid Guide is incomplete and not fully 
formatted. You can obtain a :download:`PDF version here <model_grid_guide/model_grid_description_and_guide.pdf>`.)

When creating a two-dimensional simulation model, often the most time-consuming and
error-prone task involves writing the code to set up the underlying grid. Irregular
(or ``unstructured'') grids are especially tricky to implement. Landlab's *ModelGrid*
package makes this process much easier, by providing a set of library routines for
creating and managing a 2D grid, attaching data to the grid, performing common input
and output operations, and  providing library functions that handle common numerical 
operations such as calculating a field of gradients for a particular state variable. 
By taking care of much of the overhead involved in writing grid-management code, 
*ModelGrid* is designed to help you build 2D models quickly and efficiently, freeing you
to concentration on the science behind the code.

Some of *ModelGrid's* capabilities include:

Landlab's *ModelGrid* is an open-source software package that creates and manages a regular
or irregular grid for building 2D numerical simulation models. ModelGrid is
especially useful for finite-volume (FV) and finite-difference (FD) models, but
also can be used for a variety of other applications. ModelGrid provides
efficient built-in functions for common operations in FD and FV models, such as
calculating local gradients and integrating fluxes around the perimeter of grid
cells. Staggered-grid models are especially easy to implement with ModelGrid.
A novel feature of ModelGrid is the ability to switch seamlessly between
structured and unstructured grids.

This document provides a basic introduction to building applications using
ModelGrid. It covers: (1) how grids are represented, (2) a tutorial example in
building a diffusion-model application, and (3) a guide to ModelGrid's methods
and data structures. ModelGrid is written in python. It is a component of the
Landlab modeling package.


How a Grid is Represented
=========================

.. _grid:

.. figure:: grid_schematic.png
    :scale: 50 %

    Elements of a model grid. Each grid comprises nodes, cells, faces, corners,
    patches, links, directed edges, and junctions. (Note that not all links,
    edges, and patches are shown, and only one representative cell is shaded.)


:ref:`Figure 1 <grid>` illustrates how ModelGrid represents a simulation grid. The
grid contains a set of *(x,y)* points called ``nodes``. In a typical
finite-difference model, nodes are the locations at which one tracks scalar
state variables, such as water depth, land elevation, or temperature. Each node
is associated with a polygon called a ``cell``. Each cell is bounded by a set
of line segments known as ``faces``, which it shares with its neighboring
cells.

In the simple case of a regular (raster) grid, the cells are square, the nodes
are the center points of the cells (:ref:`Figure 1a <grid>`), and the faces have
identical length (equal to the node spacing). In a Voronoi-Delaunay grid, the
cells are Voronoi polygons (also known as Theissen polygons)
(:ref:`Figure 1b <grid>`). In this case, each cell represents the surface area that
is closer to its own node than to any other node in the grid. The faces then
represent locations that are equidistant between two adjacent nodes. Other grid
configurations are possible as well. For examples, cells could be square
elements in a quad-tree grid (:ref:`Figure 1c <grid>`), or triangular elements with
nodes at their circumcenters (:ref:`Figure 1d <grid>`).

Each pair of adjacent cells is connected by a line segment called a ``link``
(:ref:`Figure 1, dashed line <grid>`). Each link connects a ``from node`` and a
``to node``, so it has direction as well as position and length. In some cases,
it may be useful to have each pair of adjacent cells connected by two vectors:
one pointing one way, and a second pointing the opposite way
\citep{guibas1985primitives,tucker2001object}. These vectors are known as
``directed edges`` (:ref:`Figure 1, gray arrows <grid>`). 

Finite-difference and finite-volume models usually need to calculate spatial
gradients in one or more scalar variables, and often these gradients are
evaluated between pairs of adjacent nodes. ModelGrid makes these calculations
easier for programmers by providing built-in functions to calculate gradients
along links, and allowing applications to associate an array of gradient values
with their corresponding links or edges.

The cell vertices are called ``corners`` (:ref:`Figure 1, solid squares <grid>`).
Each face is therefore a line segment connecting two corners. The intersection
of a face and a link (or directed edge) is known as a ``junction``
(:ref:`Figure 1, open diamonds <grid>`). Often, it is useful to calculate scalar
values (say, ice thickness in a glacier) at nodes, and vector values (say, ice
velocity) at junctions. This approach is sometimes referred to as a
staggered-grid scheme. It lends itself naturally to finite-volume methods, in
which one computes fluxes of mass, momentum, or energy across cell faces, and
maintains conservation of mass within cells
\citep[e.g.,][]{versteeg2007introduction}.

Notice that the links also enclose a set of polygons that are offset from the
cells. These secondary polygons are known as ``patches`` (:ref:`Figure 1,
dotted <grid>`). This means that any grid comprises two complementary tesselations: one
made of cells, and one made of patches. If one of these is a Voronoi
tessellation, the other is a Delaunay triangulation. For this reason, Delaunay
triangulations and Voronoi diagrams are said to be dual to one another: for any
given Delaunay triangulation, there is a unique corresponding Voronoi diagram
\citep[e.g.,][]{braun1997modelling,tucker2001object}. With ModelGrid, one can
create a mesh with either Voronoi polygons or Delaunay triangles as cells
(:ref:`Figure 1b and d <grid>`). Alternatively, with a raster grid, one simply has
two sets of square elements that are offset by half the grid spacing
(:ref:`Figure 1a <grid>`). Whatever the form of the tessellation, ModelGrid keeps
track of the geometry and topology of the grid. For example, one can call
various ModelGrid functions to obtain lists of the *(x,y)* coordinates of
nodes, corners, and junctions; get lists of neighbors for any cell; get the
endpoints of any link or directed edge, and so on. These functions are listed
and described below. 


How Boundaries are Managed
==========================

.. _raster4x5:

.. figure:: example_raster_grid.png
    :scale: 50 %

    Illustration of a simple four-row by five-column raster grid created with
    :class:`~landlab.grid.raster.RasterModelGrid`. By default, all perimeter
    nodes are tagged as active (fixed value) boundaries, and all interior cells
    are tagged as active interior. An active link is one that connects either
    two active interior cells, or one active interior and one active boundary.

.. _raster4x5openclosed:

.. figure:: example_raster_grid_with_closed_boundaries.png
    :scale: 50 %

    Illustration of a simple four-row by five-column raster grid with a
    combination of open and closed boundaries.

An important component of any numerical model is the method for handling
boundary conditions. In general, it's up to the application developer to manage
boundary conditions for each variable. However, ModelGrid makes this task a bit
easier by providing lists of nodes and links that lie along the boundary of the
grid, and those that lie in the interior. It also allows you to *de-activate*
portions of the grid perimeter, so that they effectively act as walls.

Let's look first at how ModelGrid treats its own geometrical boundaries. The
outermost elements of a grid are nodes and links (as opposed to corners and
faces). For example, :ref:`Figure 2 <raster4x5>` shows a sketch of a regular
four-row by five-column grid created by RasterModelGrid. The edges of the grid
are composed of nodes and links. Only the inner six nodes have cells around
them; the remaining 14 nodes form the perimeter of the grid.

All nodes are tagged as either *boundary* or *interior*. Those on the
perimeter of the grid are automatically tagged as boundary nodes. Nodes on the
inside are *interior* by default, but it is possible to tag some of them as
*boundary* instead (this would be useful, for example, if you wanted to
represent an irregular region inside a regular grid). In the example shown in
:ref:`Figure 2 <rester4x5>`, all the inner nodes are *active interior*, and all
perimeter nodes are *active boundary*. 

Boundary nodes are flagged as either *open* (active) or *closed*
(inactive), all links are tagged as *active* or *inactive*. An *active link*
is one that joins either two interior nodes, or an *interior* and an
*open boundary* node (:ref:`Figure 3 <raster4x5openclosed>`). You can use this
distinction in models to implement closed boundaries by performing flow
calculations only on active links, as the following simple example illustrates.

