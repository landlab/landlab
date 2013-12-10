================
Model Grid Guide
================

(NOTE: currently, this web version of the Model Grid Guide is incomplete and not fully 
formatted. You can obtain a :download:`PDF version here <model_grid_guide/model_grid_description_and_guide.pdf>`.)

ModelGrid is an open-source software package that creates and manages a regular
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


