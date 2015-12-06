What is Landlab?
================

.. image:: images/landlab_logo.jpg
   :scale: 70 %

(originally authored by Greg Tucker, December 2013)

Landlab is a Python software package that supports numerical modeling in earth science, 
and especially those fields that deal with earth-surface dynamics, including geomorphology, 
hydrology, glaciology, stratigraphy, and related areas. Landlab is a modeling environment in
which scientist can build a numerical landscape model without having to code all of the individual
parts.  Landscape models compute flows of mass, such as water, sediment,
glacial ice, volcanic material, or landslide debris, across a gridded terrain
surface. Landscape models have a number of commonalities, such as operating
on a grid of points and routing material across the grid. Scientists who want
to use a landscape model often build their own unique model from the ground
up, re-coding the basic building blocks of their landscape model rather than
taking advantage of codes that have already been written.

Landlab provides four main resources and capabilities:

(1) A library of *code resources* for building two-dimensional numerical models from scratch. The Landlab library includes a powerful "gridding engine" for creating, managing, and iteratively updating data on structured or unstructured grids. The library also includes support for input and output, including input of digital elevation models (DEMs) in ArcInfo ASCII format, handling of parameter inputs using formatted text files, and netCDF-format input and output.

(2) A set of pre-built *components*, each of which implements a numerical representation of a particular process.

(3) A *framework for building models* by assembling and linking process components.

(4) A *library of models* that have already been created from component(s).

An example of some models that can easily be built with Landlabs current capabilities are:

* A landscape evolution model using linear diffusion and the stream-power model
* A model that explores the flexural response to the growth and recession of glaciers
* An ecohydrology model in which vegetation on two sides of a valley grows and dies in response to stochastic storms and solar forcing throughout the year
* A model that routes hydrographs across a watershed based on rainfall inputs

Acknowledgements
----------------

The Landlab Team:
  - Greg Tucker (University of Colorado)
  - Nicole Gasparini (Tulane University)
  - Erkan Istanbulluoglu (University of Washington)
  - Daniel Hobley (CU)
  - Sai Nudurupati (UW)
  - Jordan Adams (TU)
  - Eric Hutton (CU)

Initial funding was provided by a grant from the National Science Foundation to Greg Tucker (ACI 1147454), 
Erkan Istanbulluoglu (ACI 1148305), and Nicole Gasparini (ACI 1147519).
