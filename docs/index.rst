.. landlab documentation master file, created by
   sphinx-quickstart on Tue Apr 23 23:40:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. Introduction to Landlab
.. =======================

.. toctree::
..   :maxdepth: 1

..   what_is_landlab

.. toctree::
..   :maxdepth: 2

..   install
..   getting_started

.. toctree::
..   :maxdepth: 1

..   dan_installs_on_linux
..   getting_example_files

Landlab Developer API
==============================

Find Landlab's `User Guide <https://github.com/landlab/landlab/wiki/User-Guide>`_ on the `Landlab Wiki <https://github.com/landlab/landlab/wiki/User-Guide>`_

Contributing to Landlab
=======================

If you're intending to make changes to the Landlab code base,
or want to develop your own components, we recommend you follow
these specialized developer install instructions.

.. toctree::
  :maxdepth: 3

  dev_guide_install
  dev_guide_components


.. toctree::
..    :maxdepth: 1

.. The Nuts and Bolts of Coding in Landlab
.. ---------------------------------------
.. toctree::
   :maxdepth: 2

..   coding_style


.. Landlab's Gridding Library
.. --------------------------

.. toctree::
..   :maxdepth: 2

..   model_grid_no_tutorials


.. Building a Model
.. ----------------

.. toctree::
..   :maxdepth: 3

..   working_with_landlab

.. toctree::
..   :maxdepth: 3

..   landlab_components


.. Landlab Grid Data Structures
.. ----------------------------
..
.. Quick links to the landlab grid data structures:
..
.. for some reason, these weren't working
.. * `RasterModelGrid <landlab.readthedocs.org/en/latest/manual_index_alt_format.html#landlab.grid.raster.RasterModelGrid.__init__>`_
.. * `HexModelGrid <landlab.readthedocs.org/en/latest/manual_index_alt_format.html#landlab.grid.raster.HexModelGrid.__init__>`_
.. * `RadialModelGrid <landlab.readthedocs.org/en/latest/manual_index_alt_format.html#landlab.grid.raster.RadialModelGrid.__init__>`_
.. * `VoronoiDelaunayGrid <landlab.readthedocs.org/en/latest/manual_index_alt_format.html#landlab.grid.raster.VonoroiDelaunayGrid.__init__>`_
..
..


.. Tutorials
.. ---------

.. These tutorials are NOT IN MODERN STYLE, as of 05/25/15
.. Thus DEJH has commented them out
.. .. toctree::
..    :maxdepth: 1
.. replaced with the notebook tutorials
.. diffusion_raster_grid_tutorial
.. overland_flow_general_tutorial
.. overland_flow_dem_tutorial

.. Start with the :ref:`10 minute Landlab introduction tutorial <getting_started>`, then choose from:

.. * A super-basic intro to Python and Numpy: http://nbviewer.ipython.org/github/landlab/tutorials/blob/master/python_intro/Python_intro.ipynb
.. * An introduction to modelling with Landlab: http://nbviewer.ipython.org/github/landlab/tutorials/blob/master/fault_scarp_notebook/landlab-fault-scarp.ipynb
.. * Using the Landlab component library: http://nbviewer.ipython.org/github/landlab/tutorials/blob/master/component_tutorial/component_tutorial.ipynb
.. * The Landlab flexure component: http://nbviewer.ipython.org/github/landlab/tutorials/blob/master/flexure/lots_of_loads.ipynb
.. * The Landlab ecohydrology components: http://nbviewer.ipython.org/github/landlab/tutorials/blob/master/ecohydrology/cellular_automaton_vegetation_DEM/cellular_automaton_vegetation_DEM.ipynb


.. CellLab-CTS
.. -----------

.. CellLab-CTS is a Landlab module for building pairwise, continuous-time stochastic (CTS) cellular automata.

.. toctree::
	:maxdepth: 1

..	celllab_manual.rst


.. Frequently Asked Questions
.. ==========================

.. toctree::
..   :maxdepth: 1

..   faq


Grids
=======================

As of Landlab version 0.2, there are four types of Landlab grid:
 - `Raster <http://landlab.readthedocs.org/en/master/landlab.grid.html#module-landlab.grid.raster>`_
 - `Voronoi-DeLaunay <http://landlab.readthedocs.org/en/master/landlab.grid.html#module-landlab.grid.voronoi>`_
 - `Hex <http://landlab.readthedocs.org/en/master/landlab.grid.html#module-landlab.grid.hex>`_
 - `Radial <http://landlab.readthedocs.org/en/master/landlab.grid.html#module-landlab.grid.radial>`_

The base class is `ModelGrid` with subclasses `RasterModelGrid` and `VoronoiDelaunayGrid`.

`VoronoiDelaunayGrid` has two further specialized subclasses: `HexModelGrid` and `RadialModelGrid`.

.. toctree::
   :maxdepth: 4

   landlab.grid


Components
=======================

.. toctree::
   :maxdepth: 4

   landlab.components


Input/Output (IO)
=======================

.. toctree::
   :maxdepth: 4

   landlab.io


Plotting and Visualization
=======================

.. toctree::
   :maxdepth: 4

   landlab.plot


Cellular Automata (CA)
=======================

.. toctree::
   :maxdepth: 4

   landlab.ca



.. Simple guides to functionality
.. ------------------------------

.. These (slightly outdated) resources provide guides to the actual functions you can find and use through Landlab.

.. toctree::
..  :maxdepth: 1

..  users_guide

.. toctree::
..  :maxdepth: 2

..  manual_index_alt_format

.. toctree::
..  :maxdepth: 1

..  standard_names

References
==========

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
