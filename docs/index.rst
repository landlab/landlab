.. landlab documentation master file, created by
   sphinx-quickstart on Tue Apr 23 23:40:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Find Landlab's `User Guide <https://github.com/landlab/landlab/wiki/User-Guide>`_ on the `Landlab Wiki <https://github.com/landlab/landlab/wiki/User-Guide>`_

==============================
Landlab Developer API
==============================

.. Contents::


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


Contributing to Landlab
=======================

If you're intending to make changes to the Landlab code base,
or want to develop your own components, we recommend you follow
these specialized developer install instructions.

.. toctree::
 :maxdepth: 3

 dev_guide_install
 dev_guide_components


References
==========

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
