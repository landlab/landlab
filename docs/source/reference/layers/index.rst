.. _api.layers:

======
Layers
======

Landlab has the ability to add layers to the grid. Two types of layers are
currently supported. First is EventLayers in which each event is preserved as
an entry into the datastructure, even if no deposition occurs. If you are
interested in chronostratigraphy, this is probably what you are interested in.
Second is MaterialLayers, in which each layer must contain some material.
If an entire layer is eroded in MaterialLayers, the layer is removed.
MaterialLayers will likely use less memory than EventLayers.

  .. toctree::
     :maxdepth: 2

     eventlayers
     materiallayers


Extensions
----------

.. toctree::
   :maxdepth: 2

   layers_ext


Module contents
---------------

.. automodule:: landlab.layers
  :members:
  :undoc-members:
  :show-inheritance:
