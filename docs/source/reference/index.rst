.. _api:

=============
API reference
=============

This page gives an overview of all public Landlab objects, functions and
methods.

Grids
-----

.. toctree::
    :maxdepth: 2

    grid/index

Layers
------

.. toctree::
    :maxdepth: 2

    layers/index

Lithology
---------

Two objects based on the EventLayers object exist to make it easier to deal
with spatially variable lithology and associated properties. The Lithology
components contain information about spatially variable lithology and connect
with the Landlab model grid so that when rock is eroded or advected upward by
rock uplift the values of rock propeties at the topographic surface are updated.

First is the Lithology component which is a generic object for variable
lithology.

.. toctree::
    :maxdepth: 1

    /generated/api/landlab.components.lithology.lithology

Second is LithoLayers which makes it easy to make layered rock.

.. toctree::
    :maxdepth: 1

    /generated/api/landlab.components.lithology.litholayers

Values
------

.. toctree::
    :maxdepth: 2

    values/index

Components
----------

.. toctree::
    :maxdepth: 2

    components/index

..
    .. toctree::
       :maxdepth: 2

       /generated/api/landlab.bmi
       /generated/api/landlab.ca
       /generated/api/landlab.cmd
       components/index
       /generated/api/landlab.core
       /generated/api/landlab.data_record
       /generated/api/landlab.field
       /generated/api/landlab.framework
       /generated/api/landlab.graph
       grid/index
       /generated/api/landlab.io
       layers/index
       /generated/api/landlab.plot
       /generated/api/landlab.utils
       values/index

Full Index
----------

.. toctree::
   :maxdepth: 2

   modules

References
----------

* :ref:`modindex`
* :ref:`search`


.. Search the Index
.. ==================

.. * :ref:`genindex`
