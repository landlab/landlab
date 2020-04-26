
.. _user_guide:

==========
User Guide
==========

The User Guide describes Landlab by topic area.

Users brand-new to Landlab should start with :ref:`10min`.

Further information on any specific method can be obtained in the
:ref:`api`.

The Landlab project creates an environment in which scientists can build a
numerical surface process model without having to code all of the individual
components. Surface process models compute flows of mass, such as water, sediment,
glacial ice, volcanic material, or landslide debris, across a gridded terrain
surface. Surface process models have a number of commonalities, such as operating
on a grid of points and routing material across the grid. Scientists who want
to use a surface process model often build their own unique model from the ground
up, re-coding the basic building blocks of their surface process model rather than
taking advantage of codes that have already been written.

A list of papers and presentations using Landlab can be found :ref:`here <papers>`.


Introduction to Python
----------------------

.. toctree::
    :maxdepth: 2

    python_numpy_intro


The Landlab Grid
----------------

.. toctree::
    :maxdepth: 2

    grid

Model with Landlab and Components
---------------------------------

.. toctree::
    :maxdepth: 2

    components
    build_a_model

Landlab and Units
-----------------

.. toctree::
    :maxdepth: 2

    units


Landlab Tutorial Library
------------------------

.. toctree::
    :maxdepth: 2

    tutorials
    teaching_tutorials

Additional resources
--------------------

.. toctree::
    :maxdepth: 2

    field_definitions
    field_io
    time_steps
    examples
    faq

Presentations, Clinics, and Classroom Use
-----------------------------------------

.. toctree::
    :maxdepth: 2

    papers_presentations
    clinics_workshops

Overland flow User Guide
------------------------

.. toctree::
    :maxdepth: 2

    overland_flow_user_guide

CellLab-CTS User Guide
----------------------

.. toctree::
    :maxdepth: 2

    cell_lab_user_guide


Major Version Transition Guides
-------------------------------

.. toctree::
    :maxdepth: 2

    landlab_zero_to_one
    landlab_one_to_two
    standard_name_changes
