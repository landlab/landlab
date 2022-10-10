
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


.. toctree::
    :caption: Introduction to Python
    :maxdepth: 2
    :hidden:

    python_numpy_intro


.. toctree::
    :caption: The Landlab Grid
    :maxdepth: 2
    :hidden:

    grid


.. toctree::
    :caption: Model with Landlab and Components
    :maxdepth: 2
    :hidden:

    components
    build_a_model


.. toctree::
    :caption: Landlab and Units
    :maxdepth: 2
    :hidden:

    units


.. toctree::
    :caption: Landlab Tutorial Library
    :maxdepth: 2
    :hidden:

    tutorials
    teaching_tutorials


.. toctree::
    :caption: Additional resources
    :maxdepth: 2
    :hidden:

    component_list
    field_definitions
    grid_summary
    time_steps
    examples
    faq


.. toctree::
    :caption: Overland flow User Guide
    :maxdepth: 2
    :hidden:

    overland_flow_user_guide


.. toctree::
    :caption: CellLab-CTS User Guide
    :maxdepth: 2
    :hidden:

    cell_lab_user_guide


.. toctree::
    :caption: Major Version Transition Guides
    :maxdepth: 2
    :hidden:

    landlab_zero_to_one
    landlab_one_to_two
    standard_name_changes
