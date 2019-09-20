
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

A list of papers and presentations using Landlab can be found on here:

.. toctree::
    :maxdepth: 2

    papers_presentations


Acknowledgements
----------------

Citing Landlab:

`Hobley, D. E. J. <http://www.earth-surf-dynam.net/5/21/2017/>`_, Adams, J. M., Nudurupati, S. S., Hutton, E. W. H., Gasparini, N. M., Istanbulluoglu, E. and Tucker, G. E., 2017, Creative computing with Landlab: an open-source toolkit for building, coupling, and exploring two-dimensional numerical models of Earth-surface dynamics, Earth Surface Dynamics, 5, p 21-46, 10.5194/esurf-5-21-2017.

BibTeX format:
::

  @article{Hobley2017,
                Author = {Hobley, D. E. J. and Adams, J. M. and Nudurupati, S. S. and Hutton, E. W. H. and Gasparini, N. M. and Istanbulluoglu, E. and Tucker, G. E.},
                Journal = {Earth Surface Dynamics},
                Year = {2017},
                Title = {Creative computing with Landlab: an open-source toolkit for building, coupling, and exploring two-dimensional numerical models of Earth-surface dynamics},
                Number = {5},
                Pages = {21-46},
                Doi = {10.5194/esurf-5-21-2017}}

If you are working with existing Landlab components, please also read the
`component documentation <http://landlab.readthedocs.io/en/latest/#components>`_
to see if there are also specific publications linked to them. An example might be
`Adams et al. (2017) <http://www.geosci-model-dev-discuss.net/gmd-2016-277/>`_
for the OverlandFlow components. Table 5 in Hobley et al. (2017) also lists
many of these papers, as of the start of 2017.

A relatively new interface also automates the process of extracting citations
for landlab. `import landlab` when you write your model, then a call to
`landlab.registry.format_citations()` will list Bibtex-formatted citations for
all Landlab components that you currently have instantiated.

The Landlab Team:
  - Greg Tucker (CU)
  - Nicole Gasparini (Tulane)
  - Erkan Istanbulluoglu (UW)
  - Daniel Hobley (Cardiff)
  - Sai S. Nudurupati (UW)
  - Jordan Adams (Tulane)
  - Eric Hutton (CU)
  - Jenny Knuth (CU)
  - Katy Barnhart (CU)
  - Margaux Mouchene (CU)
  - Christina Bandaragoda (UW)
  - Nathan Lyons (Tulane)


Installation
------------

.. toctree::
    :maxdepth: 2

    ../install/index


Updating and Troubleshooting
----------------------------

.. toctree::
    :maxdepth: 2

    updating
    troubleshooting


Introduction to Python
----------------------

.. toctree::
    :maxdepth: 2

    python_numpy_intro
    rough_guide_to_python


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


Additional resources
--------------------

.. toctree::
    :maxdepth: 2

    tutorials
    standard_names
    time_steps
    faq
    contact
    examples

Presentations, Clinics, and Classroom Use
-----------------------------------------

.. toctree::
    :maxdepth: 2

    papers_presentations
    clinic_workshops
    landlab_curriculum


Module User Guides
------------------

.. toctree::
    :maxdepth: 2

    overland_flow_user_guide
    cell_lab_user_guide


Major Version Transition Guides
-------------------------------

.. toctree::
    :maxdepth: 2

    landlab_zero_to_one


Landlab User Group
------------------

To join the Landlab User Group on Slack, send an email request to: Join
Landlab User Group
