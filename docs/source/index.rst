.. _landlab:

.. module:: landlab

==========================================================
Landlab: A modular Earth Surface Dynamics modeling library
==========================================================


Landlab is a Python-based library that allows scientists and students to
build numerical landscape models. Designed for disciplines that quantify
earth surface dynamics such as geomorphology, hydrology, glaciology, and
stratigraphy, it can also be used in related fields.

Landlab provides components to compute flows (such as water, sediment,
glacial ice, volcanic material, or landslide debris) across a gridded
terrain. With its robust, reusable components, Landlab allows scientists
to quickly build landscape model experiments and compute mass balance
across scales. Watch the webinar `Landlab Toolkit Overview <https://csdms.colorado.edu/wiki/Presenters-0407>`_ at CSDMS to learn more.

The most current source code is always available from our `git
repository <https://github.com/landlab/landlab>`_.

If you are interested in the state of current development, we compile
:ref:`ongoing development <ongoing_development>`. If you are interested in
contributing but don't know how to get started (THANK YOU!),
we compile :ref:`desired contributions <desired_contributions>` and have a
:ref:`page all about development <development>`.

Many Jupyter notebooks exist describing use of Landlab. Find an overview
:ref:`here <tutorials>`. A subset of these notebooks are designed specifically
for the classroom. Information about them and how to set them up for classroom
use is described :ref:`on this page <teaching_tutorials>`.


.. toctree::
   :maxdepth: 2

   install/index
   user_guide/index
   reference/index
   whatsnew/index
   getting_started/index
   development/index

Acknowledgements
================

Funding
-------

Landlab is funded by the US National Science Foundation. It has been supported
by the following grants:

    * A Collaborative NSF SI2-SSE proposal to
      University of Colorado (Greg Tucker,
      `1147454 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1147454&HistoricalAwards=false>`_),
      and the University of Washington (Erkan Istanbulluoglu,
      `1148305 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1148305&HistoricalAwards=false>`_)
    * A Collaborative NSF SI2-SSI proposal to
      University of Colorado (Greg Tucker and Dan Hobley,
      `1450409 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1450409&HistoricalAwards=false>`_),
      Tulane University (Nicole Gasparini,
      `1450338 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1450338&HistoricalAwards=false>`_),
      and the University of Washington (Erkan Istanbulluoglu,
      `1450412 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1450412&HistoricalAwards=false>`_).
    * A NSF EAR Postdoctoral Fellowship to Katy Barnhart
      (`1725774 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1725774&HistoricalAwards=false>`_).
    * TODO: Add CSDMS grants, add any other grant.

Citing Landlab
==============

If you use any portion of Landlab, you must cite the following paper:

`Hobley, D. E. J. <https://www.earth-surf-dynam.net/5/21/2017/>`_, Adams, J. M., Nudurupati, S. S., Hutton, E. W. H., Gasparini, N. M., Istanbulluoglu, E. and Tucker, G. E., 2017, Creative computing with Landlab: an open-source toolkit for building, coupling, and exploring two-dimensional numerical models of Earth-surface dynamics, Earth Surface Dynamics, 5, p 21-46, 10.5194/esurf-5-21-2017.

In addition, many components have an additional citation. The table below lists
these. We also provide a handy way to get a bibtex file of the components you
use. It is described below.

BibTeX format:
::

 @article{Hobley2017,
          Author = {Hobley, D. E. J. and Adams, J. M. and
                    Nudurupati, S. S. and Hutton, E. W. H.
                    and Gasparini, N. M. and Istanbulluoglu,
                    E. and Tucker, G. E.},
          Journal = {Earth Surface Dynamics},
          Year = {2017},
          Title = {Creative computing with Landlab: an open-source
                   toolkit for building, coupling, and exploring
                   two-dimensional numerical models of
                   Earth-surface dynamics},
          Number = {5},
          Pages = {21-46},
          Doi = {10.5194/esurf-5-21-2017}}

If you are working with Landlab components and utilities, many of them have
their own publication. Please cite it to acknowledge the component authors.

The following is a summary table of component and utility citations:

TODO: we need to add to citations from Table 5 in Hobley et al. (2017) to this
table.

+----------------------------------+-----------------------------------------------------------------------------------------------------+
| Component                        | Citation                                                                                            |
+==================================+=====================================================================================================+
| Cellular Automaton               | `Tucker et al. (2015) <https://www.geosci-model-dev.net/9/823/2016/>`_                              |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| Flexure                          | `Hutton and Syvitski (2008) <https://www.sciencedirect.com/science/article/pii/S0098300408000587>`_ |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| LakeMapperBarnes                 | `Barnes et al. (2014) <https://linkinghub.elsevier.com/retrieve/pii/S0098300413001337>`_            |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| LandslideProbability             | `Strauch et al. (2017) <https://www.earth-surf-dynam.net/6/49/2018/esurf-6-49-2018.html>`_          |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| LateralEroder                    | `Langston and Tucker (2018) <https://www.earth-surf-dynam.net/6/1/2018/>`_                          |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| Lithology                        | `Barnhart et al. (2018) <https://joss.theoj.org/papers/10.21105/joss.00979>`_                       |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| LithoLayers                      | `Barnhart et al. (2018) <https://joss.theoj.org/papers/10.21105/joss.00979>`_                       |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| OverlandFlow                     | `Adams et al. (2017) <https://www.geosci-model-dev-discuss.net/gmd-2016-277/>`_                     |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| SPACE                            | `Shobe et al. (2017) <https://www.geosci-model-dev.net/10/4577/2017/>`_                             |
+----------------------------------+-----------------------------------------------------------------------------------------------------+
| SpatialPrecipitationDistribution | `Singer et al. (2018) <https://www.geosci-model-dev.net/11/3713/2018/gmd-11-3713-2018.html>`_       |
+----------------------------------+-----------------------------------------------------------------------------------------------------+



A relatively new interface also automates the process of extracting citations
for landlab.

.. code-block:: python

  import landlab

  # the code of your model built with landlab goes here, then add a call to

  landlab.registry.format_citations() # will produce a Bibtex-formatted
                                      # citations for all Landlab components
                                      # that you currently have instantiated.


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

.. _contact:

Contact
=======

Questions? Feedbacks?

Need an improvement/addition to Landlab?

Want to contribute?

The recommended way to contact the Landlab development team is with a
`GitHub Issue <https://github.com/landlab/landlab/issues>`_

To keep in touch with the latest Landlab news:

-  Get the `Landlab Lookout Newsletter <https://github.us18.list-manage.com/subscribe?u=2db7cea82e3ea40fcf4c91247&id=b9bad233c7>`_
-  Landlab is on `Twitter <https://mobile.twitter.com/landlabtoolkit>`_!

During workshops and clinics, we often use the
`Landlab Slack channel <https://landlab.slack.com>`_
