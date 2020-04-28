.. _landlab:

.. module:: landlab

==========================================================
Landlab: A modular Earth Surface Dynamics modeling library
==========================================================

Landlab is an open-source Python-language package for numerical modeling of
Earth surface dynamics. It contains:

* A gridding engine which represents the model domain. Regular and irregular
  grids are supported.
* A library of process components, each of which represents a physical process
  (e.g., generation of rain, erosion by flowing water). These components have
  a common interface and can be combined based on a user's needs.
* Utilities that support general numerical methods, file input/output, and
  visualization.

In addition Landlab contains a set of Jupyter notebook tutorials providing
an introduction to core concepts and examples of use.

Landlab was designed for disciplines that quantify Earth surface dynamics such
as geomorphology, hydrology, glaciology, and stratigraphy. It can also be used
in related fields. Scientists who use this type of model often build
their own unique model from the ground up, re-coding the basic building blocks
of their landscape model rather than taking advantage of codes that have
already been written. Landlab saves practitioners from the need for this kind
of re-invention by providing standardized components that they can re-use.

Watch the webinar `Landlab Toolkit Overview <https://csdms.colorado.edu/wiki/Presenters-0407>`_
at CSDMS to learn more.

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

Landlab 2.0
-----------

In late December 2019 Landlab switched to version 2.0-beta. Landlab will be
in 2.0-beta until the Landlab 2.0 publication is finalized. Landlab dropped
support of Python 2.7 with this transition.

Supported Python Versions
-------------------------

Landlab supports Python versions >= 3.6. Landlab distributes pre-packaged
binaries through `conda-forge <https://anaconda.org/conda-forge/landlab>`_
and `PyPI <https://pypi.org/project/landlab/>`_ for versions 3.6 and 3.7
(3.8 coming soon).

Documentation Outline
---------------------

.. toctree::
   :maxdepth: 2

   install/index
   getting_started/index
   user_guide/index
   reference/index
   development/index
   whatsnew/index

Contact
-------

Questions? Feedback? Found a bug or something unexpected?

Need an improvement/addition to Landlab?

Want to contribute?

The recommended way to contact the Landlab development team is with a
`GitHub Issue <https://github.com/landlab/landlab/issues>`_.

Landlab News
------------

To keep in touch with the latest Landlab news:

-  Get the `Landlab Lookout Newsletter <https://github.us18.list-manage.com/subscribe?u=2db7cea82e3ea40fcf4c91247&id=b9bad233c7>`_
-  Landlab is on `Twitter <https://twitter.com/landlabtoolkit>`_!

Citing Landlab
--------------

If you use any portion of Landlab, you must cite the following papers:

`Hobley, D. E. J. <https://www.earth-surf-dynam.net/5/21/2017/>`__, Adams,
J. M., Nudurupati, S. S., Hutton, E. W. H., Gasparini, N. M., Istanbulluoglu,
E. and Tucker, G. E., 2017, Creative computing with Landlab: an open-source
toolkit for building, coupling, and exploring two-dimensional numerical models
of Earth-surface dynamics, Earth Surface Dynamics, 5, p 21-46,
10.5194/esurf-5-21-2017.

`Barnhart, K. R. <https://doi.org/10.5194/esurf-2020-12>`__,
Hutton, E. W. H., Tucker, G. E., Gasparini, N. M., Istanbulluoglu, E.,
Hobley, D. E. J., Lyons, N. J., Mouchene, M., Nudurupati, S. S., Adams, J. M.,
and Bandaragoda, C.: Short communication: Landlab v2.0: A software package for
Earth surface dynamics, Earth Surf. Dynam. Discuss.,
https://doi.org/10.5194/esurf-2020-12, in review, 2020.


BibTeX format:
::

 @article{Hobley2017,
          author = {Hobley, D. E. J. and Adams, J. M. and
                    Nudurupati, S. S. and Hutton, E. W. H.
                    and Gasparini, N. M. and Istanbulluoglu,
                    E. and Tucker, G. E.},
          journal = {Earth Surface Dynamics},
          year = {2017},
          title = {Creative computing with Landlab: an open-source
                   toolkit for building, coupling, and exploring
                   two-dimensional numerical models of
                   Earth-surface dynamics},
          number = {5},
          pages = {21-46},
          doi = {10.5194/esurf-5-21-2017}}

 @article{barnhart2020short,
          author = {Barnhart, K. R. and Hutton, E. W. H. and
                    Tucker, G. E. and Gasparini, N. M. and
                    Istanbulluoglu, E. and Hobley, D. E. J. and
                    Lyons, N. J. and Mouchene, M. and Nudurupati,
                    S. S. and Adams, J. M. and Bandaragoda, C.},
          title = {Short communication: Landlab v2.0: A software
                   package for Earth surface dynamics},
          journal = {Earth Surface Dynamics Discussions},
          volume = {2020},
          year = {2020},
          pages = {1--25},
          url = {https://www.earth-surf-dynam-discuss.net/esurf-2020-12/},
          doi = {10.5194/esurf-2020-12}
          }

If you are working with Landlab components and utilities, many of them have
their own publication. Please cite it to acknowledge the component authors.

Citation information for each component can be found as follows:

- Where relevant, software citation and general references, are listed in the
  Component API documentation under the References section.
- Software citations are included in component metadata. We have created a
  tool called the "Citation Registry" which creates a .bib file for software
  citations used in an application. See example usage :ref:`here <cite_as>`.

.. _contact:

The Landlab Team
----------------

The core development team is

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

Workshops and Clinics
---------------------

During :ref:`workshops and clinics <clinics_workshops>`, we sometimes use the
`Landlab Slack channel <https://landlab.slack.com>`_.

Funding
-------

Landlab is funded by the US National Science Foundation. It has been supported
by the following grants:

   * A Collaborative NSF SI2-SSE proposal to
     University of Colorado (Greg Tucker,
     `1147454 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1147454&HistoricalAwards=false>`_),
     Tulane University (Nicole Gasparini,
     `1147519 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1147519&HistoricalAwards=false>`_),
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
