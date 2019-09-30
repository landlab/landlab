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
across scales.

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
--------------

`Hobley, D. E. J. <http://www.earth-surf-dynam.net/5/21/2017/>`_, Adams, J. M., Nudurupati, S. S., Hutton, E. W. H., Gasparini, N. M., Istanbulluoglu, E. and Tucker, G. E., 2017, Creative computing with Landlab: an open-source toolkit for building, coupling, and exploring two-dimensional numerical models of Earth-surface dynamics, Earth Surface Dynamics, 5, p 21-46, 10.5194/esurf-5-21-2017.

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

If you are working with existing Landlab components, please also read the
:ref:`component documentation <api.components>`
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

.. _contact:

Contact
=======

Questions? Feedbacks?

Need an improvement/addition to Landlab?

Want to contribute?

You can contact the Landlab development team here:

-  Make an `Issue <https://github.com/landlab/landlab/issues>`_

To keep in touch with the latest Landlab news:

-  Get the `Landlab Lookout Newsletter <http://eepurl.com/dADtrT>`_
-  Join the `Landlab Slack channel <landlab.slack.com>`_
-  Landlab is on `Twitter <https://twitter.com/landlabtoolkit>`_!
