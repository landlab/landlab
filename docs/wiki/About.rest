`Landlab <http://landlab.github.io>`_ |
[[ About | About ]] |
[[ Examples | Examples ]] |
[[ User Guide | User-Guide ]] |
`Reference Manual <http://landlab.readthedocs.org/en/latest/#developer-documentation>`_ |
[[ Tutorials| Tutorials ]] |
[[ FAQs | FAQs ]]

The Landlab project creates an environment in which scientists can build a
numerical surface process model without having to code all of the individual
components. Surface process models compute flows of mass, such as water, sediment,
glacial ice, volcanic material, or landslide debris, across a gridded terrain
surface. Surface process models have a number of commonalities, such as operating
on a grid of points and routing material across the grid. Scientists who want
to use a surface process model often build their own unique model from the ground
up, re-coding the basic building blocks of their surface process model rather than
taking advantage of codes that have already been written.

A list of papers and presentations using Landlab can be found on our [[ Landlab Papers and Presentations page | Landlab-Papers-and-Presentations ]].

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

If you are working with existing Landlab components, please also read the `component documentation <http://landlab.readthedocs.io/en/latest/#components>`_ to see if there are also specific publications linked to them. An example might be `Adams et al. (2017) <http://www.geosci-model-dev-discuss.net/gmd-2016-277/>`_ for the OverlandFlow components. Table 5 in Hobley et al. (2017) also lists many of these papers, as of the start of 2017.

A relatively new interface also automates the process of extracting citations for landlab. `import landlab` when you write your model, then a call to `landlab.registry.format_citations()` will list Bibtex-formatted citations for all Landlab components that you currently have instantiated.

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
