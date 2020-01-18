.. _cite_as:

:orphan:

Landlab Citation Registry
-------------------------

A relatively new interface also automates the process of extracting citations
for Landlab components. Components with required software citations have an
attribute called ``cite_as``.

.. code-block:: python

    from landlab.components import OverlandFlow
    OverlandFlow.cite_as

This will give::

    @article{adams2017landlab,
           title={The Landlab v1.0 OverlandFlow component: a Python
                  tool for computing shallow-water flow across watersheds},
                  author={Adams, Jordan M and Gasparini, Nicole M and
                  Hobley, Daniel EJ and Tucker, Gregory E and
                  Hutton, Eric WH and Nudurupati, Sai S and
                  Istanbulluoglu, Erkan},
                  journal={Geoscientific Model Development},
                  volume={10},
                  number={4},
                  pages={1645},
                  year={2017},
                  publisher={Copernicus GmbH}
                  }

In addition, you can use the "citation registry" to make a .bib file of all
citations you've used in script.

.. code-block:: python

    import landlab
    w = landlab.registry.format_citations()
    with open("citations.bib", "w") as f:
        f.write(w)

This will produce a Bibtex-formatted citations for all Landlab components
that you currently have instantiated. This last example would have given the
main Landlab citation::

    # Citations
    ## landlab
    @article{hobley2017creative,
    AUTHOR = {
    Hobley, D. E. J. and Adams, J. M. and
    Nudurupati, S. S. and Hutton,
    E. W. H. and Gasparini, N. M. and
    Istanbulluoglu, E. and Tucker, G. E.
    },
    TITLE = {
    Creative computing with Landlab: an
    open-source toolkit for building,
    coupling, and exploring two-dimensional
    numerical models of Earth-surface
    dynamics
    },
    JOURNAL = {Earth Surface Dynamics},
    VOLUME = {5},
    YEAR = {2017},
    NUMBER = {1},
    PAGES = {21--46},
    URL = {
    https://www.earth-surf-dynam.net/5/21/2017/
    },
    DOI = {10.5194/esurf-5-21-2017}
    }
