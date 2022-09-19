.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3776837.svg
   :target: https://doi.org/10.5281/zenodo.3776837

.. image:: https://readthedocs.org/projects/landlab/badge/?version=latest
    :target: https://landlab.readthedocs.org

.. image:: https://github.com/landlab/landlab/actions/workflows/test.yml/badge.svg
    :target: https://github.com/landlab/landlab/actions/workflows/test.yml

.. image:: https://github.com/landlab/landlab/actions/workflows/flake8.yml/badge.svg
    :target: https://github.com/landlab/landlab/actions/workflows/flake8.yml

.. image:: https://github.com/landlab/landlab/actions/workflows/black.yml/badge.svg
    :target: https://github.com/landlab/landlab/actions/workflows/black.yml

.. image:: https://github.com/landlab/landlab/actions/workflows/docs.yml/badge.svg
    :target: https://github.com/landlab/landlab/actions/workflows/docs.yml

.. image:: https://coveralls.io/repos/landlab/landlab/badge.png
    :target: https://coveralls.io/r/landlab/landlab

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb

=======
Landlab
=======

What does Landlab do?
---------------------

.. start-intro

Landlab is an open-source Python-language package for numerical modeling of
Earth surface dynamics. It contains

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

.. end-intro

-----------

`Read the documentation on ReadTheDocs! <https://landlab.readthedocs.io/>`_

-----------

Installation
------------

.. start-install-release

In order to use *landlab* you will first need Python. While not
necessary, we recommend using the 
`Anaconda Python distribution <https://www.anaconda.com/distribution/>`_
as it provides a large number of third-party packages useful for
scientific computing.

To install *landlab*, simply run the following in your terminal of choice:

.. tab:: mamba

  .. code-block:: bash

    $ conda install mamba -c conda-forge
    $ mamba install landlab -c conda-forge

.. tab:: conda

  .. code-block:: bash

    $ conda install landlab -c conda-forge

.. tab:: pip

  .. code-block:: bash

    $ pip install landlab

.. end-install-release

Source code
-----------

.. start-install-source

*landlab* is actively being developed on GitHub, where the code is freely available.
If you would like to modify or contribute code, you can either clone our
repository

.. code-block:: bash

   $ git clone git://github.com/landlab/landlab.git

or download the `tarball <https://github.com/landlab/landlab/tarball/master>`_
(a zip file is available for Windows users):

.. code-block:: bash

   $ curl -OL https://github.com/landlab/landlab/tarball/master

Once you have a copy of the source code, you can install it into your current
Python environment,

.. tab:: conda

  .. code-block:: bash

     $ cd landlab
     $ conda install --file=requirements.txt
     $ pip install -e .

.. tab:: pip

  .. code-block:: bash

     $ cd landlab
     $ pip install -e .


.. end-install-source



Are there any examples of using Landlab I can look at?
------------------------------------------------------

The Landlab package contains a directory, ``landlab/notebooks``, with
Jupyter Notebooks describing core concepts and giving examples of using components.
The file ``landlab/notebooks/welcome.ipynb`` provides a table of contents to
the notebooks and is the recommended starting place.
Additionally, there are a set of notebooks curated to teach physical processes
located in the directory ``landlab/notebooks/teaching``.

Run on Binder
`````````````

To launch an instance of
Binder and `explore the notebooks click here`_.

.. _explore the notebooks click here: https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb

To launch a Binder instance that goes straight to the `teaching notebooks click here`_.

.. _teaching notebooks click here: https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/teaching/welcome_teaching.ipynb

Run on the CSDMS JupyterHub
```````````````````````````

The Landlab notebooks can also be run on the CSDMS JupyterHub.
To sign up for a free account,
`visit the CSDMS wiki`_ and follow the instructions there.
Then, click to explore the `example notebooks`_,
or to go straight to the `teaching notebooks`_.

.. _visit the CSDMS wiki: https://csdms.colorado.edu/wiki/JupyterHub
.. _example notebooks: https://jupyter.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=tree%2Flandlab%2Fnotebooks%2Fwelcome.ipynb&branch=master
.. _teaching notebooks: https://jupyter.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=tree%2Flandlab%2Fnotebooks%2Fteaching%2Fwelcome_teaching.ipynb&branch=master


License
-------

*landlab* is licensed under the MIT License.

Citing Landlab
--------------

.. start-citing-landlab

If you use any portion of Landlab, please cite the following papers:

.. tab:: APA

  `Hobley, D. E. J. <https://www.earth-surf-dynam.net/5/21/2017/>`__, Adams,
  J. M., Nudurupati, S. S., Hutton, E. W. H., Gasparini, N. M., Istanbulluoglu,
  E. and Tucker, G. E., 2017, Creative computing with Landlab: an open-source
  toolkit for building, coupling, and exploring two-dimensional numerical models
  of Earth-surface dynamics, Earth Surface Dynamics, 5(1), p 21-46,
  10.5194/esurf-5-21-2017.

  `Barnhart, K. R. <https://doi.org/10.5194/esurf-8-379-2020>`__,
  Hutton, E. W. H., Tucker, G. E., Gasparini, N. M., Istanbulluoglu, E.,
  Hobley, D. E. J., Lyons, N. J., Mouchene, M., Nudurupati, S. S., Adams, J. M.,
  and Bandaragoda, C., 2020, Short communication: Landlab v2.0: A software package for
  Earth surface dynamics, Earth Surf. Dynam., 8(2), p 379-397,
  doi:10.5194/esurf-8-379-2020.

  Hutton, E., Barnhart, K., Hobley, D., Tucker, G., Nudurupati, S., Adams, J., Gasparini, N., Shobe, C., Strauch, R., Knuth, J., Mouchene, M., Lyons, N., Litwin, D., Glade, R., Giuseppecipolla95, Manaster, A., Abby, L., Thyng, K., & Rengers, F. (2020). landlab [Computer software]. https://doi.org/10.5281/zenodo.595872

.. tab:: BibTeX

  ::

    @article{hobley2017creative,
      title={
        Creative computing with Landlab: an open-source toolkit for building,
        coupling, and exploring two-dimensional numerical models of
        Earth-surface dynamics
      },
      author={
        Hobley, Daniel EJ and Adams, Jordan M and Nudurupati, Sai Siddhartha and
        Hutton, Eric WH and Gasparini, Nicole M and Istanbulluoglu, Erkan and
        Tucker, Gregory E
      },
      journal={Earth Surface Dynamics},
      volume={5},
      number={1},
      pages={21--46},
      year={2017},
      publisher={Copernicus GmbH},
      url={https://esurf.copernicus.org/articles/5/21/2017/},
      doi={10.5194/esurf-5-21-2017}
    }

    @article{barnhart2020landlab,
      title={Landlab v2. 0: a software package for Earth surface dynamics},
      author={
        Barnhart, Katherine R and Hutton, Eric WH and Tucker, Gregory E and
        Gasparini, Nicole M and Istanbulluoglu, Erkan and Hobley, Daniel EJ and
        Lyons, Nathan J and Mouchene, Margaux and Nudurupati, Sai Siddhartha and
        Adams, Jordan M and others
      },
      journal={Earth Surface Dynamics},
      volume={8},
      number={2},
      pages={379--397},
      year={2020},
      publisher={Copernicus GmbH}
      url = {https://esurf.copernicus.org/articles/8/379/2020/},
      doi = {10.5194/esurf-8-379-2020}
    }

    @software{Hutton_landlab_2020,
    author = {Hutton, Eric and Barnhart, Katy and Hobley, Dan and Tucker, Greg and Nudurupati, Sai and Adams, Jordan and Gasparini, Nicole and Shobe, Charlie and Strauch, Ronda and Knuth, Jenny and Mouchene, Margaux and Lyons, Nathan and Litwin, David and Glade, Rachel and {Giuseppecipolla95} and Manaster, Amanda and Abby, Langston and Thyng, Kristen and Rengers, Francis},
    doi = {10.5281/zenodo.595872},
    license = {MIT},
    month = {4},
    title = {{landlab}},
    url = {https://github.com/landlab/landlab},
    year = {2020}
    }

.. end-citing-landlab

Citing Landlab Components
-------------------------

.. start-citing-components

If you are working with Landlab components and utilities, many of them have
their own publication. Please cite it to acknowledge the component authors.

Citation information for each component can be found as follows:

- Where relevant, software citation and general references, are listed in the
  Component API documentation under the References section.
- Software citations are included in component metadata. We have created a
  tool called the "Citation Registry" which creates a .bib file for software
  citations used in an application. See example usage :ref:`here <cite_as>`.

.. end-citing-components

Contact
-------

.. start-contact

The recommended way to contact the Landlab team is with a
`GitHub Issue <https://github.com/landlab/landlab/issues>`_.

* **Bug reports**: Please make an Issue describing the bug so we can address it, or work
  with you to address it. Please try to provide a `minimal, reproducible example
  <https://stackoverflow.com/help/minimal-reproducible-example>`_.
* **Documentation**: If something in our documentation is not clear to you, please make an
  issue describing the what isn't clear. Someone will tag
  the most appropriate member of the core Landlab team. We will work to clarify
  your question and revise the documentation so that it is clear for the next user.

Keep in touch with the latest *landlab* news by following us on `Twitter <https://twitter.com/landlabtoolkit>`_.

During workshops and clinics, we sometimes use the
`Landlab Slack channel <https://landlab.slack.com>`_.

.. end-contact



