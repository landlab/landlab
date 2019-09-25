.. _install:

=========================
Installation Instructions
=========================

If you want to use Landlab to explore Earth surface dynamics, this is almost
certainly the correct place to get install instruction. If you want to change
the Landlab source code go to the
:ref:`developer install instructions <developer_install>`.

In order to run Landlab you will need a python distribution. We recommend the
`Anaconda distribution <https://www.anaconda.com/distribution/>`_
(version 3.6 or higher).

You can install Landlab using either the pip or conda package management tools.
We distribute Landlab through `PyPI <https://pypi.org/project/landlab/>`_
and `conda-forge <https://anaconda.org/conda-forge/landlab>`_.

Conda instructions
------------------

Installation
````````````
In a terminal type:

.. code-block:: bash

  $ conda install landlab -c conda-forge

If you work with many different packages that require conflicting dependencies,
consider reading about (and using)
`conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_.


Updating
````````

In a terminal type:

.. code-block:: bash

  $ conda update landlab

Uninstall
`````````

In a terminal type:

.. code-block:: bash

    $ conda uninstall landlab

Pip Instructions
----------------

Installation
````````````
In a terminal type:

.. code-block:: bash

  $ pip install landlab

Updating
````````

In a terminal type:

.. code-block:: bash

  $ pip update landlab

Uninstall
`````````

In a terminal type:

.. code-block:: bash

    $ pip uninstall landlab


Additional Resources
--------------------

Over time the core Landlab development team has written resources for
different issues related to installation. These are below. Note that some are
very out of date. For example, one of them describes an installation of
python 2.7. We do not recommend installing a python 2.x version.

.. toctree::
   :maxdepth: 2

   installing_python
   rough_guide_to_python
   anaconda_install
   installing_grass
   installing_windows_compiler
   correcting_paths
   troubleshooting
