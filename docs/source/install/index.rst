.. _install:

=========================
Installation Instructions
=========================

.. include:: ../../../README.rst
  :start-after: .. start-install-release
  :end-before: .. end-install-release

From source
-----------

.. include:: ../../../README.rst
  :start-after: .. start-install-release
  :end-before: .. end-install-release

Updating and uninstalling
-------------------------

You can update an existing *landlab* installation to a newer version
(if available) by running the following,

.. tab:: conda

  .. code-block:: bash

    $ conda update landlab

.. tab:: pip

  .. code-block:: bash

    $ pip -U install landlab


To uninstall *landlab*,

.. tab:: conda

  .. code-block:: bash

    $ conda uninstall landlab

.. tab:: pip

  .. code-block:: bash

    $ pip uninstall landlab
    

.. _conda_environment:

Conda Environment
-----------------

We have specified a conda environment which will install Landlab and everything
else you need to run the Landlab Notebooks. First get a compressed copy of the
latest release `here <https://github.com/landlab/landlab/releases>`_, or you can
clone the github repo.

Next you will need to have conda installed on your machine. In a terminal
window/command prompt, navigate to the Landlab directory. Then create the
environment with the following command.

.. code-block:: bash

   $ conda env create --file=environment.yml

This will create a Landlab environment called ``landlab_notebooks``. Activate
that environment so that you will be using that version of python and all of
the dependencies you just installed.

.. code-block:: bash

   $ conda activate landlab_notebooks

You will need to activate this environment every time you want to use it. See
:ref:`tutorials <tutorials>` for more about Landlab Notebooks.

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
   correcting_paths
   troubleshooting
