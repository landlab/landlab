.. _install:

=========================
Installation Instructions
=========================

In order to use *landlab* you will first need Python >=3.8. While not
necessary, we recommend using the 
`Anaconda Python distribution <https://www.anaconda.com/distribution/>`_
as it provides a large number of third-party packages useful for
scientific computing.

To install *landlab*, simply run the following in your terminal of choice:

.. tab:: conda

  .. code-block:: bash

    $ conda install landlab -c conda-forge

.. tab:: pip

  .. code-block:: bash

    $ pip install landlab

or, if you simply can't wait for the latest release, you can install *landlab*
directly from GitHub,

.. code-block:: bash

   $ pip install git+https://github.com/landlab/landlab


Source code
-----------

*landlab* is actively being developed on GitHub, where the code is freely-available.
If you would like to modify or contribute code, you can either clone our
repository

.. code-block:: bash

   $ git clone git://github.com/landlab/landlab.git

or download the `tarball <https://github.com/landlab/landlab/tarball/master>`_
(a zip file is available for Windows users):

.. code-block:: bash

   $ curl -OL https://github.com/landlab/landlab/tarball/master

Once you have a copy of the source code, you can install it into your current
Python,

.. tab:: conda

  .. code-block:: bash

     $ cd landlab
     $ conda install --file=requirements.txt
     $ pip install -e .

.. tab:: pip

  .. code-block:: bash

     $ cd landlab
     $ pip install -e .


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
