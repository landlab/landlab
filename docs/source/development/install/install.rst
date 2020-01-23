.. _dev_install_install:

===========================
Install Landlab from Source
===========================

Before You Begin
----------------

1.  Ensure you have installed with Xcode from the
    Apple app store (macs) or :ref:`installed a working C++ compiler on your
    machine (PCs) <compile_in_windows>` before proceeding
    with your developer install. **You should also update your Python
    distribution!** For Anaconda, use

    .. code-block:: bash

        $ conda update --all

    (two dashes), and then separately,

    .. code-block:: bash

        $ conda update setuptools

    (the second being essential!) from your terminal.

2.  Ensure that you have removed other versions of Landlab from your computer.

    If you've previously used pip to install Landlab, we recommend
    you take that version off first. At a command prompt, use the command

    .. code-block:: bash

        $ pip uninstall landlab

    If you have used ``conda`` to install a prebuilt version of Landlab, you
    should uninstall that too.

    .. code-block:: bash

        $ conda uninstall landlab

    If you're not sure whether you have or not in the past, there's no harm
    doing both of these uninstall commands.

Install
-------

Now that you have a working copy of the Landlab code on you computer,
you need to install it. To install Landlab in developer mode, navigate
to the root Landlab folder (it will be landlab with a small ``l`` and
will contain the file ``setup.py``) and run the commands below.

Landlab has a number of dependencies to run, test, and develop with. These are
described in more detail :ref:`here <dependencies>`. We have created a conda
environment file which contains everything you will need for development.

Read more about
`conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_.
on the conda webpage.

The next step is it create this environment, called *landlab_dev*.

.. code-block:: bash

   $ conda env create --file=environment-dev.yml

The conda environment described by ``environment-dev.yml`` contains the minimal
set of dependencies necessary to run the Landlab tests and notebooks, and keep
the codebase clean and tidy. It may not include some of your favorite
development tools (e.g., spyder). See below for how to install additional
packages into the conda environment.

In addition, this environment does not have everything needed to build the
documentation. These requirements are specified in the file
``landlab/docs/environment.yml``.

Activate that environment so that you will be using that version of python and
all of the dependencies you just installed.

.. code-block:: bash

   $ conda activate landlab_dev

Install Landlab on your computer in such a way that Python always
imports Landlab from the working copy you just cloned. This ensures that
any changes you make to your copy of the code is seen by Python the
*next* time you import Landlab.

.. code-block:: bash

   $ python setup.py develop

Conda Environment Tips
----------------------

*   In order to use the ``landlab_dev`` environment created during installation,
    you will need to activate it every time you open a new terminal instance.
    Use the command: ``conda activate landlab_dev``.
*   If you use python tools for your work that are not in the conda environment
    but you have previously installed them on your computer (e.g., spyder),
    you will need to add them to the environment. Use the standard terminal
    calls to accomplish this.

.. code-block:: bash

   $ conda install name_of_package
   $ conda install name_of_package -c name_of_channel

Uninstall
---------

To uninstall your development version of Landlab (again from the root
``landlab/`` folder) run the following command:

.. code-block:: bash

   $ python setup.py develop -u

With Landlab uninstalled, you will no longer be able to import Landlab
from outside the root folder of your working copy.
