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

    If you’ve previously used pip to install Landlab, we recommend
    you take that version off first. At a command prompt, use the command

    .. code-block:: bash

        $ pip uninstall landlab

    If you have used ``conda`` to install a prebuilt version of Landlab, you
    should uninstall that too.

    .. code-block:: bash

        $ conda uninstall landlab

    If you’re not sure whether you have or not in the past, there’s no harm
    doing both of these uninstall commands.

Install
-------

Now that you have a working copy of the Landlab code on you computer,
you need to install it. To install Landlab in developer mode, navigate
to the root Landlab folder (it will be landlab with a small ``l`` and
will contain the file ``setup.py``) and run the following commands:

.. code-block:: bash

   $ conda env create --file=environment-dev.yml
   $ conda activate landlab_dev
   $ pip install -e .

This first command installs all of the dependencies required by Landlab
into a new environment called *landlab_dev*. Read more about
`conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_.
on the conda webpage. The second command
activates that environment so that you will be using that version of
python and all of the dependencies you just installed. The third command
installs Landlab on your computer in such a way that Python always
imports Landlab from the working copy you just cloned. This ensures that
any changes you make to your copy of the code is seen by Python the
*next* time you import Landlab.

Conda Environment Tips
----------------------

*   In order to use this environment, you will need to activate it every time
    you open a new terminal instance.
*   If you use python tools for your work that are not in the conda environment
    but you have previously installed them on your computer (e.g., spyder),
    you will need to add them to the environment. Use the standard
    `conda install name_of_package` or
    `conda install name_of_package -c name_of_channel`
    terminal calls to acomplish this.

Uninstall
---------

To uninstall your development version of Landlab (again from the root
``landlab/`` folder) run the following command:

.. code-block:: bash

   $ pip unintall landlab

With Landlab uninstalled, you will no longer be able to import Landlab
from outside the root folder of your working copy.
