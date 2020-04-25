.. _anaconda_install:

================
Anaconda Install
================

Install Anaconda
----------------

This should work for `Anaconda <https://www.anaconda.com/distribution/>`_ users
with Windows 7+, Mac OS 10.6+, or
Ubuntu Linux (only the latest version has been tested).

Once you have a
:ref:`full Python distribution on your machine<python_installation>`,
it is vital to
check that it has been successfully set as the default copy of Python on
your system. Open a command prompt (Terminal on a Mac, or Command Prompt
on a PC) and type the lines below (note the `$` indicates that you are on
a command line):

MAC:

.. code-block:: bash

 $ which python
 $ which ipython

WINDOWS:

.. code-block:: bash

    $ where python
    $ where ipython

In each case, both commands should return the same path, and it should
clearly refer to Anaconda. Details will depend on your
operating system but it could look something like this:

.. code-block:: bash

    /anaconda/bin/python


If you don't see reference to your newly installed distribution (i.e.,
`/anaconda`), :ref:`click here to resolve the problem<correcting_install_paths>`.

Make sure you have the latest version installed (close anaconda before
doing this):

.. code-block:: bash

    $ conda update –all

Once the path to both `python` and `ipython` point to your new distribution,
open the Python editor in Anaconda called Spyder.

On the Spyder toolbar, go to `Tools → Open` command prompt to open the
command line.

Alternatively you can open a standard terminal window, such as an xterm
(X11.app) or terminal window (Terminal.app) on a Mac, or a command
prompt on a Windows machine. If you do use a standard terminal and run
into problems, make sure you have
:ref:`resolved your path issues<correcting_install_paths>`.

.. _anaconda_install_landlab:

Now to install Landlab!
-----------------------

You can either install with the conda or the pip package managers. Conda
is recommended, as it reduces the chances of versioning conflicts. Try
to remember which you choose to avoid confusion when updating later! (If
you installed landlab prior to May 19th 2016, you will have used pip).

Type either (for conda install):

.. code-block:: bash

    $ conda install landlab -c landlab -c conda-forge

... or for pip install:

.. code-block:: bash

    $ pip install landlab

.. _test_landlab_install:

Test Landlab install
--------------------

Once Landlab has been successfully installed, on the Python shell line,
check to make sure it is up-to-date (note that those are double
underscores around version; also note that you may need to close and
reopen Anaconda before typing the below commands):

.. code-block:: bash

    $ import landlab
    $ landlab.__version__


The version number should be greater than 1. You can check the version
number of the most recent release `here <https://github.com/landlab/landlab/releases>`_.

Install/Test problems
---------------------

If you are having problems when installing, testing or running Landlab,
please visit our :ref:`Troubleshooting<troubleshooting>` page.
