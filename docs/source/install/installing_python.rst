.. _python_installation:

===================
Python Installation
===================

If you have used a decent amount of scientific software on  your machine before, it is
likely that have already installed "pieces" of Python
onto your system. Nevertheless, we *strongly* recommend that if you haven't before,
you download a complete Python distribution, as this will ensure that all of the Python
modules that Landlab needs to run are definitely present. Most of the bug reports we
get about problems installing Landlab relate to conflicts with older or incomplete versions of Python.
Common symptoms are running the python setup commands at the end of this
file, but then not being able to load landlab. Try :ref:`the Test <test_landlab_install>`
to ensure your current distribution with all its dependencies is the default.
If you suspect a path or version problem after you've installed one
of the distributions described below, click :ref:`here <correcting_install_paths>`.

On all platforms (Linux, Windows 7 or greater, and MacOS X), we recommend a
preassembled scientific python distribution such as `Continuum IO's
Anaconda <https://www.anaconda.com/distribution/>`_. Its conda package manager
makes controlling your Python packages easier. (It is also possible to use `Enthought's
Canopy <https://assets.enthought.com/downloads/>`_, but be aware you will need to sign
up for an academic license with Enthought to take full advantage of its features.)

Download and follow the appropriate instructions
for your operating system/distribution. These collections already include compatible
(and in some cases accelerated) versions of all of Landlab's dependencies. **If the
distribution asks if you want to set it as the default Python for your system, say yes.**
Note that both Canopy and Anaconda also provide a front-end, or GUI environment, from
which you can work, making coding, running code, and debugging relatively easy.

**Note if using Anaconda:**
*There have been documented issues with resolution with default inline plotting within the Spyder IDE iPython console. To generate dynamic plots (e.g. Matlab-like plots), we recommend you change the graphics settings in Spyder after installation by following this work flow:*

In *Spyder -> Preferences -> iPython console -> Graphics -> Graphics Backend -> Automatic -> Apply -> OK -> Make sure to restart Spyder to update the preferences.*

On Linux systems, you can also install Python and the Landlab dependencies
from your package manager. If you're running Linux but aren't that familiar
with handling Python packages in it, :ref:`this <rough_guide>`
might help.

(Landlab uses `setuptools <https://pypi.org/project/setuptools/>`_ for
packaging and is configured to automatically download and install the most
up-to-date version of its dependencies from `PyPI
<https://pypi.org/>`_, if a satisfactory version is not already
installed.)

The Test
--------

Once you have a full Python distribution on your machine, it is vital to check that
it has been successfully set as the default copy of Python on your system.

On Linux or Mac, open a Terminal window and and type the lines below (note the ``>`` indicates that you are on a command line)::

  > which python
  > which ipython

This will show your default path to Python and iPython. In each case the path to the file you asked for should be the same and it should clearly refer to Anaconda (or Canopy). Details will depend on your
operating system. For instance, Dan's Macbook Pro shows ``python`` on the path::

   /anaconda/bin/python

On a PC, run ``python`` from your Command Prompt window. You should see a reference to anaconda when it is running.

If you *don't* see reference to your newly installed distribution or the file is not found, click :ref:`here <correcting_install_paths>` to resolve the problem. 

Once you have installed a complete Python distribution on your machine, follow :ref:`these instructions <install>` to install Landlab.
