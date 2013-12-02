============
Installation
============

Dependencies
============

The Landlab has the following dependencies:

- Python 2.7
- Numpy 1.7 or greater
- Scipy 0.12 or greater

.. note::

  Although not supported, The landlab can be used with Python 3.X by simply
  running 2to3 on the source.


Installing Python and The Landlab Dependencies
==============================================

On all platforms (Linux, Windows 7 or greater, and MacOS X), we recommend a
preassembled scientific python distribution, such as `Continuum IO's Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ or `Enthought's Canopy
<https://www.enthought.com/products/canopy/>`_ (we prefer to use Canopy but
any of these should be fine). These collections already include compatible
(and in some cases accelerated) versions of all of landlab's dependencies.
Download and follow the appropriate instructions for your operating
system/distribution.

On Linux systems, you can also install Python and the landlab dependencies
from your package manager.

landlab uses `setuptools <https://pypi.python.org/pypi/setuptools>`_ for
packaging and is configured to automatically download and install the most
up-to-date version of its dependencies from `PyPI
<https://pypi.python.org/pypi>`_, if a satisfactory version is not already
installed.


Installing The Landlab
======================

Install with ``pip``
--------------------

Use this method if you would like to install the landlab onto your machine
so you can use it as-is. You will be able to import landlab from the Python
command line, and write scripts that use landlab modules but will not have
access to the source. If you would would like to see the code or make tweaks
to it (large or small) you should get a copy of the source code (see
below on how to do this).

The most recent stable release of landlab is available at the `Python Package
Index <https://pypi.python.org/pypi>`_ and can be installed by running::

    pip install TheLandlab

This will install Landlab as well as any prerequisite packages. To upgrade
landlab using pip::

    pip install TheLandlab --upgrade


Installing from a source distribution
-------------------------------------

This is the recommended way to install the landlab from source.

.. note::

  This section assumes that you have a version of `Subversion (svn)
  <http://mercurial.selenic.com/>`_ installed for your operating system. 
  Subversion is a version control system (VCS).  For an introduction to
  Subversion, see `<http://svnbook.red-bean.com/>`_.

1. Checkout a version of landlab using Subversion::

    svn checkout https://csdms.colorado.edu/svn/TheLandlab/trunk landlab

2. From the root directory of the landlab package (this is the directory
   that contains the file, ``setup.py``::

    python setup.py install

This will put the landlab package into a system-wide location so that you can
import landlab from any python session. It will also install any required
packages.

.. note::

  If you are developing landlab you will probably not want to do this!
  If you do, you would have to run this command everytime you make a change to
  the code. Instead, you want to install the package in "development mode". See
  below in the "Development Environment" section for details on how to do this.


Install a Development Version
=============================

This method with get a copy of the landlab source code, and install it in
such a way that you can make changes to the code and the changes will
immediately be seen by Python.

From the base directory of your landlab repository (the folder that contains
setup.py), run the following::

  pip install -e .

This sets things up so that Python looks into your source directory for
landlab modules, rather than copying your source into the Python site-packages
directory. This ensures that any changes you make to sourece files will be
seen by Python.


Test Your Installation
======================

Once you have installed landlab through any of the above methods, you
optionally run some tests to see if your installation is working (or rather,
if it isn't working). From the Python command line, run::

  >>> import landlab
  >>> landlab.test()

If this results in any errors, please report them to the landlab team.

