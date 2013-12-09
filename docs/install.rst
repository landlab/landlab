.. _install:

==================
Installing Landlab
==================

Dependencies
============

The Landlab has the following dependencies:

- Python 2.7
- Numpy 1.7 or greater
- Scipy 0.12 or greater

.. note::

  Although not supported, landlab can be used with Python 3.X by simply
  running `2to3 <http://docs.python.org/2/library/2to3.html>`_ on the source.


Install Python
==============

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

Landlab uses `setuptools <https://pypi.python.org/pypi/setuptools>`_ for
packaging and is configured to automatically download and install the most
up-to-date version of its dependencies from `PyPI
<https://pypi.python.org/pypi>`_, if a satisfactory version is not already
installed.


Install Landlab
===============

There are several ways to install Landlab. Which method you use will 
depend on how you want to use landlab.

- If you want to use Landlab as-is, and not see or make changes to the
  code, you should follow the :ref:`Basic Installation <basic-install>`
  instructions.

- If you would like to see the Landlab source code and/or make changes or
  add to the code, you should follow :ref:`Install as a Developer
  <developer-install>` instructions.

.. _basic-install:

Basic Installation
------------------

Use this method if you would like to install the landlab onto your machine
so you can use it as-is. You will be able to import landlab from the Python
command line, and write scripts that use landlab modules but will not have
access to the source. If you would would like to see the code or make tweaks
to it (large or small) you should get a copy of the source code (see
below for a description of :ref:`how to do this <source-install>`).

The most recent stable release of Landlab is available at the `Python Package
Index <https://pypi.python.org/pypi>`_ and can be installed by running::

    pip install TheLandlab

This will install Landlab as well as any prerequisite packages. To upgrade
an existing Landlab installation::

    pip install TheLandlab --upgrade


.. _developer-install:

Install as a Developer
----------------------

Use this method if you want to see and/or mess with the Landlab source code.

.. note::

   This section assumes that you have a version of `Subversion (svn)
   <http://subversion.apache.org/>`_ installed for your operating system. 
   Subversion is a version control system (VCS).  For an introduction to
   Subversion, see the Subversion `quick start guide
   <http://subversion.apache.org/quick-start>`_


1. Checkout a version of Landlab using Subversion::

    svn checkout https://csdms.colorado.edu/svn/TheLandlab/trunk landlab

2. From the base directory of your landlab repository (the folder that
   contains ``setup.py``), run the following::

    pip install -e .

This sets things up so that Python looks into your source directory for
landlab modules, rather than copying your source into the Python site-packages
directory. This ensures that any changes you make to sourece files will be
immediately seen by Python.


.. _source-install:

Install the Latest Version
--------------------------

This is the recommended way to install the latest version of landlab from
source.

.. note::

   This section assumes that you have a version of `Subversion (svn)
   <http://subversion.apache.org/>`_ installed for your operating system. 
   Subversion is a version control system (VCS).  For an introduction to
   Subversion, see the Subversion `quick start guide
   <http://subversion.apache.org/quick-start>`_


1. Checkout a version of Landlab using Subversion::

    svn checkout https://csdms.colorado.edu/svn/TheLandlab/trunk landlab

2. From the root directory of the landlab package (this is the directory
   that contains the file, ``setup.py``)::

    python setup.py install

This will put the landlab package into a system-wide location so that you can
import landlab from any python session. It will also install any required
packages.

.. note::

  If you are developing landlab you will probably not want to do this!
  If you do, you would have to run this command everytime you make a change to
  the code. Instead, you want to install the package in "development mode". See
  the :ref:`Install as a Developer <developer-install>` section for details
  on how to do this.


Test Your Installation
----------------------

Once you have installed Landlab through one of the above methods, you can
optionally run some tests to see if your installation is working (or rather,
if it isn't working). From the Python command line, run::

  >>> import landlab
  >>> landlab.test()

If this results in any errors, please report them to the `landlab team <huttone@colorado.edu>`_.

