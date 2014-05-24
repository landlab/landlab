.. _install:

==================
Installing Landlab
==================

Dependencies
============

Landlab has the following dependencies:
- Python 2.7
- Numpy 1.7 or greater
- Scipy 0.12 or greater

If you don't already have these packages installed on your computer, simply install one 
of the preassembled scientific Python collections described below under *Installing
Python*.

.. note::

  Although not supported, Landlab can be used with Python 3.X by simply
  running `2to3 <http://docs.python.org/2/library/2to3.html>`_ on the source.


Installing Python
=================

On all platforms (Linux, Windows 7 or greater, and MacOS X), we recommend a
preassembled scientific python distribution, such as `Continuum IO's Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ or `Enthought's Canopy
<https://www.enthought.com/products/canopy/>`_ (we prefer to use Canopy but
any of these should be fine). Download and follow the appropriate instructions 
for your operating system/distribution. These collections already include compatible
(and in some cases accelerated) versions of all of landlab's dependencies.

On Linux systems, you can also install Python and the Landlab dependencies
from your package manager. If you're running Linux but aren't that familiar
with handling Python packages in it, :ref:`this <dan_installs_on_linux>`
might help.

(Landlab uses `setuptools <https://pypi.python.org/pypi/setuptools>`_ for
packaging and is configured to automatically download and install the most
up-to-date version of its dependencies from `PyPI
<https://pypi.python.org/pypi>`_, if a satisfactory version is not already
installed.)


Install Landlab
===============

There are several ways to install Landlab. Which method you use will 
depend on how you want to use Landlab: whether you want to view the source code, or 
simply run it as is.

- If you would like to see the Landlab source code and/or make changes or
  add to the code, you should follow :ref:`The Developer's Guide
  <dev_guide>`.
- If you want to use Landlab as-is, and not see or make changes to the
  code, you should follow the :ref:`Installing with pip <basic-install>`
  instructions.

.. _source-install:

Installing from source
----------------------

From Git
>>>>>>>>

This is the recommended way to install from source, as it will make it easiest
to keep up with the latest bug fixes.

.. note::

    If you are planning on developing for landlab (making changes to the code)
    please see our :ref:`Developers' Guide <dev_guide>` for installation instructions.

.. note::

    The following instructions assume you have a working version of `Git
    <http://git-scm.com/>`_ installed on your system. Git is a
    distributed version control system (DVCS) and source code management
    system. For an introduction to Git and DVCS, see the official
    `git documentation <http://git-scm.com/documentation>`_.


1. Clone landlab from the master repository, hosted on `github.com <http://www.github.com>`_::

    git clone https://github.com/landlab/landlab.git

2. From the root directory of your landlab clone (the folder that contains
   `setup.py`)::

    python setup.py install


From source tarball
>>>>>>>>>>>>>>>>>>>

1. Download the `latest tarball <https://github.com/landlab/landlab/archive/master.zip>`_
   from the `landlab github page <https://github.com/landlab/landlab/>`_.

2. From the root directory where your unpacked Landlab, run::

    python setup.py install


.. _basic-install:

Installing with pip
-------------------

.. note::

  If you are developing landlab you will probably not want to do this!
  If you do, you would have to run this command everytime you make a change to
  the code. Instead, you want to install landlab in "development mode". See
  the :ref:`The Developer's Guide <dev_guide>` section for details
  on how to do this.

Use this method if you would like to install the landlab onto your machine
so you can use it as-is. You will be able to import landlab from the Python
command line, and write scripts that use landlab modules but will not have
access to the source. If you would would like to see the code or make tweaks
to it (large or small) you should install landlab from source (see
:ref:`Installing from source <source-install>`, or the
:ref:`Developer's Guide <dev_guide>`).

The most recent stable release of Landlab is available at the `Python Package
Index <https://pypi.python.org/pypi>`_ and can be installed by running::

    pip install TheLandlab

This will install Landlab as well as any prerequisite packages. To upgrade
an existing Landlab installation::

    pip install TheLandlab --upgrade


