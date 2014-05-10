landlab
=======

The Landlab project creates an environment in which scientists can build a
numerical landscape model without having to code all of the individual
components. Landscape models compute flows of mass, such as water, sediment,
glacial ice, volcanic material, or landslide debris, across a gridded terrain
surface. Landscape models have a number of commonalities, such as operating on
a grid of points and routing material across the grid. Scientists who want to
use a landscape model often build their own unique model from the ground up,
re-coding the basic building blocks of their landscape model rather than
taking advantage of codes that have already been written.

Install landlab
===============

Installing with ``pip``
-----------------------

The most recent stable version of Landlab is available on `PyPI
<https://pypi.python.org/pypi>`_ and can be installed with the following
command::

    $ pip install TheLandlab

Use pip to upgrade an existing landlab installation,

    $ pip install TheLandlab --upgrade

This will install Landlab as well as any prerequisite packages (required packages
are listed in setup.py).

Installing from source
----------------------

From Git
>>>>>>>>

This is the recommended way to install from source, as it will make it easiest
to keep up with the latest bug fixes.

.. note::

    If you are planning on developing for landlab (making changes to the code)
    please see our developers' guide for installation instructions.

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

1. Download the latest tarball from `github <https://github.com/landlab/landlab/archive/master.zip>`_.

2. From the root directory where your unpacked Landlab, run::

    python setup.py install


:package name: TheLandlab
:version: 0.1
:release date: 2013-03-24
:authors:
  Greg Tucker,
  Nicole Gasparini,
  Erkan Istanbulluoglu,
  Daniel Hobley,
  Sai Nudurupati,
  Jordan Adams,
  Eric Hutton

:url: http://csdms.colorado.edu/trac/landlab

:license: MIT

