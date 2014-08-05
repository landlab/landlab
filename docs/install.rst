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

If you don't already have these packages installed on your computer, simply
install one of the preassembled scientific Python collections described below
under *Installing Python*.

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


Installing Landlab
==================

The recommended way to install landlab is from source, as it will make it
easiest to keep up with the latest bug fixes.

.. note::

    The following instructions assume you have a working version of `Git
    <http://git-scm.com/>`_ installed on your system. Git is a
    distributed version control system (DVCS) and source code management
    system. For an introduction to Git and DVCS, see the official
    `git documentation <http://git-scm.com/documentation>`_.


.. _source-install:

Installing from source
----------------------

.. _gui-install:

With GitHub GUI
>>>>>>>>>>>>>>>

#. Install the `GitHub app 
   <https://help.github.com/articles/set-up-git>`_. Follow the directions for
   installing the native app for your operating system.
     * `Mac <https://mac.github.com>`_
     * `Windows <https://windows.github.com>`_
     * Linux: Follow the command-line :ref:`installation instructions
       <command-line-install>`.
#. With your browser, go to the `landlab page
   <https://github.com/landlab/landlab>`_ on GitHub and click the "Clone in
   Desktop" button.
#. From the root directory of your landlab clone (the folder that contains
   `setup.py`)::

    python setup.py install

.. _command-line-install:

With Git
>>>>>>>>

.. note::

    If you are planning on developing for landlab (making changes to the code)
    please see our :ref:`Developers' Guide <dev_guide>` for installation instructions.

#. Clone landlab from the master repository, hosted on `github.com <http://www.github.com>`_::

    git clone https://github.com/landlab/landlab.git

#. From the root directory of your landlab clone (the folder that contains
   `setup.py`)::

    python setup.py install


