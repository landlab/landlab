.. _install:

==================
Installing Landlab
==================

Dependencies
============

Landlab has the following dependencies:

- Python 2.7
- Numpy 1.8 or greater
- Scipy 0.12 or greater

If you don't already have these packages installed on your computer, simply
install one of the preassembled scientific Python collections described below
under :ref:`Installing Python <python-install>`.

.. note::

  Although not supported, Landlab can be used with Python 3.X by simply
  running `2to3 <http://docs.python.org/2/library/2to3.html>`_ on the source.

.. _python-install:

Installing Python
=================

.. note::

    If you have used a decent amount of scientific software on  your machine before, it is 
    likely that some of this software will have already installed some "pieces" of Python
    onto your system. Nevertheless, we *strongly* recommend that if you haven't before, 
    you download a whole Python distribution, as this will ensure that all of the Python 
    modules that Landlab needs to run are definitely present. Most of the bug reports we
    get about problems installing Landlab relate to conflicts with old versions of Python
    on machines. Common symptoms are running the python setup commands at the end of this
    file, but then not being able to load landlab.
    If you suspect this might be happening to you after you've installed one
    of the distributions described below, click :ref:`here <correcting_python_version>`.

On all platforms (Linux, Windows 7 or greater, and MacOS X), we recommend a
preassembled scientific python distribution, such as `Continuum IO's Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ or `Enthought's Canopy
<https://www.enthought.com/products/canopy/>`_ (we prefer to use Anaconda but
either of these should be fine). Download and follow the appropriate instructions 
for your operating system/distribution. These collections already include compatible
(and in some cases accelerated) versions of all of landlab's dependencies. When the
distribution asks if you want to set it as the default Python for your system, say yes.  Note that both Canopy and Anaconda also provide a front-end, or GUI environment, from which you can work, making coding, running code, and debugging relatively easy.

On Linux systems, you can also install Python and the Landlab dependencies
from your package manager. If you're running Linux but aren't that familiar
with handling Python packages in it, :ref:`this <dan_installs_on_linux>`
might help.

(Landlab uses `setuptools <https://pypi.python.org/pypi/setuptools>`_ for
packaging and is configured to automatically download and install the most
up-to-date version of its dependencies from `PyPI
<https://pypi.python.org/pypi>`_, if a satisfactory version is not already
installed.)

Once you have a full Python distribution on your machine, it is vital to check that
it has been successfully set as the default copy of Python on your system. Open a command
prompt (Terminal on a Mac, or Command Prompt on a PC) and type the lines below (note the ``>`` indicates that you are on a command line)::

  > which python
  > which ipython 

In each case, path should be the same (except the (i)python at the 
end), and it should clearly refer to Canopy or Anaconda. Details will depend on your
operating system. For instance, Dan's Macbook Pro gives::

    /Users/danhobley/Library/Enthought/Canopy_64bit/User/bin/python

If you *don't* see reference to your newly installed distribution, click :ref:`here 
<correcting_python_version>` to resolve the problem.

.. _landlab-install:

Installing Landlab
==================

Here we describe how to install the latest release package of Landlab.  Note that this method of installation hides the code behind Landlab.  If you are an experienced Landlab user and want to actually edit existing Landlab code and add to the Landlab repository, please follow the installation instructions :ref:`here 
<dev_guide>`.

We here assume that you have read :ref:`the previous section <python-install>` and you have now installed a Python front-end  on your computer (which should have also installed a Python distribution) and that your default Python path is set correctly (more on Python path :ref:`here <correcting_python_version>`).

Quick Install Instructions (For Experienced Python Users)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you are new to Python, you probably should see instructions :ref:`here for Anaconda users <landlab-install_with_anaconda>` and :ref:`here for Canopy users <landlab-install_with_canopy>`.  Otherwise, if you don't like details, continue!

- Open a terminal (or the command prompt) and type the following::

   > pip install --upgrade pip
   > pip install landlab

.. _landlab-install_with_anaconda:

Installing Using using Anaconda  - Recommended Method
++++++++++++++++++++++++++++++++++++++++++++++++++++++

This should work for Anaconda users with Windows 7, Mac OS 10.6+, or Ubuntu Linux (only the latest version has been tested).

-	Open the Python editor in Anaconda called Spyder.

-	On the Spyder toolbar, go to **Tools → Open** command prompt to open the command line.  Alternatively you can open a standard terminal window, such as an xterm (X11.app) or terminal window (Terminal.app) on a mac, which are both found in the Applications/Utilities directory.  If you do use a standard terminal, make sure you have :ref:`resolved your path issues <correcting_python_version>`).

- To ensure that your version of *pip* (a package installer) is up-to-date, enter the following command::

  > pip install -upgrade pip
  
- Once the correct version is installed, now install **netCDF4**.  (Note the ``conda`` command below handles Anaconda-supported package installation and updates)::

  > conda install netCDF4

- On the Python shell line in Anaconda, check the install of **netCDF4** to make sure it is up-to-date (note that those are double underscores around version)::

  >>> netCDF4.__version__

As of May 2015 this should return ``1.1.8``

- Now to install Landlab! Enter the following command ::

  > pip install landlab

- Once Landlab has been successfully installed, on the python sheel line, check to make sure it is up-to-date (note that those are double undersocres around version)

  >>> import landlab
  >>> landlab.__version__

The version number is changing rapidly at this point, but it should be something higher than 0.1.12.  If you are having problems with Landlab, check with the Landlab development team to make sure you have the latest version.

.. _landlab-install_with_canopy:

Installing using Enthought Canopy
+++++++++++++++++++++++++++++++++

This should work for Canopy users with Windows 7 or Mac OS 10.6 and above.

- Open the Python editor by clicking on the Canopy icon.

-	On the “Welcome to Canopy” window, log in to your Enthought Account. This will give you access to the package manager and required subpackages. 

- On the Canopy toolbar, go to Tools → Package Manager to install required dependencies.

- In the Package Manager, search for and install the **pip 6.1.1-1** and **netCDF4 1.1.7.1-2** libraries.

-	Once **pip** and **netCDF** are installed, go to the Canopy editor window. On the toolbar, go to Tools → Canopy Terminal to open the command line.  Alternatively you can open a standard terminal window, such as an xterm (X11.app) or terminal window (Terminal.app) on a mac, which are both found in the Applications/Utilities directory.  If you do use a standard terminal, make sure you have :ref:`resolved your path issues <correcting_python_version>`).

- Now to install Landlab! On the command line, enter the following command::

  > pip install landlab
  
- Once Landlab has been successfully installed, on the Python shell line in the Canopy editor window, check to make sure it is up-to-date (note that those are double undersocres around version)

  >>> import landlab
  >>> landlab.__version__
  
The version number is changing rapidly at this point, but it should be something higher than 0.1.12.  If you are having problems with Landlab, check with the Landlab development team to make sure you have the latest version.


Developer Installation - Installing from Source Code
++++++++++++++++++++++++++++++++++++++++++++++++++++

This is recommended only for users who have gotten a feel for Landlab and want to keep up with the absolute latest Landlab developments and contribute codes back to the Landlab repository.  If this is not you, please follow the standard installation instructions :ref:`above <landlab_install>`.  Otherwise, if you are ready to become a Landlab developer, follow :ref:`these directions <dev_guide>`.
