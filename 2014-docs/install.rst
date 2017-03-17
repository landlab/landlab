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
- netCDF4 (will run without, but recommended)
- Cython (only required if building Landlab from source code)

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
preassembled scientific python distribution. In all ases, we recommend `Continuum IO's 
Anaconda <https://store.continuum.io/cshop/anaconda/>`_, as its conda package manager 
makes controlling your Python packages easier. (It is also possible to use `Enthought's 
Canopy <https://www.enthought.com/products/canopy/>`_, but be aware you will need to sign 
up for an academic license with Enthought to take full advantage of its features.) 
Download and follow the appropriate instructions 
for your operating system/distribution. These collections already include compatible
(and in some cases accelerated) versions of all of landlab's dependencies. When the
distribution asks if you want to set it as the default Python for your system, say yes.  
Note that both Canopy and Anaconda also provide a front-end, or GUI environment, from 
which you can work, making coding, running code, and debugging relatively easy.

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
end), and it should clearly refer to Anaconda (or Canopy). Details will depend on your
operating system. For instance, Dan's Macbook Pro gives::

    /anaconda/bin/python

If you *don't* see reference to your newly installed distribution, click :ref:`here 
<correcting_python_version>` to resolve the problem.

.. _landlab-install:

Installing Landlab
==================

.. note::

    If you already have a Landlab install on your machine, see 
    :ref:`updating landlab, <LL_update>` below.

.. note::

    If you already had a Python distribution on your machine, but it's a bit old, 
    remember to update both the distribution itself and its internal packages
    before attempting a Landlab install, to make sure the necessary dependencies
    are up to date. Do this from the command prompt
    for Anaconda, using: *conda update --all* (two dashes), then also
    *conda update setuptools* (or through the GUI in Canopy). 
    Also update or install netCDF4 through conda if you need to.


Here we describe how to install the latest release package of Landlab.  Note that this method of installation hides the code behind Landlab.  If you are an experienced Landlab user and want to actually edit existing Landlab code and add to the Landlab repository, please follow the developers' installation instructions :ref:`here 
<dev_guide>`.

We here assume that you have read :ref:`the previous section <python-install>` and you have now installed a Python front-end  on your computer (which should have also installed a Python distribution) and that your default Python path is set correctly (more on Python path :ref:`here <correcting_python_version>`).

The instructions below describe the installation of Landlab with :ref:`Anaconda <landlab-install_with_anaconda>`on both PCs and Macs.  If you're running Linux, it's likely you've already got your system the way you like it, and you'll already know how to get Landlab running on your machine using only the :ref:`fast install directions <landlab-install_quickly>`.  If, however, you want a bit more advice on beating your Linux system into shape with regards to running Python and getting Landlab, you can follow :ref:`this link <dan_installs_on_linux>`.

.. _landlab-install_quickly:

Quick Landlab Install Instructions (For Experienced Python Users)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you are new to Python, you probably should see instructions :ref:`here for Anaconda users <landlab-install_with_anaconda>`.  Otherwise, if you don't like details, continue!

- Make sure your Python distribution is up to date, especially setuptools.

- Open a terminal (or the command prompt) and type the following::

   > pip install --upgrade pip
   > pip install landlab

.. _landlab-install_with_anaconda:

Installing Landlab Using using Anaconda  - Recommended Method
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This should work for Anaconda users with Windows 7+, Mac OS 10.6+, or Ubuntu Linux (only the latest version has been tested).

-	Open the Python editor in Anaconda called Spyder.

-	On the Spyder toolbar, go to **Tools → Open** command prompt to open the command line.  Alternatively you can open a standard terminal window, such as an xterm (X11.app) or terminal window (Terminal.app) on a Mac, or a command prompt on a Windows machine.  If you do use a standard terminal and run into problems, make sure you have :ref:`resolved your path issues <correcting_python_version>`).

- To ensure that your version of **pip** (a package installer) is up-to-date, enter the following command::

  > pip install --upgrade pip

- Next, make sure the necessary dependencies are up-to-date. The following conda command will update all Anaconda packages (Note the ``conda`` command below handles Anaconda-supported package installation and updates)::

  > conda update --all
  
- Installing also requires a fully up-to-date version of setuptools, which (irritatingly) is not updated by the *--all* call above. So also run::

  > conda update setuptools

- Once the Anaconda packages are updated and the correct version of pip is installed, now install **netCDF4**::

  > conda install netCDF4

- Now to install Landlab! Enter the following command::

  > pip install landlab

- Once Landlab has been successfully installed, on the **Python shell line**, check to make sure it is up-to-date (note that those are double underscores around version; also note that you may need to close and reopen Anaconda before typing the below commands)::

  >>> import landlab
  >>> landlab.__version__

The version number is changing rapidly at this point, but it should be something higher than 0.1.18.  If you are having problems with Landlab, check with the Landlab development team to make sure you have the latest version.

..
    .. _landlab-install_with_canopy:

    Installing Landlab using Enthought Canopy
    +++++++++++++++++++++++++++++++++++++++++

    This should work for Canopy users with Windows 7+ or Mac OS 10.6 and above.

    - Open the Python editor by clicking on the Canopy icon.

    -	On the “Welcome to Canopy” window, log in to your Enthought Account. This will give you access to the package manager and required subpackages. 

    - On the Canopy toolbar, go to **Tools → Package Manager** to install required dependencies.

    - In the Package Manager, search for and install the **pip** and **netCDF4** libraries.

    -	Once **pip** and **netCDF** are installed, go to the Canopy editor window. On the toolbar, go to **Tools → Canopy Terminal** to open the command line.  Alternatively you can open a standard terminal window, such as an xterm (X11.app) or terminal window (Terminal.app) on a Mac, or a command prompt on a Windows machine.  If you do use a standard terminal and run into problems, make sure you have :ref:`resolved your path issues <correcting_python_version>`).

    - Now to install Landlab! On the command line, enter the following command::

      > pip install landlab
  
    - Once Landlab has been successfully installed, on the Python shell line in the Canopy editor window, check to make sure it is up-to-date (note that those are double undersocres around version)

      >>> import landlab
      >>> landlab.__version__
  
    The version number is changing rapidly at this point, but it should be something higher than 0.1.12.  If you are having problems with Landlab, check with the Landlab development team to make sure you have the latest version.


Developer Installation - Installing from Landlab Source Code
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This is recommended only for users who have gotten a feel for Landlab and want to keep up with the absolute latest Landlab developments and contribute codes back to the Landlab repository.  If this is not you, please follow the standard installation instructions :ref:`above <landlab_install>`.  Otherwise, if you are ready to become a Landlab developer, follow :ref:`these directions <dev_install>`.


.. _LL_update:

Updating Landlab
================

As Landlab is still relatively early in its development cycle, the code will update 
fairly often and new release versions will become available.
To take advantage of new features and new library additions, we recommend you
**update Landlab** fairly frequently.

.. note::

    Whenever you update Landlab, we **strongly** recommend you also update your
    Python package! For Anaconda, use the conda package manager from a
    command prompt:

    > conda update --all #(two dashes)

    (From Canopy, use the GUI to update all the available new packages listed.)

If you installed Landlab through the instructions on this page, this is trivial.
Simply use pip again to update, like so::

    > pip install landlab --upgrade


However, if you have ever used another method to install Landlab on your machine,
this might not be adequate (i.e., pip will give you error messages).
The first thing to do in such a case is to try a full uninstall and reinstall::

    > pip uninstall landlab
    > pip install landlab


Still having problems? This probably means that some time early in our 
development cycle you installed Landlab with one of our old procedures. The clue
will be that you still have a (very out of date!) copy of the Landlab code
base somewhere on your machine. Another possibility is that you've previously
tried a :ref:`developer install <dev_install>`.
This procedure will also work in this case.

Try this:

In a terminal, navigate to the top level directory of
that old code, the one that contains the file *setup.py*.
This is likely to be *your_home_dir*/landlab, if you installed with git
and left all the defaults as is.
Then::

    > pip uninstall landlab #just to be on the safe side, may get errors again
    > python setup.py develop -u

This should remove the install, **if** you installed as a developer.

Still getting error messages? This means we're going to have to excise the
old Landlab install "by hand". You're looking to remove any reference to
Landlab that lives inside *your_python_install*/lib/python2.7/site-packages.
**Do this only after you've exhausted other possibilities, above**, as
packages like pip will get annoyed with you if you start manually deleting
their files if they installed them in the first place. To minimize the risk,
onc again make sure you have just run::

    > pip uninstall landlab

Then find your Python directory with::

    > which python

Find that folder, ignoring everything after and including the subfolder 
*bin*. Instead, go to *your_install*/lib/python2.7/site-packages. In here,
you should find one (or more) folders referrring to landlab, e.g.,
*landlab* or *landlab.egg-link*, or some other reference to 
*landlab.egg*. Delete these. Leave everything else as it is!

Now try another pip install::

    > pip install landlab

This should now take. *Still* having problems? This is probably multiple
versions of Python on your machine interfering with each other. Solve
that problem first, then return to trying to install Landlab.
See :ref:`here <correcting_python_version>` for some help. 
