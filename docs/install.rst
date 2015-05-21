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
distribution asks if you want to set it as the default Python for your system, say yes.

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


Installing Landlab
==================

Classroom Version
+++++++++++++++++

If you are new to Landlab and just want to try it out, we recomend installing the latest release package. This is a quick and easy way to get a non-updateable snapshot of Landlab.  Installing the release package is ideal for users in a classroom or Landlab clinic.  The code behind Landlab will not be visable to the user with this installation method.

MORE STUFF HERE FROM JORDAN!

Developer Installation - Installing from Source Code
++++++++++++++++++++++++++++++++++++++++++++++++++++

If you have already gotten a feel for Landlab and want to commit for the longterm, we recommended installing landlab from the source, as it will make it easiest to keep up with the latest bug fixes and contribute codes back to the Landlab repository.

.. note::

    The following instructions assume you have a working version of `Git
    <http://git-scm.com/>`_ installed on your system. Git is a
    distributed version control system (DVCS) and source code management
    system. For an introduction to Git and DVCS, see the official
    `git documentation <http://git-scm.com/documentation>`_. Installing the
    Github graphical user interface (see below) will give you the necessary
    git tools.


.. _source-install:

The Landlab code lives in the `Github <https://github.com>`_ online code repository. For install, 
you have two choices. Firstly, you can sign up to Github as a user of that website, 
download their third party GUI, and use that to get a copy of the code. 
Alternatively, you can manage the code acquisition directly through the command line 
on your machine using the Git text interface. Using Github takes slightly longer, 
but their graphical interface is arguably more straightforward - especially for updating
Landlab once you have it installed.

.. _gui-install:

With GitHub GUI
>>>>>>>>>>>>>>>

#. Go to the `Github webpage <https://github.com>`_. The homepage will prompt you to sign
   up for an account. Do so! Choose the free plan. Note down your user name and password.
#. Install the `GitHub app 
   <https://help.github.com/articles/set-up-git>`_. Follow the directions for
   installing the native app for your operating system.
     * `Mac <https://mac.github.com>`_
     * `Windows <https://windows.github.com>`_
     * Linux: Follow the command-line :ref:`installation instructions
       <command-line-install>`.
#. Open the app. You need to provide it with your user name and password to allow it to
   interact smoothly with the website. You should be prompted to do so when it boots up
   for the first time. If not, go to Preferences and enter your sign-in details. Click 
   through the remainder of the options, skipping the "add repositories" step.
#. With your browser, go to the `landlab page
   <https://github.com/landlab/landlab>`_ on GitHub and click the "Clone in
   Desktop" button (midway down the right hand side of the page). This will automatically
   cause your machine to switch back to the Github app and begin the download process. 
   Pick a location to store the Landlab files on your hard drive, and click through.
   Download will begin.
#. Now, leave the Github app and open a command prompt (PC) or Terminal (Mac/Unix). 
   Navigate to the root directory of your Landlab download (reminder: change directory
   in a prompt/terminal using the command ``cd``, then the name of the subfolder; 
   ``cd ..`` takes you up one folder level). This root directory will contain a file
   called `setup.py` (check with ``dir`` (PC) or ``ls`` (Mac/Linux)).
   From this directory, type at the prompt::

  ``> python setup.py develop``

.. note::
    
    This command tells your install of Python that `landlab` is a Python module that 
    you have now installed on your system, and where to look for the files it needs
    to run. Using the keyword `develop` warns Python that the code you have saved 
    on your disc might change from time to time. This
    means that should you so desire, you can make changes to the code, add 
    functionality, add your own modules, or otherwise tinker with the .py files you
    will find in the directories that Github has placed on your system. Importantly,
    however, it also allows to you quickly and easily use Github to download more
    up-to-date versions of Landlab - which may contain bug fixes, etc. For more on
    updating your installation of Landlab, click :ref:`here <updating_landlab>`.
        
    
#. Finally, test everything worked. From the same command line, type::
    
>>> python
    
   An interactive Python window will open in the command line; the prompt will look like
   ``>>>``. From here, enter::
    
       >>> import landlab
    
   If you are returned to the >>> prompt after a few moments, everything is fine. If you
   see an error message, you might have some problems with your install. See the 
   :ref:`install FAQ page <install_FAQ>` for a list of known install issues, and their 
   solutions. 
   
   Leave the Python shell by typing::
   
       >>> exit()
      

.. _command-line-install:

With Git
>>>>>>>>

.. note::

    This assumes that you already have Git on your machine. To check, open a command 
    prompt and type ``git``. If you have it, you will see usage instructions. If you
    don't, you will see an error message.

#. Using the command prompt, clone landlab from the master repository. This is 
   hosted on `github.com <http://www.github.com>`_. The files will be added inside 
   whichever directory you are in when you enter this command.::

    > git clone https://github.com/landlab/landlab.git

#. Navigate From the root directory of your landlab clone (the folder that contains
   `setup.py`). From your likely current location this will probably just be 
   ``cd landlab``. From here, enter::

    > python setup.py develop

#. Finally, test everything worked. From the same command line, type::
    
      > python
    
   An interactive Python window will open in the command line; the prompt will look like
   ``>>>``. From here, enter::
    
      >>> import landlab
    
   If you are returned to the >>> prompt after a few moments, everything is fine. If you
   see an error message, you might have some problems with your install. See the 
   :ref:`install FAQ page <install_FAQ>` for a list of known install issues, and their 
   solutions. 
   
   Leave the Python shell by typing::
   
      >>> exit()

You can find more details about installing Landlab as a developer :ref:`here 
<dev_guide>`.
