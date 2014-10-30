.. _correcting_python_version:

==================================
Dealing with Python on your system
==================================

Known issues, and their symptoms
================================

Based on our experience of guiding users through the install process so far, we have 
found that the vast majority of technical issues arise due to conflicts between various
different copies of and versions of Python that tend to accumulate on laptops used 
frequently for scientific computing. Very many scientific computing applications are
built onto a Python framework (notably various flavours of GIS software), and in some
cases this can cause your machine to become confused over which version it should be using
at a given time.

Other causes of this problem may occur if you have ever messed about with your `.profile`,
`.bash_profile`, or `.bashrc` files (Mac/Linux) or your Environment Variables (PC). In
such situations, you're likely to need :ref:`The Hard Way <the_hard_way>`.


SYMPTOMS:
>>>>>>>>>

#. Typing ``which python`` and ``which ipython`` gives versions of Python that appear to 
   point to different locations.
#. The install procedure appeared to go fine, including the test import in the Python
   shell described at the end of the instructions, but now you can't import Landlab from
   the Python prompt provided inside the Canopy or Anaconda app.
#. As above, but you can't import Landlab from inside an iPython shell or from an iPython
   Notebook.
#. Landlab appears to import OK (though there may be warnings or errors at the time of
   import), but when you try to use it you see lots of errors and warnings indicating it
   can't find Python modules it needs (e.g., Numpy, SciPy).

You may have some subset of the above symptoms occurring together, but one is normally
enough to indicate Python conflict is the problem.


Resolving the problem:
>>>>>>>>>>>>>>>>>>>>>>

There are two ways to fix this problem; and easy way, and a hard way. Try the easy one 
first!

.. _the_easy_way:

The Easy Way
------------

Both Canopy and Anaconda perform tests to see if you have other versions of Python
already on your system when you install. This is why they ask you if you want to set
Canopy/Anaconda as the default Python package at install time. It's possible you clicked
`No` when you should have clicked `yes`, or for some other reason the software failed to
set the first time.

If you are trying to install Canopy, you're in luck; you can also do this after install. 
Go to `Preferences`. Under `General`, you will see text telling you Canopy is not the 
default environment, with a button to push saying `Set as default`. Press it. 

Unfortunately, in Anaconda, you have to set this up more yourself. See :ref:`The Hard Way 
<the_hard_way>`.

.. _the_hard_way:

The Hard Way
------------

If the above didn't work, you're going to have to doctor your system environment variables
yourself. What you do is system-dependent:

**MAC/LINUX**. Open a terminal. It should automatically open in your home directory, but 
just to be on the safe side, type::

    cd ~
    
Now get a list of all the files in this folder, including the hidden ones::

    ls -a

Have a thorough look at this list. You're looking for files called `.profile` and
`.bash_profile`. Which you have will depend on your system configuration. If you have 
*both*, this is probably the root cause of your problem; `.bash_profile` overrides
`.profile` if both are present. In this case, open `.profile` and copy everything you
find inside across to `.bash_profile`. If you've already tried The Easy Way, profile
will probably contain several (repeated?) entries clearly referring to Canopy or "EPD"
having tried to reset itself as default. You can safely delete any repeated text in this
file. (Remember, # is the "comment out" symbol in the Bash terminal). Make sure you
don't actually delete anything else, and we recommend you don't actually remove the 
.profile file. Use your favorite Unix editor to make the changes - almost all machines
have `nano`, which is super basic, but works. Once you're done, save (Ctrl-X, say yes to
save in `nano`), then quit and restart any apps that are using Python - including the
terminal app itself. Test if this solved your problem.

If you only have one of `.bash_profile` or `.profile`, open it in an editor and have a 
look. You're looking for any obvious references to Python on any of the uncommented-out
lines. You'll probably see several lines that look something like::

    PATH="Users/yourusername/Library/...:${PATH}"
    export PATH

These lines are modifying your environment property PATH, by appending the location of
whatever thing it is the line refers to to the start of PATH, which is a big long 
composite string containing lots of addresses. The key here is to understand that *things
nearer the front of the PATH list get precedence*. So if you see lots of `export PATH`
commands in this file that all refer to Python, your system is seeing *only the last
entry*. So find any references to the distribution you want, and duplicate them **at the
end of this file**. Alternatively, you can simply put the comment out symbol # in front
of the lines that are referring to the path you see when you type ``which python``.

In the rare case that you don't see any reference to Anaconda at all (with Canopy, The
Easy Way will have already embedded the necessary lines in this file), you will need to
add the lines yourself. Add the following to the end of the file::

    export PATH=~/anaconda/bin:$PATH

(Assuming you installed Anaconda in the default location).


**PC**. On a PC, the same principle of modifying your environment variables applies, but
you access them differently. Go to the Control Panel, then System. On Windows 8, you then
want Advanced System Settings, though this will be similar on older OSes. Go to 
Advanced, then to the `Environment Variables...` button. Under User Variables, see if 
there is an entry called PATH. If there is, we will modify it. If there isn't, we will
create one. It is **VERY IMPORTANT** that you do not modify any existing text, 
*especially* under `System Variables` below.

As is the situation for Mac, above, the system reads these PATH strings from left to 
right, and stops once it has found what it is looking for. It also reads User before
System variables. Hence, we want to add new 
strings to the left hand (start) of the existing text, if there is any.

First, scan the existing string(s) (including under System) to see if there is any 
reference to the Python distribution you are trying to set as default already there.
e.g., my User PATH (running Anaconda cleanly) currently reads::

    C:\Users\Dan\AppData\Local\Continuum\Anaconda;C:\Users\Dan\AppData\Local\Continuum\Anaconda\Scripts

If you find a reference or references like this to the version you're currently trying to
run, copy the text, and add it (repeated) at the start of the User string. Copy this
syntax - semicolons separate paths.

Note that modifying the User Variables will only affect the current user account. Add the
text - carefully!! - to the System Variables if you want the changes for all users.






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
    of the distributions described below, click `here <correcting_python_version>`.

On all platforms (Linux, Windows 7 or greater, and MacOS X), we recommend a
preassembled scientific python distribution, such as `Continuum IO's Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ or `Enthought's Canopy
<https://www.enthought.com/products/canopy/>`_ (we prefer to use Canopy but
any of these should be fine). Download and follow the appropriate instructions 
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
prompt (Terminal on a Mac, or Command Prompt on a PC) and type ``which python``, then
``which ipython``. In each case, path should be the same (except the (i)python at the 
end), and it should clearly refer to Canopy or Anaconda. Details will depend on your
operating system. For instance, Dan's Macbook Pro gives:

    /Users/danhobley/Library/Enthought/Canopy_64bit/User/bin/python

If you *don't* see reference to your newly installed distribution, click `here 
<correcting_python_version>` to resolve the problem.


Installing Landlab
==================

The recommended way to install landlab is from source, as it will make it
easiest to keep up with the latest bug fixes.

It is also possible to install a "classroom version" of Landlab, which will give you
(or more likely, your putative undergraduate class) a single, quick and easy, but
non-updateable snapshot of Landlab. For details on installing the classroom version,
click `here <classroom_install>`.

.. note::

    The following instructions assume you have a working version of `Git
    <http://git-scm.com/>`_ installed on your system. Git is a
    distributed version control system (DVCS) and source code management
    system. For an introduction to Git and DVCS, see the official
    `git documentation <http://git-scm.com/documentation>`_. Installing the
    Github graphical user interface (see below) will give you the necessary
     git tools.


.. _source-install:

Installing from source
----------------------

Our code lives in the `Github <https://github.com>`_ online code repository. For install, 
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
#. Now, leave the Github app and open a command prompt(PC) or Terminal(Mac/Unix). 
   Navigate to the root directory of your Landlab download (reminder: change directory
   in a prompt/terminal using the command ``cd``, then the name of the subfolder; 
   ``cd ..`` takes you up one folder level). This root directory will contain a file
   called `setup.py` (check with ``dir`` (PC) or ``ls`` (Mac/Linux)).
   From this directory, type at the prompt::

    python setup.py develop

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
        updating your installation of Landlab, click `here <updating_landlab>`.
        
    
    Finally, test everything worked. From the same command line, type::
    
      python
    
    An interactive Python window will open in the command line; the prompt will look like
    ``>>>``. From here, enter::
    
      import landlab
    
    If you are returned to the >>> prompt after a few moments, everything is fine. If you
    see an error message, you might have some problems with your install. See the `install
    FAQ page <install_FAQ>` for a list of known install issues, and their solutions. 
    
    Leave the Python shell by typing::
   
      exit()
      

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

    git clone https://github.com/landlab/landlab.git

#. Navigate From the root directory of your landlab clone (the folder that contains
   `setup.py`). From your likely current location this will probably just be 
   ``cd landlab``. From here, enter::

    python setup.py develop

#. Finally, test everything worked. From the same command line, type::
    
      python
    
   An interactive Python window will open in the command line; the prompt will look like
   ``>>>``. From here, enter::
    
      import landlab
    
   If you are returned to the >>> prompt after a few moments, everything is fine. If you
   see an error message, you might have some problems with your install. See the `install
   FAQ page <install_FAQ>` for a list of known install issues, and their solutions. 
   
   Leave the Python shell by typing::
   
      exit()


