.. _developer_install:

Developer Installation
======================

If you intend to modify the Landlab code base, including design of new
components, you will probably want to install Landlab directly from the
code. This way, when you modify the code, the changes you have made will
be reflected immediately in your Landlab library. This section provides
information on how to do this.

The process breaks down to three parts:

1. [[Fork|Installing-Landlab-from-source-code-(“developer-install”)#1-getting-the-code-from-github]]
   our code repository on `Github <http://github.com>`__ to give you a
   copy of the code on your own machine.
2. [[Install|Installing-Landlab-from-source-code-(“developer-install”)#2-installing-landlab-in-developer-mode]]
   the Landlab library from the code; updating the code will
   automatically update your install.
3. Periodically,
   [[update|Installing-Landlab-from-source-code-(“developer-install”)#3-updating-your-landlab-install]]
   your local Landlab code from Github to keep up to date with bug
   fixes, updates etc.

**Note:** *For dev work, we actively recommend Anaconda over the
Enthought Python Distribution, especially on Windows machines. This is
because it ships with a working compiler already associated with Python,
whereas the EPD does not. On a Mac, this is less important as the Xcode
app (available through iTunes) gives you the necessary compilers
instead—install it now if you don’t have it! If you choose to use the
EPD on a Windows machine, however,*\ `you’ll need to install separately
either Visual Basic or MinGW and successfully associate them with your
Python
install <http://landlab.readthedocs.org/en/latest/compilers_in_windows.html#compile-in-windows>`__\ *.
Email the development team if you’re really struggling. But unless
you’re really invested in Canopy and the EPD, uninstalling it and
replacing with Anaconda is probably the more stress-free way to go.*

*Either way, you’ll need a working C++ compiler running alongside Python
to be able to perform a full developer install. You’ll see errors
referring to [[ Cython \| Python,-NumPy,-SciPy,-Cython#cython]] if you
don’t have working compiler when calling ``python setup.py develop``
(see [[Section
2|Installing-Landlab-from-source-code-(“developer-install”)#2-installing-landlab-in-developer-mode]]).*

1. Getting the code from Github
===============================

Landlab development takes place in your own *fork* of the main Landlab
repository. A fork is a *mirror* of the repository and is hosted on your
personal GitHub account. You will use this fork for developing new
Landlab features. Your changes will migrate to the core repository (for
review and merging) by requesting that the main repository “pull” in
your changes. This is known as a pull request and is facilitated through
the GitHub website.

How to create a fork
--------------------

You will only need to do this once for each project to which you want to
contribute. Github has some great documentation on `how to create a
fork <https://help.github.com/articles/fork-a-repo/>`__. We outline
below the basic steps as applied to Landlab.

Create a Github account
~~~~~~~~~~~~~~~~~~~~~~~

1. You can create a GitHub account by going to the `GitHub
   website <http://github.com>`__.
2. If you haven’t already, `Install
   Git <https://help.github.com/articles/set-up-git/>`__. Note that if
   you are using a mac with OS lower than 10.7, we have found it
   difficult to setup git, and you may want to upgrade your OS.
   Alternatively, you can do this in a more user-friendly fashion by
   installing the `Github graphical user
   interface <https://desktop.github.com>`__ for your OS, which comes
   packaged with git for command line by default, and also lets you
   perform many git functions (updating your code, pushing changes back
   to github…) without having to deal with the command line interface.
3. Configure your account to allow write access. If you choose to
   install the git GUI, then it will set-up an SSH key for you. If you
   install on the command line, you might need some help with this, see
   `Generating SSH Keys on
   GitHub <https://help.github.com/articles/generating-an-ssh-key/>`__.
4. If you want to ensure the ssh key is set correctly, go to your home
   page on github and hit the account settings (wrench and screwdriver
   button in upper right corner). On this page hit the SSH keys tab on
   the left. This should show that you have a key for whatever computer
   you are currently working on. Note that you may have more than one
   key if you have installed git on more than one computer with the same
   user account.

Creating your own fork of Landlab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following steps will create a fork of the Landlab repository under
your github account.

1. Sign in to your GitHub account.
2. Go to the `Landlab home page <https://github.com/landlab/landlab>`__
   on GitHub.
3. Click the fork button in the upper-right corner of the page.

Once completed, you will be redirected to the home page for your own
copy of Landlab (github.com/UserName/landlab).

Cloning your fork to your computer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can clone the fork (that lives on the Github website) locally to
your computer either using the Github app (the GUI), or directly from
the command line using git. If you’ve never used git before, the app is
probably the way to go.

**A. Using the GUI/app:**

1. Sign in to github on the GUI.
2. Hit on your account on the left side of the GUI.
3. This will show you your fork of Landlab. Hit the clone option next to
   the fork. This will download the Landlab package to your local
   computer. If you are on a windows machine, this will put Landlab in
   the Documents/GitHub folder. If you are on a mac, you are given the
   option of where to put the downloaded Landlab package.

Alternatively, you can also clone your fork in the GUI directly through
the website. Navigate to the page for your fork on the web
(UserName/landlab) and hit the “Clone in Desktop” button
(monitor-with-downward-arrow) on the top/right hand side of the page.
Make sure you have the GUI installed and set up on your machine before
you try this for the most pain-free results. Also ensure you clone your
fork, not the main development version of Landlab (i.e., the version at
github.com/your-user-name/landlab, not the version at
github.com/landlab/landlab).

**B. Using the command line:**

Use the following commands from a terminal of your choice.

::

   > git clone git@github.com:your-user-name/landlab.git
   > cd landlab
   > git remote add upstream git://github.com/landlab/landlab.git

2. Installing Landlab in developer mode
=======================================

**Before you start:** Ensure you have installed with Xcode from the
Apple app store (macs) or [[installed a working C++ compiler on your
machine (PCs) \| Installing-Compilers-on-Windows ]] before proceeding
with your developer install. **You should also update your Python
distribution!** For Anaconda, use

``conda update --all``

(two dashes), and then separately,

``conda update setuptools``

(the second being essential!) from your terminal.

**Note:** \_This assumes **you have never put Landlab on your machine
before**. If you’ve previously used pip to install Landlab, we recommend
you take that version off first. At a command prompt, use the command

::

   > pip uninstall landlab

If you have used ``conda`` to install a prebuilt version of Landlab, you
should uninstall that too.

::

   > conda uninstall landlab

If you’re not sure whether you have or not in the past, there’s no harm
doing both of these uninstall commands.

Now that you have a working copy of the Landlab code on you computer,
you need to install it. To install Landlab in developer mode, navigate
to the root Landlab folder (it will be landlab with a small ``l`` and
will contain the file ``setup.py``) and run the following commands:

::

   > conda env create --file=environment-dev.yml
   > conda activate landlab_dev
   > pip install -e .

This first command installs all of the dependencies required by Landlab
into a new environment called *landlab_dev*. The second command
activates that environment so that you will be using that version of
python and all of the dependencies you just installed. The third command
installs Landlab on your computer in such a way that Python always
imports Landlab from the working copy you just cloned. This ensures that
any changes you make to your copy of the code is seen by Python the
*next* time you import Landlab.

To uninstall your development version of Landlab (again from the root
``landlab/`` folder) run the following command:

::

   > pip unintall landlab

With Landlab uninstalled, you will no longer be able to import Landlab
from outside the root folder of your working copy.

Testing your install
--------------------

In order to test your installation you’ll need to install the
```pytest`` <https://docs.pytest.org/en/latest/>`__ package that is used
to run the tests.

::

   > conda install pytest

Once ``pytest`` has been installed navigate to the main Landlab
directory (the one with ``setup.py`` in it) and type into a terminal:

::

   > pytest

This command will collect and run all of the tests. If you want to only
want to test one part of Landlab (perhaps a component you are working
on), you would run:

::

   > pytest path\to\directory\you\want\to\test

You may also want to see the code coverage of different parts of
Landlab. To do this, you’ll first need to install
```pytest-cov`` <https://pytest-cov.readthedocs.io/en/latest/readme.html>`__.

::

   > conda install pytest-cov

Then execute

::

   > pytest landlab --doctest-modules --cov=landlab --cov-report term-missing

from the main Landlab directory.

This will run the tests and print the coverage statistics (including the
missing line numbers) to the terminal.

As above, you can also run the coverage tools for a more specific
directory. For example, to run them for your current directory you could
execute

::

   > pytest . --doctest-modules --cov=. --cov-report term-missing

or to run them for a specific directory (for example, the
erosion_deposition submodule) stored as an environment variable you
would do the following:

::

   > TEST_DIR=landlab/components/erosion_deposition/
   > pytest $TEST_DIR --doctest-modules --cov=$TEST_DIR --cov-report term-missing

3. Updating your Landlab install
================================

It is very important to regularly update your code to keep up with bug
fixes, new features and improvements!

See `the Update
page <https://github.com/landlab/landlab/wiki/Updating-Landlab>`__ for
instructions.

Working with your local version of Landlab
==========================================

Obviously, feel free to just dive into modifying the code, but your life
in the future will be a bit easier if you follow some basic
recommendations for good work flow with git forks and branches. Even if
you have a working knowledge of using git in a collaborative project, we
highly recommend that you review [[this section|Developing with github
and git]] of the documentation to get a sense of how to track
modifications to your version of Landlab in a way that makes it easy to
(a) get updates to Landlab made by the development team and other
contributors and (b) contribute improvements and new features you
develop back to the community. For information about our in-house code
formatting conventions and standards, see [[here|Style-conventions]].

Troubleshooting
===============

I updated my working version and now it is broken. What do I do?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One possibility is that the landlab requirements changed between when
you originally installed landlab and when you updated landlab. To
address this, re-run the following lines and then test the installation.

::

   > conda install --yes --file=requirements.txt
   > python setup.py develop

What do I do if my pull request cannot be automatically merged?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the latest upstream/master and go to the master branch. Remember,
*do not develop here*. Always develop in a feature branch. Merge the
lastest upstream master with your master:

::

   > git fetch upstream
   > git checkout master
   > git merge upstream/master

Go to the branch on which you are developing and merge the lastest
upstream master with your branch:

::

   > git checkout <branch_name>
   > git merge upstream/master

Fix the conflicts. Do this by hand or with a merge editor. This is where
you decide how to integrate the conflicting changes. Since only you know
what and why you made the changes you did, this can only be done by you:

::

   > git merge tool

After everything has been fixed, commit the changes and push the changes
to the repository. The pull request will automatically be updated:

::

   > git commit
   > git push

Most of these steps have equivalents in the Github app. Use the
“changes” pane to identify where conflicts exist in your version, then
resolve them one by one. When you’re done, commit then sync the
un-conflicted version’s changes as if they were any other.

I’m seeing errors about Cython when I try to run my code/import Landlab. It used to be fine.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Very occasionally, local code updates or rebasing can break the compiled
code that lives in your local developer’s install. *Provided you used to
have a fully working Landlab install*, you can fix this by just calling
again from the main Landlab local folder

::

   > python setup.py develop

as described above in the main text. If this is happening when you call
this install function rather than when you try to actually run some
code/import Landlab, see immediately below.

I see errors about Cython when I try to *install* Landlab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see errors referring to Cython when you try to run
``python setup.py develop``, it indicates you have a problem with your
local compilers. This can happen both the first time you ever try this,
or also subsequently, apparently at random. On a Mac, check first that
you have the free Apple Xcode app (get it from the app store). If you do
have it already, typically this means Xcode has updated itself (this can
happen automatically without your knowledge!) and needs you to
re-authorize its permissions. Open the Xcode app manually, follow the
instructions it will give you, then try the install for Landlab again.
On a PC? Try updating Anaconda.

I’m still confused
~~~~~~~~~~~~~~~~~~

If you are having problems when installing, testing or running Landlab,
please visit our `Troubleshooting
page <https://github.com/landlab/landlab/wiki/Troubleshooting>`__.

The Landlab development team will be happy to hear from you. We
recommend that you either post a question to the [[ Landlab User Group
\| User-Guide#landlab-user-group]], or [[create an new issue
request|https://github.com/landlab/landlab/issues/new]], and we’ll try
to resolve your problem. When reporting your problem (in either place)
we recommend that you provide a minimal, complete, and verifiable
example which will help the development team and involved users
reproduce your problem and determine a solution. [[This page from Stack
Overflow|https://stackoverflow.com/help/mcve]] provides some background
on how to make a minimal, complete, and verifiable example.
