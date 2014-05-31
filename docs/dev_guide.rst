.. _dev_guide:

=========================
landlab Developers' Guide
=========================

landlab development takes place in your own *fork* of the main landlab
repository. A fork is a *mirror* of the repository and is hosted on your
personal GitHub account. You will use this fork for developing new landlab
features. Your changes will migrate to the core repository (for review and
merging) by requesting that the main repository "pull" in your changes. This
is known as a pull request and is facilitated through the GitHub website.

How to create a fork
====================

You will only need to do this once for each project to which you want to
contribute. Github has some great documentation on
`how to create a fork <https://help.github.com/articles/fork-a-repo>`_. We
outline below the basic steps as applied to landlab.

Create a GitHub account
-----------------------

1. You can create a GitHub account by going to the
   `GitHub website <https://github.com>`_.

2. If you haven't already, `Install Git 
   <https://help.github.com/articles/set-up-git>`_.  Note that if you are using
   a mac with OS lower than 10.7, we have found it difficult to setup git, and
   you may want to upgrade your OS. 

3. Configure your account to allow write access. If you choose to install the
   git GUI, then it will set-up an SSH key for you.  If you install on the
   command line, you might need some help with this, see `Generating SSH Keys
   <https://help.github.com/articles/generating-ssh-keys>`_ on GitHub.

4. If you want to ensure the ssh key is set correctly, go to your home page on
   github and hit the account settings (wrench and screwdriver button in upper
   right corner).  On this page hit the SSH keys tab on the left.  This should
   show that you have a key for whatever computer you are currently working on.
   Note that you may have more than one key if you have installed git on more
   than one computer with the same user account.


Creating your own fork of landlab
---------------------------------

The following steps will create a fork of the landlab repository under your
github account.

  1. Sign in to your GitHub account.  
  2. Go to the `landlab home page <https://github.com/landlab/landlab>`_ on
     GitHub.
  3. Click the *fork* button in the upper-right corner of the page.

Once completed, you will be redirected to the home page for your own copy
of the landlab.


Cloning your fork to your computer
----------------------------------

This is done from the GUI by:

1. Sign in to git on the GUI.
2. Hit on your account on the left side of the GUI.
3. This will show you your fork of landlab.  Hit the clone option next to the fork.  This will download the landlab package to your local computer.  If you are on a windows machine, this will put landlab in the Documents/GitHub folder.  If you are on a mac, you are given the option of where to put the downloaded landlab package.

This is done from the command line with the following commands::

  git clone git@github.com:your-user-name/landlab.git
  cd landlab
  git remote add upstream git://github.com/landlab/landlab.git


.. _developer-install:

Installing landlab in developer mode
------------------------------------

Now that you have a working copy of landlab on you computer, you need to
install it. To install landlab in developer mode run the following command
from the root landlab folder (the one that contains `setup.py`)::

  python setup.py develop

This installs landlab on your computer in such a way that Python always
imports landlab from the working copy you just cloned. This ensures that any
changes you make to your copy of the code is seen by Python the *next* time
you import landlab.

To uninstall your development version of landlab (again from the root landlab
folder) run the following command::

  python setup.py develop -u

With landlab uninstalled, you will not longer be able to import landlab
from outside to root folder of your working copy.

To check you have correctly installed landlab, run the landlab tests.


Fetching updates to the trunk
-----------------------------

From time to time you should fetch commits to the trunk that you don't have
in your working copy. You do this with the following command::

  git fetch upstream


Making a new branch
-------------------

Before making any changes to your code, you should create a new branch.

Update your mirror with any upstream changes you don't have::

  git fetch upstream

Make the new branch::

  git branch name-of-branch upstream/master
  git checkout name-of-branch

You will probably want to choose a descriptive name for your new branch so that
you and others will remember what it is you are intending to do with your
branch (for example, `bugfix-for-that-major-problem`, or
`add-that-cool-feature`).

If you want to keep your branches on you public GitHub page for landlab (you
probably do) you need to tell git to push changes to your github repo. This
is done with the following command::

  git push --set-upstream origin name-of-branch

On your landlab GitHub page you will now be able to toggle between your
various branches to see the code you have committed.


Testing the landlab installation
================================

The easiest way to run the landlab tests is to do so from inside the Python
interpreter::

  >>> import landlab
  >>> landlab.test()

This will run a series of tests and print our the result of each test. If
there are any failures, you can report that at the `landlab issue tracker <https://github.com/landlab/landlab/issues>`_.


Coding Style
============

* Please stick to the coding style described by `PEP8
  <http://www.python.org/dev/peps/pep-0008/>`_.

* Class and function docstrings should follow the `numpydoc conventions
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

