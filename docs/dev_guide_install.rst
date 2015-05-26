.. _dev_guide_install:

==================================
Installing Landlab for Development
==================================

Landlab development takes place in your own *fork* of the main Landlab
repository. A fork is a *mirror* of the repository and is hosted on your
personal GitHub account. You will use this fork for developing new Landlab
features. Your changes will migrate to the core repository (for review and
merging) by requesting that the main repository "pull" in your changes. This
is known as a pull request and is facilitated through the GitHub website.


How to create a fork
====================

You will only need to do this once for each project to which you want to
contribute. Github has some great documentation on
`how to create a fork <https://help.github.com/articles/fork-a-repo>`_. We
outline below the basic steps as applied to Landlab.


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


Creating your own fork of Landlab
---------------------------------

The following steps will create a fork of the Landlab repository under your
github account.

1. Sign in to your GitHub account.  
2. Go to the `landlab home page <https://github.com/landlab/landlab>`_ on
   GitHub.
3. Click the *fork* button in the upper-right corner of the page.

Once completed, you will be redirected to the home page for your own copy
of Landlab.


Cloning your fork to your computer
----------------------------------

This is done from the GUI by:

1. Sign in to git on the GUI.
2. Hit on your account on the left side of the GUI.
3. This will show you your fork of Landlab.  Hit the clone option next to the
   fork.  This will download the Landlab package to your local computer.  If
   you are on a windows machine, this will put Landlab in the Documents/GitHub
   folder.  If you are on a mac, you are given the option of where to put the
   downloaded Landlab package.

This is done from the command line with the following commands::

  > git clone git@github.com:your-user-name/landlab.git
  > cd landlab
  > git remote add upstream git://github.com/landlab/landlab.git

You can also clone your fork in the GUI directly through the website; navigate
to the page for your fork on the web (UserName/landlab) and hit the "Clone in
Desktop" button on the right hand side of the page. Make sure you have the GUI
installed and set up on your machine before you try this for the most 
pain-free results.


.. _developer-install:

Installing Landlab in developer mode
------------------------------------

Now that you have a working copy of the Landlab code on you computer, you need to
install it. To install Landlab in developer mode run the following command
from the root Landlab folder (it will be landlab with a small *l* and will contain `setup.py`)::

  > python setup.py develop

This installs Landlab on your computer in such a way that Python always
imports Landlab from the working copy you just cloned. This ensures that any
changes you make to your copy of the code is seen by Python the *next* time
you import Landlab.

To uninstall your development version of Landlab (again from the root landlab
folder) run the following command::

  > python setup.py develop -u

With Landlab uninstalled, you will not longer be able to import Landlab
from outside to root folder of your working copy.

To check you have correctly installed Landlab, run the Landlab tests.
Do this by importing landlab in an interactive Python shell, then calling
*landlab.test()*.


Fetching updates to the trunk
-----------------------------

From time to time you should fetch commits to the trunk that you don't have
in your working copy. You do this with the following command::

  > git fetch upstream


Making a new branch
-------------------

Before making any changes to your code, you should create a new branch.

Update your mirror with any upstream changes you don't have::

  > git fetch upstream

Make the new branch::

  > git branch name-of-branch upstream/master
  > git checkout name-of-branch

You will probably want to choose a descriptive name for your new branch so that
you and others will remember what it is you are intending to do with your
branch (for example, `bugfix-for-that-major-problem`, or
`add-that-cool-feature`).

If you want to keep your branches on you public GitHub page for Landlab (you
probably do) you need to tell git to push changes to your github repo. This
is done with the following command::

  > git push --set-upstream origin name-of-branch

On your Landlab GitHub page you will now be able to toggle between your
various branches to see the code you have committed.

The GUI also offers fairly simple, intuitive, and powerful control over 
branching and merging, for those uninclined to use the command line
tools.


Testing the Landlab installation
================================

The easiest way to run the Landlab tests is to do so from inside the Python
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
  
* Further, Landlab-specific advice for developing your own components can be found
  in the :ref:`component development guide <dev_guide_components>`.


Testing
=======

Before merging any changes into the Landlab trunk, *all* unit tests (including
doctests) should be passing. In addition, any new features added to Landlab
should have an associated set of unit tests to verify that the new features
are working properly.

Landlab uses `Travis <https://travis-ci.org>`_ for continuous integration
testing. The `landlab page on Travis <https://travis-ci.org/landlab/landlab>`_
shows the latest testing results. A new set of tests are executed whenever
any changes are pushed to the Landlab repository and with every pull request.
We currently run test suites for Python versions 2.6, 2.7, 3.3, and 3.4 (for
Python 3 versions we use the `2to3 <https://docs.python.org/2/library/2to3.html>`_
tool).

Once you send a pull request from GitHub, you will be taken to the Landlab
pull request page and all unit tests are run. You will see the status
of the unit tests next to your latest commit description. If you see a green
check, all tests passed and your changes can be merged! However, if you see
an ex there was a problem running the tests. If you believe your changes are
responsible for the failures, please fix them until the tests pass. Note that
you do not need to send a new pull request after committing for fixes. They
will be added to the current pull request and the tests automatically rerun.

You can also run unit tests locally with `nose
<https://nose.readthedocs.org>`_. From the top-level Landlab folder (landlab with a lowercase *l* and the folder
that contains `setup.py`) run::

  > nosetests


Troubleshooting
===============

What do I do if my pull request cannot be automatically merged?
---------------------------------------------------------------

Get the latest upstream/master and go to the `master` branch. Remember, *do not develop here*.
Always develop in a feature branch. Merge the lastest upstream master with your master::
  > git fetch upstream
  > git checkout master
  > git merge upstream/master

Go to the branch on which you are developing and merge the lastest upstream master with your
branch::
  > git checkout <branch_name>
  > git merge upstream/master

Fix the conflicts. Do this by hand or with a merge editor. This is where you decide how to
integrate the conflicting changes. Since only you know what and why you made the changes
you did, this can only be done by you::
  > git mergetool

After everything has been fixed, commit the changes and push the changes to the repository.
The pull request will automatically be updated::
  > git commit
  > git push


I'm still confused
------------------

The Landlab development team will be happy to hear from you. Email one of us 
or create an issue request and we'll try to resolve your problem.

