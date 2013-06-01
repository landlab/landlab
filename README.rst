===========
The Landlab
===========

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

Install
=======

Installing a release
--------------------

We don't have any releases of The Landlab yet so this doesn't work righ now.
However, once we have a stable release, this is how users will install our
software on their machine::

    $ pip install TheLandlab

This will install Landlab as well as any prerequisite packages (required packages
are listed in setup.py).

Installing from a source distribution
-------------------------------------

To install the landlab package from a source distribution (like what you get with
a Subversion checkout), run the setup.py script with the install argument::

    $ python setup.py install

This will put the landlab package into a system-wide location so that you can
import landlab from any python session. It will also install any required
packages.

Note that if you are developing landlab you will probably not want to do this!
If you do, you would have to run this command everytime you make a change to
the code. Instead, you want to install the package in "development mode". See
below in the "Development Environment" section for details on how to do this.


----------------------
Developing The Landlab
----------------------

This section describes one possible workflow when developing The Landlab.

Using Subversion
================

The standard resource for Subversion "Version Control with Subversion" book.
It's available online and is completely *free*!

http://svnbook.red-bean.com/


Checkout a copy of the source
-----------------------------

Use Subversion to checkout the latest version of the source code::

    $ svn co https://csdms.colorado.edu/svn/TheLandlab/trunk landlab

Alternatively, if you already have a working copy of the repository you may want
to update it with changes that others have committed to the repository::

    $ svn update

Review your changes
-------------------

Now that you've made changes to some files, you'll probably want to have a
look at what you've done. To see what it is you've changed since your last
commit (or since you last checkout), use the status command::

    $ svn status

Use the Subversion diff command to see the differences between your copy and the
copy that you checked out::

    $ svn diff <changed_file>

where *<changed_file>* is the name of the file you've made changes to. If you
don't list any file names, this command will print the differences for *all*
files that have changes.

Commit your changes to the repository
-------------------------------------

Once you have made changes to your copy of the source and are happy with the
changes, you can commit them back to the repository::

    $ svn commit

This will commit all changes that you have made under the current directory. If
you only want to commit changes to a file or two, you can list the file names
separately on the command line::

    $ svn commit <one_file> <another_file>


---------------------------
The Development Environment
---------------------------

Once I have a working copy of The Landlab source code, I use the pip command to
install a development version of the code. If I'm in the base landlab folder
(the folder that contains setup.py), I run the following::

    $ pip install -e .

This sets up python so that it knows where the landlab package is when try to
import it - regardless of what directory you are in. This allows python commands
like::

    >>> import landlab
    >>> from landlab import craters

to work. If you didn't do this you might start getting errors that contain 
something like::

    ImportError: No module named landlab

To uninstall your development version of landlab::

    $ pip uninstall TheLandlab


------------------
Running Unit Tests
------------------

Immediatly after update your working copy of the code (or checking out a new
version) I will normally run the unit tests for the package to make sure nothing
is broken. You can do this with setup.py::

    $ python setup.py test

You should also probably do this before commiting changes to the repository to
make sure you didn't break things.


------------
Coding Style
------------

Because Python is so flexible style-wise, please try to stick to the coding
style described by PEP8,

http://www.python.org/dev/peps/pep-0008/

An easy way to make sure that you've done this is by running the pep8 command
on each file that you edit. If you don't have pep8 installed, you will have to
install it with::

    $ pip install pep8

Now you can run it on a Python source file. For instance::

    $ pep8 craters.py

At first, this will probably return lots of problems with you source file. Don't
worry though, it won't take long to get used to the coding style and be able to
write compatible code straigt away. If we stick to this it will make it much
easier to read the code written by any one of us.


-----------------------
Build API Documentation
-----------------------

You can build documentation for the LandLab API using Sphinx. Once you have set
up your envrionment to properly import landlab, you can generate the necessary
sphinx files with::

    $ python setup.py build_sphinx

This will put a bunch of files in the docs folder. The HTML documentation will
be under the docs/_build/html/ folder. Pointing your browser to index.html
under this folder will give you the top-level page for the documentation. This
entire folder is relocatable, so if you would like your documentation elsewhere
you can easily move the folder around.

If you have added, removed, or renamed files you may need to regenerate some of
the sphinx files and rebuild the api docs. If you have Sphinx installed, you can
do this with (from the directory that contains setup.py)::

    $ sphinx-apidoc -o docs landlab


Happy Landlab-ing!
==================
