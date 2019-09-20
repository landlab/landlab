`Landlab <http://landlab.github.io>`__ \| [[About \| About]] \|
[[Examples \| Examples]] \| [[User Guide \| User-Guide]] \| `Reference
Manual <http://landlab.readthedocs.org/en/latest/#developer-documentation>`__
\| [[Tutorials\| Tutorials ]] \| [[FAQs \|FAQs]]

.. raw:: html

   <p>

This should work for Anaconda users with Windows 7+, Mac OS 10.6+, or
Ubuntu Linux (only the latest version has been tested).

.. raw:: html

   </p>

.. raw:: html

   <p>

Once you have a full Python distribution on your machine, it is vital to
check that it has been successfully set as the default copy of Python on
your system. Open a command prompt (Terminal on a Mac, or Command Prompt
on a PC) and type the lines below (note the > indicates that you are on
a command line):

.. raw:: html

   </p>

MAC:

.. raw:: html

   <p>

 > which python

.. raw:: html

   </p>

.. raw:: html

   <p>

 > which ipython

.. raw:: html

   </p>

.. raw:: html

   <p>

WINDOWS:

.. raw:: html

   <p>

 > where python

.. raw:: html

   </p>

.. raw:: html

   <p>

 > where ipython

.. raw:: html

   </p>

.. raw:: html

   <p>

In each case, both commands should return the same path, and it should
clearly refer to Anaconda (or Canopy). Details will depend on your
operating system but it could look something like this:

.. raw:: html

   </p>

.. raw:: html

   <p>

 /anaconda/bin/python

.. raw:: html

   </p>

.. raw:: html

   <p>

If you don’t see reference to your newly installed distribution (i.e.,
/anaconda), click here to resolve the problem.

.. raw:: html

   </p>

.. raw:: html

   </p>

Make sure you have the latest version installed (close anaconda before
doing this):

.. raw:: html

   </p>

 conda update –all

.. raw:: html

   <p>

.. raw:: html

   <p>

Once the path to both python and ipython point to your new distribution,
open the Python editor in Anaconda called Spyder.

.. raw:: html

   </p>

.. raw:: html

   <p>

On the Spyder toolbar, go to Tools → Open command prompt to open the
command line.

.. raw:: html

   <p>

.. raw:: html

   </p>

Alternatively you can open a standard terminal window, such as an xterm
(X11.app) or terminal window (Terminal.app) on a Mac, or a command
prompt on a Windows machine. If you do use a standard terminal and run
into problems, make sure you have resolved your path issues.

.. raw:: html

   </p>

.. raw:: html

   <h3>

Now to install Landlab!

.. raw:: html

   </h3>

.. raw:: html

   <p>

You can either install with the conda or the pip package managers. Conda
is recommended, as it reduces the chances of versioning conflicts. Try
to remember which you choose to avoid confusion when updating later! (If
you installed landlab prior to May 19th 2016, you will have used pip).

.. raw:: html

   </p>

.. raw:: html

   <p>

Type either (for conda install):

.. raw:: html

   <p>

 > conda install landlab -c landlab -c conda-forge

.. raw:: html

   </p>

.. raw:: html

   </p>

.. raw:: html

   <p>

…or (for pip install, not recommended for entry level users):

.. raw:: html

   <p>

 > pip install landlab

.. raw:: html

   </p>

.. raw:: html

   </p>

.. raw:: html

   <h3>

Test Landlab install

.. raw:: html

   </h3>

.. raw:: html

   <p>

Once Landlab has been successfully installed, on the Python shell line,
check to make sure it is up-to-date (note that those are double
underscores around version; also note that you may need to close and
reopen Anaconda before typing the below commands):

.. raw:: html

   </p>

.. raw:: html

   <p>

 > import landlab

.. raw:: html

   </p>

.. raw:: html

   <p>

 > landlab.__version_\_

.. raw:: html

   </p>

.. raw:: html

   <p>

The version number should be greater than 1. You can check the version
number of the most recent release here.

.. raw:: html

   </p>

.. raw:: html

   <h3>

Install/Test problems

.. raw:: html

   </h3>

.. raw:: html

   <p>

If you are having problems when installing, testing or running Landlab,
please visit our Troubleshooting page.

.. raw:: html

   </p>
