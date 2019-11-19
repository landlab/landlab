.. _rough_guide:

Rough Guide to Python and Landlab Install
==========================================

*Handling Python installs on Linux if you are non-expert can be a confusing
and frustrating experience. However, one of us (DEJH) has recently gone
through the pain, and offers these thoughts on how to get a Python install
capable of running Landlab onto your machine. (As long as your machine is
Debian...)*

NB: A version of Canopy does exist for linux, but DEJH couldn't get this
to install properly. If you can, this is probably the easier way to go!


What you already have, apt-get, and versioning
----------------------------------------------

Your unix box almost certainly already has a version on Python on it.
However, DEJH found it wasn't the right version in his case. The problem
is that the stable Python (and Python package) versions that are
maintained through the Debian package repository that *apt-get* looks at
tend not to be up to date. This is problematic, as Landlab needs at least
Python 2.7, Numpy 1.8, and Scipy 0.12.
*NOTE: this is very very old. As of Jan 2020, Landlab will require Python 3 or
higher*.
So first, check your version:

.. code-block:: bash

    $ python --version

Note the two dashes.
If you got 2.7 or higher, congratulations! You win. No further changes to
python installs are needed. But if you didn't, first check whether you can
simply update through the package manager:

.. code-block:: bash

    $ sudo apt-get install -u python
    $ python --version

If you're still seeing a sub-2.7 version number (as I was...), you're going
to need to download the python source and install it manually. Go to
`the python download site <https://www.python.org/downloads/>`_, select the
link to the highest version
number of Python 2.X available (2.7.6 at time of writing), and download it
as a gzipped source tarball. Then:

.. code-block:: bash

    $ tar -xzf the_downloaded_file.tgz
    $ cd the_new_directory
    $ ./configure
    $ make
    $ sudo make install

You should now have Python 2.7! Close the terminal, reopen, and check:

.. code-block:: bash

    $ python2.7

BUT, it's highly likely if you just run "python", you won't get this
version. Check by

.. code-block:: bash

    $ python --version

again. You can resolve this by messing with your ~/.bashrc profile. Open
that file in your favourite text editor, and add at the end:

.. code-block:: bash

    $ alias python="python2.7"

Restart the terminal again, and check the version number now. Should be
right! (You will probably also want to do something similar with ipython,
but DEJH didn't explore that...) 

Installing pip
--------------

You now have the right version of python, but if you try, e.g., *import
numpy*, you'll notice you don't have any packages. This is where pip -
the powerful command line Python package manager - comes in.

Get pip by going to `the pip site
<https://pip.pypa.io/en/latest/installing/>`_
and downloading get_pip.py, which links from that page. Navigate to
the folder you downloaded it into, and simply run

.. code-block:: bash

    $ sudo python get_pip.py

This should give you a trouble free install of pip.

Once you have it, make sure you're fully up-to-date:

.. code-block:: bash

    $ pip install --upgrade pip

NB: DO NOT TRY TO INSTALL PIP WITH APT-GET. Pip binds to your Python
install, and the binding probably won't take properly unless you
install through your version of Python, as described here. Note that
there is a copy of pip you can get with apt-get, but you don't
want it.


Downloading the packages
------------------------

Now you have pip and it's bound correctly to your Python install,
adding packages should be trouble free:

.. code-block:: bash

    $ sudo pip install numpy
    $ sudo pip install scipy
    $ sudo pip install matplotlib
    $ sudo pip install sympy
    $ sudo pip install netCDF4

Note in future, you can update these packages to new versions by:

.. code-block:: bash

    $ sudo pip install --upgrade [package_name]

Now test the versions like this:

.. code-block:: bash

    $ python
    >>> import numpy
    >>> numpy.__version__

And everything should now be great. You can now continue to install
Landlab as you would in the main instructions. e.g., if you have
a clone or downloaded copy of Landlab you want to install in
developer mode, just navigate to the download's top level directory
and run

.. code-block:: bash

    $ python setup.py develop

and test:

.. code-block:: bash

    $ python
    >>> import landlab
    >>> landlab.test()

Or alternatively, just grab the release version using pip, as in the
main instructions:

.. code-block:: bash

    $ pip install landlab

& again, test as above.
