.. _install:

=================
Developer Install
=================

.. important::

  The following commands will install *landlab* into your current environment. Although
  not necessary, we **highly recommend** you install landlab into its own
  :ref:`virtual environment <virtual_environments>`.


If you will be modifying code or contributing new code to *landlab*, you will first
need to get *landlab*'s source code and then install *landlab* from that code.

Source Install
--------------

.. start-install-source

*Landlab* is actively being developed on GitHub, where the code is freely available.
If you would like to modify or contribute code, you can either clone our
repository

.. code-block:: bash

   $ git clone git://github.com/landlab/landlab.git

or download the `tarball <https://github.com/landlab/landlab/tarball/master>`_
(a zip file is available for Windows users):

.. code-block:: bash

   $ curl -OL https://github.com/landlab/landlab/tarball/master

Once you have a copy of the source code, you can install it into your current
Python environment,

.. tab:: conda

  .. code-block:: bash

     $ cd landlab
     $ conda install --file=requirements.txt
     $ pip install -e .

.. tab:: pip

  .. code-block:: bash

     $ cd landlab
     $ pip install -e .

.. end-install-source

Developer Tools
---------------

Once you start developing with *landlab*, there are a number of other packages you
may find useful to install. These packages are used for, among other things,
testing *landlab*, and ensuring your code complies with *landlab*'s development
standards.

.. tab:: conda

  .. code-block:: bash

    $ conda install --file requirements-dev.txt --file requirements-testing.txt

.. tab:: pip

  .. code-block:: bash

    $ pip install -e ".[dev,testing]"
