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

.. tab:: ssh

    .. code-block:: bash

       git clone git@github.com:landlab/landlab.git

.. tab:: https

    .. code-block:: bash

       git clone https://github.com/landlab/landlab.git

or download a `zip file <https://github.com/landlab/landlab/archive/refs/heads/master.zip>`_:

.. code-block:: bash

   curl -OL https://github.com/landlab/landlab/archive/refs/heads/master.zip

Once you have a copy of the source code, you can install it into your current
Python environment by first installing *Landlab* dependencies and then building
and installing *Landlab*.


Install dependencies
````````````````````

*Landlab*'s dependencies are listed in *requirements.in*.

.. tab:: mamba

  .. code-block:: bash

     cd landlab
     mamba install --file=requirements.in -c nodefaults -c conda-forge --override-channels

.. tab:: conda

  .. code-block:: bash

     cd landlab
     conda install --file=requirements.in -c nodefaults -c conda-forge --override-channels

.. tab:: pip

  .. code-block:: bash

     cd landlab
     pip install -r requirements.in


Build and install *Landlab*
```````````````````````````
*Landlab*'s build process includes compiling Python extensions, which requires
you to have a C++ compiler installed. *Linux* will usually already have one,
on *Mac* you can use *XCode*, and on *Windows* you will need to install *MSVC*.
For help on installing *MSVC*, you may want to refer to the *conda-forge* page
on `compiling code on Windows <https://conda-forge.org/docs/maintainer/knowledge_base.html#notes-on-native-code>`__
or the `Python wiki page for Windows compilers <https://wiki.python.org/moin/WindowsCompilers>`__.


If you are using *conda*/*mamba*, set up your compilers to build libraries
compatible with other installed packages,

.. tab:: mamba

  .. code-block:: bash

     mamba install compilers -c nodefaults -c conda-forge --override-channels

.. tab:: conda

  .. code-block:: bash

     conda install compilers -c nodefaults -c conda-forge --override-channels


With compilers set up and dependencies installed, build and install *Landlab*,

.. code-block:: bash

   pip install -e .



.. end-install-source

Developer Tools
---------------

Once you start developing with *Landlab*, we recommend that you use `nox`_  to
automate common tasks such as, for example, running the tests, building the docs, and
finding lint.

.. _nox: https://nox.thea.codes/en/stable/

.. code-block:: bash

  pip install nox

The following list shows how to use `nox`_ for some of the more common tasks:

* Run the tests:

  .. code-block:: bash

     nox -s test
* Run the tests on the notebooks:

  .. code-block:: bash

     nox -s test-notebooks
* Build the docs:

  .. code-block:: bash

     nox -s build-docs
* Run the linters:

  .. code-block:: bash

     nox -s lint
* To get a complete list of the available targets:

  .. code-block:: bash

     nox -l
