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

   git clone git://github.com/landlab/landlab.git

or download a `zip file <https://github.com/landlab/landlab/archive/refs/heads/master.zip>`_:

.. code-block:: bash

   curl -OL https://github.com/landlab/landlab/archive/refs/heads/master.zip

Once you have a copy of the source code, you can install it into your current
Python environment,

.. tab:: mamba

  .. code-block:: bash

     cd landlab
     mamba install --file=requirements.in
     pip install -e .

.. tab:: conda

  .. code-block:: bash

     cd landlab
     conda install --file=requirements.in
     pip install -e .

.. tab:: pip

  .. code-block:: bash

     cd landlab
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

The following list show how to use `nox`_ for some of the more common tasks:

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
