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

.. include:: ../../../README.rst
    :start-after: .. start-install-source
    :end-before: .. end-install-source

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
