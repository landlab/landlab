.. _basic_install:

=======
Install
=======

.. important::

  The following commands will install *landlab* into your current environment. Although
  not necessary, we **highly recommend** you install landlab into its own
  :ref:`virtual environment <virtual_environments>`.

.. start-install-release

In order to use *landlab* you will first need Python. While not
necessary, we recommend using the
`Anaconda Python distribution <https://www.anaconda.com/distribution/>`_
as it provides a large number of third-party packages useful for
scientific computing.

To install *Landlab*, simply run the following in your terminal of choice:

.. tab:: mamba

  .. code-block:: bash

    mamba install landlab -c nodefaults -c conda-forge --override-channels

.. tab:: conda

  .. code-block:: bash

    conda install landlab -c nodefaults -c conda-forge --override-channels

.. tab:: pip

  .. code-block:: bash

    pip install landlab

.. end-install-release

If you would like the very latest development version of *landlab* or want to modify
or contribute code to the *landlab* project, you will need to do a
:ref:`developer installation <install>` of *landlab* from source.
