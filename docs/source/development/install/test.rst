.. _testing:

=========================================
Verify that Landlab is Working with Tests
=========================================

Testing Your Installation
-------------------------

In order to test your installation you'll need to install the
`pytest <https://docs.pytest.org/en/latest/>`_ package that is used to run the
tests. This is included in the conda environment in the
:ref:`developer_install <dev_install_install>`.

Alternatively, you can install pytest with conda:

.. code-block:: bash

   $ conda install pytest

Once ``pytest`` has been installed navigate to the main Landlab
directory (the one with ``setup.py`` in it) and type into a terminal:

.. code-block:: bash

   $ pytest

This command will collect and run all of the tests.

More Advanced Testing Options
`````````````````````````````

If you want to only want to test one part of Landlab:

.. code-block:: bash

   $ pytest path/to/directory/you/want/to/test

Landlab's code (and thus docstring tests) are located in the directory
``landlab/landlab`` and the Landlab unit tests are located in the directory
``landlab/tests``. The directory structure is the same for both directories, so
to run the doctests and the unit tests for a particular part of landlab, one
can use the ``*`` wildcard. For example, to test one particular component:

.. code-block:: bash

   $ pytest */components/<name of component directory>

Testing the Notebooks
`````````````````````

Because the Landlab notebooks take a lot longer to test, the default
configuration of pytest does not evaluate them. To test the entire package,
including the notebooks, evaluate:

.. code-block:: bash

   $ pytest --run-notebook

Testing with Coverage
`````````````````````

You may also want to see the code coverage of different parts of
Landlab. To do this, you'll first need to install
`pytest-cov <https://pytest-cov.readthedocs.io/en/latest/readme.html>`_.

.. code-block:: bash

   $ conda install pytest-cov

Then execute

.. code-block:: bash

   $ pytest landlab --doctest-modules --cov=landlab --cov-report term-missing

from the main Landlab directory.

This will run the tests and print the coverage statistics (including the
missing line numbers) to the terminal.

As above, you can also run the coverage tools for a more specific
directory. For example, to run them for your current directory you could
execute

.. code-block:: bash

   $ pytest . --doctest-modules --cov=. --cov-report term-missing

or to run them for a specific directory (for example, the
erosion_deposition submodule) stored as an environment variable you
would do the following:

.. code-block:: bash

   $ TEST_DIR=landlab/components/erosion_deposition/
   $ pytest $TEST_DIR --doctest-modules --cov=$TEST_DIR --cov-report term-missing
