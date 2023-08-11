.. _dev_ci:

================================
Continuous Integration Practices
================================

Before merging any changes into the Landlab trunk, *all* unit tests (including
doctests) should be passing. In addition, any new features added to Landlab
should have an associated set of unit tests to verify that the new features
are working properly.

Landlab uses `Travis <https://travis-ci.org>`_ for continuous integration
testing on OSX and Linux. The `landlab page on Travis <https://travis-ci.org/landlab/landlab>`_
shows the latest testing results. A new set of tests are executed whenever
any changes are pushed to the Landlab repository and with every pull request.
We currently run test suites for Python versions 3.6 and 3.7 (soon to include 3.8).

Continuous integration for Windows is done on
`Appveyor <https://ci.appveyor.com>`_ and also tests for the same
Python versions as OSX and Linux.

Once you send a pull request from GitHub, you will be taken to the Landlab
pull request page and all unit tests are run. You will see the status
of the unit tests next to your latest commit description. If you see a green
check, all tests passed and your changes can be merged! However, if you see
an ex there was a problem running the tests. If you believe your changes are
responsible for the failures, please fix them until the tests pass. Note that
you do not need to send a new pull request after committing for fixes. They
will be added to the current pull request and the tests automatically rerun.

You can also run unit tests locally in the top `landlab` directory and typing

.. code-block:: bash

    $ pytest

Note that this will test whatever version of landlab you have installed,
which may or may not be the one you are working on in your current working
directory. These test will not work with numpy 1.14.

TODO: Add information about coverage here.
