.. _dependencies:

==================================
How Landlab Specifies Dependencies
==================================

The Landlab package dependencies are found in the top-level directory file
``requirements.txt``.

The ``setup.py`` file looks to this file and reads it when compiling the
package. In this way, we only state the dependencies in one location.

In addition to the core set of dependencies, development with Landlab, and
common use of Landlab (e.g., running the notebooks) may have additional
dependencies. These dependencies are described in the following files:

- ``requirements-notebooks.txt`` indicates dependencies for running the notebooks.
- ``requirements-testing.txt`` indicates dependencies for running tests.
- ``requirements-dev.txt`` indicates dependencies for development (less building the documentation).

We provide two convenience environment files that address the two most common
use cases.

- ``environment.yml`` specifies an environment which installs a Landlab binary
  and all of the notebook dependencies.
- ``environment-dev.yml`` specifies an environment with all development
  dependencies.

If a developer wants to build the documentation locally there are some
additional dependencies. We do not include these in ``requirements-dev.txt``
because it is uncommon for developers to need to build the docs locally.

An environment specifying documentation-building requirements can be found at
``landlab/docs/environment.yml``.
