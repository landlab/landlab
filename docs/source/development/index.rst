.. _development:

===============
Developer Guide
===============

If you're intending to make changes to the Landlab code base, or want to
develop your own components, this set of pages will provide information you
need.

Key Development Recommendations
-------------------------------

The most important things to remember are to:

* Start by creating a fork of Landlab's repository;
* Make changes to the source code on a development branch, not the default
  `master` branch; and
* Keep your fork's `master` and development branches up to date with changes
  in the main Landlab repository.

Supported Python Versions
-------------------------
Python 3.10, 3.11, and 3.12

If you need to introduce a new dependency, that dependency must be compatible
with Python 3.10+ and be available on Linux, Mac, and Windows.

.. toctree::
   :maxdepth: 2
   :hidden:

   contribution/index
   practices/index
   package_organization

Quick Links For Package Maintenance
-----------------------------------

There are a few pages that are particularly important for the maintenance of
the package. These are:

* :ref:`directory organization <organization>`
* :ref:`dependency organization <dependencies>`, and
* :ref:`release workflow <dev_releases>`

If package maintainers change any of these, the prior pages likely need
updating.
