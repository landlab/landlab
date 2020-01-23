.. _development:

====================
Guide for Developers
====================

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
  * Consider creating a small Pull Request to update the page on
    :ref:`ongoing development <ongoing_development>` so that others know what
    you are working on.

Supported Python Versions
-------------------------
Python 3.6, 3.7, and 3.8

If you need to introduce a new dependency, that dependency must be compatible
with Python  3.6+.

Installation, Contribution, and Development Practices
-----------------------------------------------------

.. toctree::
   :maxdepth: 2

   install/index
   contribution/index
   practices/index
   package_organization

Quick Links For Package Maintenance
-----------------------------------

There are a few pages that are particularly important for the maintenance of
the package. These are:

* :ref:`directory organization <organization>`
* :ref:`testing protocol <testing>`
* :ref:`dependency organization <dependencies>`, and
* :ref:`release workflow <dev_releases>`

If package maintainers change any of these, the prior pages likely need
updating.
