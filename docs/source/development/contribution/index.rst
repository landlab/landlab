.. _dev_contributing:

===========================================
Guidelines for Contributing Code to Landlab
===========================================

Please review the :ref:`development practices <development_practices>`.

.. toctree::
   :maxdepth: 2

   develop_a_component
   recommendations
   ongoing_development
   desired_contributions

Contributions that change high level Landlab organization should update
:ref:`this page <organization>` in the documentation.

Publication of Landlab Contributions
------------------------------------

A number of researchers have used the Journal of Open Source Software (JOSS)
and Geoscientific Model Development as a journal outlet to publish their
contribution to Landlab. There are a few aspects of the submission workflow for
JOSS that are non-standard and are described below.

.. toctree::
   :maxdepth: 2

   joss_workflow

In addition, the JOSS editorial board and the Landlab core development team are
presently (Jan 2019) working on defining some guidelines regarding what
contributions to Landlab fit the scope of a JOSS submission. These will be
summarized here when finalized. If you have any questions regarding whether
your potential JOSS submission is appropriate, the best thing to do is to
make a pre-submission inquiry with JOSS.

News Entries
------------

The ``CHANGES.rst`` file is managed using `towncrier`_ and all non trivial changes
must be accompanied by a news entry.

To add an entry to the news file, first you need to have created an issue
describing the change you want to make. A Pull Request itself *may* function as
such, but it is preferred to have a dedicated issue (for example, in case the
PR ends up rejected due to code quality reasons).

Once you have an issue or pull request, you take the number and you create a
file inside of the ``news/`` directory, named after that issue number with a
"type" of ``component``, ``notebook``, ``feature``, ``bugfix``, ``docs``, or ``misc``
associated with it.

If your issue or PR number is ``1234`` and this change is fixing a bug,
then you would create a file ``news/1234.bugfix.rst``. PRs can span multiple
categories by creating multiple files (for instance, if you added a new component
and an associated notebook that demonstrates how to use it, you would create
``news/NNNN.component.rst`` and ``news/NNNN.notebook.rst``).

If a PR touches multiple issues/PRs, you may create a file for each of them
with the exact same contents and Towncrier will deduplicate them.


.. _`towncrier`: https://pypi.org/project/towncrier/
