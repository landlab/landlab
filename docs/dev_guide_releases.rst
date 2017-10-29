.. _dev_releases:

===============================
How to create a Landlab release
===============================

New releases are built and uploaded to
`Anaconda.org <https://anaconda.org/landlab/landlab>`_ whenever a new tag
that starts with the letter ``v`` is
`created and pushed to <https://git-scm.com/book/en/v2/Git-Basics-Tagging>`_
`GitHub <https://github.com/landlab/landlab>`_. As an example, the following
will cause a new release to be built:

.. code-block:: bash

    $ git tag v0.1.1 # Create the tag locally
    $ git push --tags # Push the tag to the remote

A new release is created (``v0.1.1``) and the tag pushed to GitHub.
`Travis-CI <https://travis-ci.org/landlab/landlab>`_ notices the tagged commit,
and after building and testing the package, creates a fresh new package that
is uploaded to `Anaconda.org <https://anaconda.org/landlab/landlab>`_.

.. note::

  Although you can create such a tag on any branch, releases should **only**
  come from the ``release`` branch. Make sure that when you create a tag
  you are doing so on ``release`` (and all your changes are committed).

A couple notes about creating a new version:

1. Landlab follows `Semantic Versioning <http://semver.org/>`_
   rules for version assignment and formatting. Please stick to them.

2. The version given in the tag name must match that in
   ``landlab/__init__.py``. The version must also be changed in
   ``.conda-recipe/meta.yaml``. However, ``meta.yaml`` doesn't like dashes
   so if your version contains a dash just leave it out in this file
   (for example ``1.0.0-beta.6`` becomes ``1.0.0beta.6``).

3. If you mess up (forget to update all the version strings scattered
   throughout the code, for example), you can always `delete the tag and
   recreate it <https://git-scm.com/docs/git-tag>`_. To do this, you'll
   need to delete both the remote tag and the local tag.

   .. code-block:: bash

      $ git push --delete origin <tagname> # Delete the tag on the remote repository
      $ git tag --delete <tagname> # Delete the tag from the local repository

   where ``<tagname>`` is the name of your tag (``v0.1.1``, for example).

4. If your new tag was successfully pushed to GitHub, you will be able to see
   it with the rest of the
   `releases <https://github.com/landlab/landlab/releases>`_ and
   `tags <https://github.com/landlab/landlab/tags>`_.

5. To see if your new release was created successfully, you can do one or all
   of the following:

   *  Check the logs for the build of your tagged commit on
      `Travis-CI <https://travis-ci.org/landlab/landlab>`_.
   *  Check `Anaconda.org <https://anaconda.org/landlab/landlab>`_ to see
      if your release appears there.
   *  Check if `conda` can see your new release with
      ``conda search landlab -c landlab``. See the
      `conda docs <http://conda.pydata.org/docs/using/index.html>`_
      for a description of ``conda`` and how to use it, or you can always use
      ``conda -h`` from the command line.

The Release Checklist
=====================
1. Make sure you are on the ``release`` branch.

   .. code-block:: bash

      $ git checkout release
2. Make sure all the version strings match and use
   `Semantic Versioning <http://semver.org/>`_.

   *  ``landlab/__init__.py``
   *  ``.conda-recipe/meta.yaml``
3. Commit your changes.
4. Create a tag for this release that matches the string in ``__init__.py``
   but that starts with the letter ``v``.

   .. code-block:: bash

      $ git tag v0.1.1
5. Push your tag to the remote.

   .. code-block:: bash

      $ git push --tags

Helpful links
=============

1. `Using conda <http://conda.pydata.org/docs/using/index.html>`_: What
   `conda` is and how to use it.
2. `git tags <https://git-scm.com/book/en/v2/Git-Basics-Tagging>`_: What git
   tags are and how to create them.
3. `The git tag command <https://git-scm.com/docs/git-tag>`_: A description
   of all of the options for the `git tag` command (including `git tag
   --delete`).
4. `landlab on Travis <https://travis-ci.org/landlab/landlab>`_: The latest
   Travis builds of landlab.
5. `landlab on Anaconda <https://anaconda.org/landlab/landlab>`_: The
   conda packages for landlab releases.
6. `Semantic Versioning <http://semver.org/>`_: Rules for assigning and
   formatting versions.
