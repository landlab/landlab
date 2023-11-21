.. _landlab_component_dev_page:

=====================================
Develop your own component or utility
=====================================

Landlab grows and improves thanks to user contributions. We encourage you to
develop your own component or utility!

Thank you to all who have contributed to landlab!

In addition to components, Landlab has many utilities useful in creating
components and doing model output post-processing. While this document
primarily discusses creating, documenting, and testing a components, most of
the steps are the same for developing a utility. Unlike a component, a utility
can just be a function. But like a component, we expect contributed utilities
to follow the landlab standard practices outlined in our documentation.

Once you have installed Landlab (:ref:`developer install <install>`)
and :ref:`created your own branch <landlab_develop_with_git>`, you can start
writing a Python script for your component.

See `this tutorial <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/making_components/making_components.ipynb>`_
for instructions on the structure and content of your component code. See also
`this example pull request <https://github.com/landlab/landlab/pull/678>`_, which
shows you the common set of files (such as ``__init__.py`` and documentation
files) that get added or modified when a component is added to Landlab.

See the :ref:`Standard Naming conventions <component_standard_names>` for good practice
on parameters and variables naming.

The following pages describe the software development practices that Landlab
strives to follow. Our goal is to make the capabilities of Landlab
well-documented to support new users while not enforcing substantial burdens on
community contributors. If you have any questions about the process after you
have finished reading the documentation, consider making an
`Issue <https://github.com/landlab/landlab/issues/>`_ to ask the
development team for help.


We recommend that you review the Landlab development practices:

.. toctree::
   :maxdepth: 2

   ../practices/index

Files structure
---------------
For your new component, you should create a folder in
``landlab/landlab/components`` caled ``<my_component_name>`` that contains:

- Your Python script `my_component_name.py`
- `_init_.py` which is structured as

  .. code-block:: python

      from .my_component_name import MyComponent

      __all__ = [
          "MyComponent",
      ]

  Where

      - `'.my_component_name'` is the name of the python script.
      - `'MyComponent'` is as defined in the _name header of your python script

  See `this tutorial <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/making_components/making_components.ipynb>`_
  on making a component for additional document requirements.

  In addition there are a number of recommendations and requirements for
  Landlab components :ref:`summarized here <dev_component_rules>`.

- a folder in the ``landlab/tests/components/<my_component_name>`` directory containing
  unit tests. The unit tests are run every time changes are pushed to the
  Landlab repository. They should go through every line of your code (e.g.
  test every possible scenario in if/else loops, exceptions, etc.). See
  `the tutorial on making a component <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/making_components/making_components.ipynb>`_
  for instructions about making docstring tests and the next section for more
  information about making the unit tests.

- a document in ``docs/source/reference/components`` called ``my_component_name.rst``.
  Look at other documents in that folder to get a sense of the typical format.
  This document is what will put your component's documentation on the
  ReadTheDocs page. See more below.

Once everything is working, you can :ref:`create a pull request <landlab_develop_with_git>`
to have your branch merged into the master so that your component can be
included in the Landlab library and used by others.

This will trigger :ref:`continuous integration testing <dev_ci>` of your branch
(doc tests, unit tests, and lint) to ensure its compatibility on all supported
environments. You can find the results of these tests on the GitHub page of
your pull request. If the tests fail, edit your files and commit your changes
to re-run the tests (you don't need to make another pull request).

Getting your component into the documentation
---------------------------------------------
Landlab uses the third party Sphinx code documentation tool to automatically
build the :ref:`Reference section <api>` that list our user-facing components.
This means your new component won't appear on the webpages unless you also make
some changes to files you'll find in `landlab/docs/source`.

You need to modify `landlab/docs/source/reference/components/index.rst`, and
also create a new file in the folder
`landlab/docs/source/reference/components`, called
`[short_name_for_your_component].rst`.

The best advice for both of these is to follow an existing example.

For the new `.rst` file, use e.g. `diffusion.rst` as a template. The first line
with the path specification needs to be changed to give the same name as the
`.py` file in which your component lives; the rest of the code text stays the
same.

For the update to `index.rst`, just copy what has been done for the others,
where the path specification now points at the new `.rst` file you made, i.e.,
`[short_name_for_your_component]` (leaving off the `.rst`).

Note your component won't appear on the user-facing part of the website until
it's included in a Landlab release.

Your component is accepted to Landlab. What's next?
---------------------------------------------------
Congrats on all your hard work! Once you know your component has been accepted
and is included in a Landlab release, please
`add it to the CSDMS Model Repository <https://csdms.colorado.edu/wiki/Contribute_model>`_.
You can link directly to the source code on Github. Just fill out the questionnaire.

Have you written up some tutorials or Jupyter notebooks to help teach new users
about your component? Consider submitting a tutorial along with your component.

If you've presented a poster or submitted a paper about your Landlab component,
advertise your work on the :ref:`Landlab Papers and Presentations <papers>`
page. Make your request to a member of the Landlab development team via a
GitHub Issue. If you'd like, also attach an abstract or poster PDF.

I'm still confused
------------------
The Landlab development team will be happy to hear from you.
`Create an issue <https://github.com/landlab/landlab/issues>`_ and we'll try to
resolve your problem.
