.. _landlab_component_dev_page:

Develop your own component or utility
-------------------------------------


Landlab grows and improves thanks to user contributions. We encourage you to develop your own component or utility!

Thank you to all who have contributed to landlab!
:heart::tada::heart::heart::tada::heart::tada::heart::heart::tada::heart::

In addition to components, Landlab has many utilities useful in creating components and doing model output post-processing. While this document primarily discusses creating, documenting, and testing a components, most of the steps are the same for developing a utility. Unlike a component, a utility can just be a function. But like a component, we expect contributed utilities to follow the landlab standard practices outlined here.

Once you have installed Landlab (`developer install <https://github.com/landlab/landlab/wiki/Installing-Landlab-from-source-code-(%22developer-install%22)>`_) and `created your own branch <https://github.com/landlab/landlab/wiki/Developing-with-github-and-git>`_, you can start writing a Python script for your component.

See `this tutorial <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/making_components/making_components.ipynb>`_ for instructions on the structure and content of your component code. See also `this example pull request <https://github.com/landlab/landlab/pull/678>`_, which shows you the common set of files (such as ``__init__.py`` and documentation files) that get added or modified when a component is added to Landlab.

See `Standard Naming conventions <https://github.com/landlab/landlab/wiki/Standard-names>`_ for good practice on parameters and variables naming.

If you have any questions about the process after you have finished reading the documentation, consider making an `Issue
<https://github.com/landlab/landlab/issues/new/>`_ to ask the development team for help.

Recommendations on coding style
-------------------------------

See `this page <https://github.com/landlab/landlab/wiki/Style-conventions#style-enforcement>`_ for recommendations on coding style.

- Please stick to the coding style described by `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_.
- Class and function docstrings should follow the `numpydoc conventions <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
- Further, Landlab-specific advice for developing your own components can be found in the `tutorial <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/making_components/making_components.ipynb>`_.

If you want to check how well you are doing, please look at our `Landscape page <https://landscape.io>`_. Landscape will grade the health of the landlab code with every push to GitHub.

Files structure
---------------
For your new component, you should create a folder in landlab/landlab/components that contains:

- Your Python script `my_component_name.py`
- `_init_.py` which is structured as::


   from .my_component_name import MyComponent

   __all__ = ['MyComponent', ]


`‘.my_component_name’` is the name of the python script.
`‘MyComponent’` is as defined in the _name header of your python script

- a folder called ‘tests’ containing unit tests. The unit tests are run every time changes are pushed to the Landlab repository. They should go through every line of your code (e.g. test every possible scenario in if/else loops, exceptions, etc.). See `the tutorial on making a component <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/making_components/making_components.ipynb>`_ for instructions about making docstring tests and the next section for more information about making the unit tests.

Once everything is working, you can `create a pull request <https://github.com/landlab/landlab/wiki/Developing-with-github-and-git>`_ to have your branch merged into the master so that your component can be included in the Landlab library and used by others.

This will trigger integration `testing <http://landlab.readthedocs.io/en/latest/dev_guide_install.html#testing>`_ of your branch (doc tests and unit tests) to ensure its compatibility on all supported environments. You can find the results of these tests on the GitHub page of your pull request. If the tests fail, edit your files and commit your changes to re-run the tests (you don’t need to make another pull request).

Writing docstring and unit tests for your component (or utility)
----------------------------------------------------------------

All contributed code should be well tested. This should be done through both doctests and standard unit tests using `pytest <https://docs.pytest.org/en/latest/>`_.

Doctests should be short, easy-to-read tests that are instructive to a user. These tests are included as examples in the Landlab `Reference Manual <http://landlab.readthedocs.io/en/release/>`_.

Unit tests should be more extensive than doctests and give your new code thorough testing. Ideally your tests will cover what happens within every ``if``, ``elif``, or ``else``, and every ``try`` or ``except`` block. These test will also verify that if a bad parameter value is passed, that the correct type of error is raised.

Additionally, unless there is a specific reason your component or utility can only work with one Landlab grid type, the tests should verify that it can work with multiple model grid types (verifying that the model works with both ``HexModelGrid`` and ``RasterModelGrid`` is a good place to start).

Your unit tests should verify that the component or utility you are creating does exactly what it is expected to do. This means you will probably want to create a very small (e.g. 5x5) model grid and hand calculate what the correct answer is. Then use *assertions* to ensure that your code reproduces that answer exactly. It is very important to not just test that the code reproduces the first answer you get. Instead you should construct the test so that you **know** what the right answer is.

The `numpy testing <https://docs.scipy.org/doc/numpy-1.13.0/reference/routines.testing.html>`_ functions are useful for making comparison between actual and expected results (e.g. to assert that one array is equal to another array). The `pytest testing tools <https://docs.pytest.org/en/latest/assert.html>`_ are useful for things like asserting that providing a particular set of values to a function or component will raise a specific type of error.

Unit tests must be discoverable by `pytest <https://docs.pytest.org/en/latest/>`_. This means that the unit tests should be in folders called ``test`` within the component or utility folder, in ``.py`` files that start with the name `test` and within functions with names that begin with the word ``test``.

Thus, a file to provide the unit tests for your component would be called ``test_my_component_name.py`` file, located in the directory ``landlab\components_my_component_name\tests\``. The inside of it might look like::

   # test_my_component_name.py
   # numpy.testing and pytest are two modules commonly used
   # for testing whether your code behaves as expected.
   # import what you need from landlab

   def test_something_about_my_component():
       """Make a one-line docstring that describes your unit test."""
       # do things to set up for your test like make a model grid.

       # make your test and assert that you get the right answer

Getting your component into the documentation
---------------------------------------------
Landlab uses the third party Sphinx code documentation tool to automatically build the Reference Manual webpages that list our user-facing components. This means your new component won't appear on the webpages unless you also make some changes to files you'll find in `landlab/docs`.

You need to modify `index.rst`, and also create a new file in the folder, called `landlab.components.[short_name_for_your_component].rst`.

The best advice for both of these is to follow an existing example.

For the new `.rst` file, use e.g. `landlab.components.diffusion.rst` as a template. The first line with the path specification needs to be changed to give the same name as the `.py` file in which your component lives; the rest of the code text stays the same.

For the update to `index.rst`, just copy what has been done for the others, where the path specification now points at the new `.rst` file you made, i.e., `landlab.components.[short_name_for_your_component]` (leaving off the `.rst`).

Note your component won't appear on the user-facing part of the website until it's included in a Landlab release.

Your component is accepted to Landlab. What's next?
---------------------------------------------------
Congrats on all your hard work! Once you know your component has been accepted and is included in a Landlab release, please `add it to the CSDMS Model Repository <http://csdms.colorado.edu/wiki/Contribute_model>`_. You can link directly to the source code on Github. Just fill out the questionnaire.

Have you written up some tutorials or Jupyter notebooks to help teach new users about your component? Consider submitting to the `Landlab tutorials <https://github.com/landlab/tutorials>`_ repository. Contact a Landlab developer for more information on how to share your teaching tools.

If you've presented a poster or submitted a paper about your Landlab component, advertise your work on the `Landlab Papers and Presentations <https://github.com/landlab/landlab/wiki/Landlab-Papers-and-Presentations>`_ wiki page. E-mail your request to a member of the Landlab development team. If you'd like, also attach an abstract or poster PDF.

I’m still confused
------------------
The Landlab development team will be happy to hear from you. Email one of us or `create an issue <https://github.com/landlab/landlab/issues>`_ and we’ll try to resolve your problem.
