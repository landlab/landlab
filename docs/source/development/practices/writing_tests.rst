.. _writing_tests:

================================================================
Writing docstring and unit tests for your component (or utility)
================================================================

All contributed code should be well tested. This should be done through both
doctests and standard unit tests using `pytest <https://docs.pytest.org/en/latest/>`_.

All public functions, classes, methods, etc. must have a docstring that follows
the `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_
conventions. Docstring tests should be short, easy-to-read tests that are
instructive to a user. These tests are included as examples in the Landlab
:ref:`Reference Manual <landlab>`.

Every ``.py`` file must contain a module-level docstring at the top of the file
that describes what the purpose of the file is.

Unit tests should be more extensive than doctests and give your new code
thorough testing. Ideally your tests will cover what happens within every
``if``, ``elif``, or ``else``, and every ``try`` or ``except`` block. These
test will also verify that if a bad parameter value is passed, that the correct
type of error is raised.

Additionally, unless there is a specific reason your component or utility can
only work with one Landlab grid type, the tests should verify that it can work
with multiple model grid types (verifying that the model works with both `
`HexModelGrid`` and ``RasterModelGrid`` is a good place to start).

Your unit tests should verify that the component or utility you are creating
does exactly what it is expected to do. This means you will probably want to
create a very small (e.g. 5x5) model grid and hand calculate what the correct
answer is. Then use *assertions* to ensure that your code reproduces that
answer exactly. It is very important to not just test that the code reproduces
the first answer you get. Instead you should construct the test so that you
**know** what the right answer is.

The `numpy testing <https://docs.scipy.org/doc/numpy-1.13.0/reference/routines.testing.html>`_
functions are useful for making comparison between actual and expected results
(e.g. to assert that one array is equal to another array). The
`pytest testing tools <https://docs.pytest.org/en/latest/assert.html>`_ are
useful for things like asserting that providing a particular set of values to
a function or component will raise a specific type of error.

Unit tests must be discoverable by `pytest <https://docs.pytest.org/en/latest/>`_.
This means that the unit tests should be in folders within the ``test``
directory, in ``.py`` files that start with the name `test`
and within functions with names that begin with the word ``test``.

Thus, a file to provide the unit tests for your component would be called
``test_my_component_name.py`` file, located in the directory
``tests\components\my_component_name\``. The inside of it might look like:

.. code-block:: python

    # test_my_component_name.py
    # numpy.testing and pytest are two modules commonly used
    # for testing whether your code behaves as expected.
    # import what you need from landlab


    def test_something_about_my_component():
        """Make a one-line docstring that describes your unit test."""
        # do things to set up for your test like make a model grid.

        # make your test and assert that you get the right answer
