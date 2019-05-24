#! /usr/bin/env python
"""
Unit tests for ComponentMaker
"""
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.components import ComponentMaker

# nose.tools and numpy.testing are both good resources for functions that help
# with unit tests. This is not all of the functions that are avaliable, but a
# good smattering of functions to start with. Read their docs!


# it is important that each function that is tested starts with the word test
# it is important that this file has a name that starts with the word test
# it is important this file is in a folder that starts with the word test.
# these things help nose find the tests!

# Ideally all options (e.g. if/elif/else blocks) are tested by your unit tests
# and docstring tests.


def test_simple():
    """This is a statment about what the test does. It should be short and informative."""
    pass


def test_spam_type():
    """Test that passing a bad variable type to spam raises a ValueError"""
    mg = RasterModelGrid(3, 3)
    with pytest.raises(ValueError):
        ComponentMaker(mg, spam=1.0, eggs=1.0)


def test_eggs_type():
    """Test that passing a bad variable type to eggs raises a ValueError"""
    mg = RasterModelGrid(3, 3)
    with pytest.raises(ValueError):
        ComponentMaker(mg, spam=True, eggs=False)


def test_that_component_works_with_non_raster_grids():
    """Test that your component can work with non raster grids like HexModelGrid."""
    # This is especially important if your component behaves differently when
    # flow directing is done with our without the raster diagonals (called d8s
    # in landlab).
    pass


def test_correct_solution_found():
    """Test that your component creates the correct solution."""
    # Ideally you component (or functions used by your component) have known
    # behavior. That is, given a set of inputs, you know what the output is.

    # You should have a series of tests that use small (probably no bigger than
    # 5x5 for your sanity) model grids for which you are able to hand calculate
    # what the answer is.

    # you then have your component calculate the answer, and create a value or
    # array that has the known correct values, and then assert that these two
    # things are equal.
    pass
