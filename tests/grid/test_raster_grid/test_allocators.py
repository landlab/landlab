import numpy as np
from pytest import approx

from landlab import RasterModelGrid

field_dtype = float


def test_zeros(graph_element):
    grid = RasterModelGrid((4, 5))
    number_of_elements = grid.number_of_elements(graph_element)
    assert np.all(
        grid.zeros(at=graph_element, dtype=field_dtype)
        == approx(np.zeros(number_of_elements, dtype=field_dtype))
    )


def test_add_zeros(graph_element):
    grid = RasterModelGrid((4, 5))
    number_of_elements = grid.number_of_elements(graph_element)
    rtn_values = grid.add_zeros("name", at=graph_element, dtype=field_dtype)
    assert rtn_values is grid.field_values(graph_element, "name")
    assert np.all(rtn_values == approx(np.zeros(number_of_elements, dtype=field_dtype)))


def test_ones(graph_element):
    grid = RasterModelGrid((4, 5))
    number_of_elements = grid.number_of_elements(graph_element)
    assert np.all(
        grid.ones(at=graph_element)
        == approx(np.ones(number_of_elements, dtype=field_dtype))
    )


def test_add_ones(graph_element):
    grid = RasterModelGrid((4, 5))
    number_of_elements = grid.number_of_elements(graph_element)
    rtn_values = grid.add_ones("name", at=graph_element, dtype=field_dtype)
    assert rtn_values is grid.field_values(graph_element, "name")
    assert np.all(rtn_values == approx(np.ones(number_of_elements, dtype=field_dtype)))


def test_empty(graph_element):
    grid = RasterModelGrid((4, 5))
    number_of_elements = grid.number_of_elements(graph_element)
    assert grid.empty(at=graph_element, dtype=field_dtype).size == number_of_elements


def test_add_empty(graph_element):
    grid = RasterModelGrid((4, 5))
    rtn_values = grid.add_empty("name", at=graph_element, dtype=field_dtype)
    assert rtn_values is grid.field_values(graph_element, "name")
