import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid

ELEMENTS = ["node", "cell", "link", "face"]
# ELEMENTS += ['core_node', 'core_cell', 'active_link', 'active_face']
TYPES = ["float", "int", "bool"]


def generate_zeros_tests():
    for element in ELEMENTS:
        for type in TYPES:

            def _test():
                rmg = RasterModelGrid((4, 5))
                number_of_elements = rmg.number_of_elements(element)
                assert_array_equal(
                    rmg.zeros(centering=element),
                    np.zeros(number_of_elements, dtype=np.float),
                )

            _test.description = "%s.test_zeros_%s_%s" % (__name__, type, element)
            yield _test


def generate_add_zeros_tests():
    for element in ELEMENTS:
        for type in TYPES:

            def _test():
                rmg = RasterModelGrid((4, 5))
                number_of_elements = rmg.number_of_elements(element)
                rtn_values = rmg.add_zeros(element, "name")
                assert rtn_values is rmg.field_values(element, "name")
                assert_array_equal(
                    rtn_values, np.zeros(number_of_elements, dtype=np.float)
                )

            _test.description = "%s.test_add_zeros_%s_%s" % (__name__, type, element)
            yield _test


def generate_ones_tests():
    for element in ELEMENTS:
        for type in TYPES:

            def _test():
                rmg = RasterModelGrid((4, 5))
                number_of_elements = rmg.number_of_elements(element)
                assert_array_equal(
                    rmg.ones(centering=element),
                    np.ones(number_of_elements, dtype=np.float),
                )

            _test.description = "%s.test_zeros_%s_%s" % (__name__, type, element)
            yield _test


def generate_add_ones_tests():
    for element in ELEMENTS:
        for type in TYPES:

            def _test():
                rmg = RasterModelGrid((4, 5))
                number_of_elements = rmg.number_of_elements(element)
                rtn_values = rmg.add_ones(element, "name")
                assert rtn_values is rmg.field_values(element, "name")
                assert_array_equal(
                    rtn_values, np.ones(number_of_elements, dtype=np.float)
                )

            _test.description = "%s.test_add_zeros_%s_%s" % (__name__, type, element)
            yield _test


def generate_empty_tests():
    for element in ELEMENTS:
        for type in TYPES:

            def _test():
                rmg = RasterModelGrid((4, 5))
                number_of_elements = rmg.number_of_elements(element)
                assert rmg.empty(centering=element).size == number_of_elements

            _test.description = "%s.test_zeros_%s_%s" % (__name__, type, element)
            yield _test


def generate_add_empty_tests():
    for element in ELEMENTS:
        for type in TYPES:

            def _test():
                rmg = RasterModelGrid((4, 5))
                number_of_elements = rmg.number_of_elements(element)
                rtn_values = rmg.add_empty(element, "name")
                assert rtn_values is rmg.field_values(element, "name")
                assert_array_equal(rtn_values.size, number_of_elements)

            _test.description = "%s.test_zeros_%s_%s" % (__name__, type, element)
            yield _test
