import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_equal

from landlab import RasterModelGrid


ELEMENTS = ['node', 'cell', 'link', 'face']
ELEMENTS += ['active_' + name for name in ELEMENTS]
TYPES = ['float', 'int', 'bool']

def generate_zeros_tests():
    for element in ELEMENTS:
        for type in TYPES:
            def _test():
                rmg = RasterModelGrid(4, 5)
                number_of_elements = rmg.number_of_elements(element)
                assert_array_equal(rmg.zeros(centering=element),
                                   np.zeros(number_of_elements, dtype=np.float))
            _test.description = '%s.test_zeros_%s_%s'% (__name__, type, element)
            yield _test


def generate_ones_tests():
    for element in ELEMENTS:
        for type in TYPES:
            def _test():
                rmg = RasterModelGrid(4, 5)
                number_of_elements = rmg.number_of_elements(element)
                assert_array_equal(rmg.ones(centering=element),
                                   np.ones(number_of_elements, dtype=np.float))
            _test.description = '%s.test_zeros_%s_%s'% (__name__, type, element)
            yield _test


def generate_empty_tests():
    elements = ['node', 'cell', 'link', 'face']
    elements += ['active_' + name for name in elements]

    types = ['float', 'int', 'bool']

    for element in ELEMENTS:
        for type in TYPES:
            def _test():
                rmg = RasterModelGrid(4, 5)
                number_of_elements = rmg.number_of_elements(element)
                assert_equal(rmg.empty(centering=element).size,
                             number_of_elements)
            _test.description = '%s.test_zeros_%s_%s'% (__name__, type, element)
            yield _test

