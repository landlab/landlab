from numpy.testing import assert_array_equal


class NumpyArrayTestingMixIn(object):
    def assertArrayEqual(self, actual, expected):
        try:
            assert_array_equal(actual, expected)
        except AssertionError as error:
            self.fail(error)
