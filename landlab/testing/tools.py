from nose.tools import assert_true


def assert_is(expr1, expr2, msg=None):
    assert_true(expr1 is expr2, msg=msg)


def assert_is_not(expr1, expr2, msg=None):
    assert_true(expr1 is not expr2, msg=msg)


def assert_is_none(expr, msg=None):
    assert_true(expr is None, msg=msg)


def assert_is_not_none(expr, msg=None):
    assert_true(expr is not None, msg=msg)


def assert_in(first, second, msg=None):
    assert_true(first in second, msg=msg)


def assert_is_not_in(first, second, msg=None):
    assert_true(first not in second, msg=msg)


def assert_is_not(expr1, expr2, msg=None):
    assert_true(expr1 is not expr2, msg=msg)


def assert_is_instance(obj, cls, msg=None):
    assert_true(isinstance(obj, cls), msg=msg)


def assert_not_is_instance(obj, cls, msg=None):
    assert_true(not isinstance(obj, cls), msg=msg)
