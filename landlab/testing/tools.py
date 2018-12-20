import os
import shutil
import tempfile
from distutils.dir_util import mkpath

import numpy as np


class cd(object):
    """Context that changes to a new directory.

    Examples
    --------
    >>> import os, tempfile
    >>> from landlab.testing.tools import cd

    Create a temporary directory for testing.

    >>> test_dir = os.path.realpath(tempfile.mkdtemp())

    Withing the context, we're in the new working directory, after exiting
    the context we're back where we started.

    >>> this_dir = os.getcwd()
    >>> with cd(test_dir) as _:
    ...     wdir = os.getcwd()
    >>> test_dir == wdir
    True
    >>> os.getcwd() == this_dir
    True

    If the new working directory does not exists, create it.

    >>> new_dir = os.path.join(test_dir, 'testing.d')
    >>> os.path.exists(new_dir)
    False
    >>> with cd(new_dir) as _:
    ...     wdir = os.getcwd()
    >>> os.path.exists(new_dir)
    True
    >>> wdir == new_dir
    True
    >>> os.getcwd() == this_dir
    True
    """

    def __init__(self, path_to_dir):
        self._dir = path_to_dir

    def __enter__(self):
        self._starting_dir = os.path.abspath(os.getcwd())
        if not os.path.isdir(self._dir):
            mkpath(self._dir)
        os.chdir(self._dir)
        return os.path.abspath(os.getcwd())

    def __exit__(self, ex_type, ex_value, traceback):
        os.chdir(self._starting_dir)


class cdtemp(object):
    """Context that creates and changes to a temporary directory.

    Examples
    --------
    >>> import os
    >>> from landlab.testing.tools import cdtemp

    Change to the newly-created temporary directory after entering the
    context. Upon exiting, remove the temporary directory and return to the
    original working directory.

    >>> this_dir = os.getcwd()
    >>> with cdtemp() as tdir:
    ...     wdir = os.getcwd()
    >>> this_dir == os.getcwd()
    True
    >>> os.path.exists(wdir)
    False
    """

    def __init__(self, **kwds):
        self._kwds = kwds
        self._tmp_dir = None

    def __enter__(self):
        self._starting_dir = os.path.abspath(os.getcwd())
        self._tmp_dir = tempfile.mkdtemp(**self._kwds)
        os.chdir(self._tmp_dir)
        return os.path.abspath(self._tmp_dir)

    def __exit__(self, ex_type, ex_value, traceback):
        os.chdir(self._starting_dir)
        shutil.rmtree(self._tmp_dir)


def assert_array_is_int(x):
    assert x.dtype == np.int32 or x.dtype == np.int64


def assert_close(val1, val2, msg=None):
    assert np.allclose(val1, val2), msg


def assert_is(expr1, expr2, msg=None):
    assert expr1 is expr2, msg


def assert_is_not(expr1, expr2, msg=None):
    assert expr1 is not expr2, msg


def assert_is_none(expr, msg=None):
    assert expr is None, msg


def assert_is_not_none(expr, msg=None):
    assert expr is not None, msg


def assert_in(first, second, msg=None):
    assert first in second, msg


def assert_is_not_in(first, second, msg=None):
    assert first not in second, msg


def assert_is_instance(obj, cls, msg=None):
    assert isinstance(obj, cls), msg


def assert_not_is_instance(obj, cls, msg=None):
    assert not isinstance(obj, cls), msg


def assert_list_equal(list1, list2, msg=None):
    assert isinstance(list1, list), msg
    assert isinstance(list2, list), msg
    for a, b in zip(list1, list2):
        assert a == b, msg


def assert_tuple_equal(tuple1, tuple2, msg=None):
    assert isinstance(tuple1, tuple), msg
    assert isinstance(tuple1, tuple), msg
    for a, b in zip(tuple1, tuple2):
        assert a == b, msg


def assert_dict_equal(dict1, dict2, msg=None):
    assert isinstance(dict1, dict), msg
    assert isinstance(dict2, dict), msg
    assert dict1 == dict2, msg


def assert_set_equal(set1, set2, msg=None):
    assert isinstance(set1, set), msg
    assert isinstance(set2, set), msg
    assert set1 == set2, msg
