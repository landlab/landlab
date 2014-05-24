#! /usr/bin/env python

import unittest

from landlab.field.grouped import (_is_regular_method_name,
                                   _get_method_from_class,
                                   _prepend_arg_list_with_dict_value,
                                   _regular_method_names)


class TestRegularMethodName(unittest.TestCase):
    def test_special_name(self):
        self.assertFalse(_is_regular_method_name('__init__'))

    def test_regular_name(self):
        self.assertTrue(_is_regular_method_name('do_something'))

    def test_private_name(self):
        self.assertTrue(_is_regular_method_name('_do_something_privately'))


class TestRegularMethods(unittest.TestCase):
    def test_no_regulars(self):
        class A(object):
            pass
        self.assertSetEqual(set(), _regular_method_names(A))

    def test_one_regular(self):
        class A(object):
            def do_something(self):
                pass
        self.assertSetEqual(set(['do_something']), _regular_method_names(A))

    def test_multiple_regular(self):
        class A(object):
            def do_something(self):
                pass
            def do_something_else(self):
                pass
            def stop(self):
                pass
        self.assertSetEqual(set(['do_something', 'do_something_else', 'stop']),
                            _regular_method_names(A))

    def test_staticmethod(self):
        class A(object):
            @staticmethod
            def do_something(self):
                pass
        self.assertSetEqual(set(), _regular_method_names(A))

    def test_attribute(self):
        class A(object):
            do_something = 'nothing'
        self.assertSetEqual(set(), _regular_method_names(A))


class TestGetMethodFromClass(unittest.TestCase):
    def test_exists(self):
        class A(object):
            def do_something(self):
                pass
        self.assertEqual(A.do_something,
                         _get_method_from_class(A, 'do_something'))

    def test_absent(self):
        class A(object):
            pass
        with self.assertRaises(AttributeError):
            _get_method_from_class(A, 'do_something')

    def test_is_not_method(self):
        class A(object):
            do_something = False
        with self.assertRaises(TypeError):
            _get_method_from_class(A, 'do_something')

class TestPrependArgList(unittest.TestCase):
    def test_(self):
        class A(dict):
            pass
        def func():
            pass
        setattr(A, 'do_something', _prepend_arg_list_with_dict_value(A, func))
        self.assertTrue(hasattr(A, 'do_something'))
