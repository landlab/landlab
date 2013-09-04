#! /usr/bin/env python

import types
import inspect

from landlab.field import ScalarDataFields


class ModelDataFields(dict):
    def __init__(self):
        super(ModelDataFields, self).__init__()

    @property
    def groups(self):
        return set(self.keys())

    def new_field_location(self, group, size):
        if group not in self:
            self[group] = ScalarDataFields(size)
            setattr(self, 'at_' + group, self[group])
        else:
            raise ValueError('ModelDataField already contains %s' % group)

    def field_values(self, group, field):
        return self[group][field]

    def field_units(self, group, field):
        return self[group].units[field]


def get_method_from_class(clazz, method_name):
    method = getattr(clazz, method_name)
    if inspect.ismethod(method):
        return method
    else:
        raise TypeError('%s is not a bound method' % method_name)


def is_regular_method_name(name):
    return not (name.startswith('__') and name.endswith('__'))


def regular_method_names(clazz):
    methods = set()
    for (name, _) in inspect.getmembers(clazz, inspect.ismethod):
        if is_regular_method_name(name):
            methods.add(name)
    return methods


def prepend_arg_list_with_dict_value(destination_class, method):
    def func_with_dict_value_first(self, key, *args, **kwds):
        return method(self[key], *args, **kwds)
    func_with_dict_value_first.__doc__ = method.__doc__
    func_with_dict_value_first.__name__ = method.__name__
    return types.MethodType(func_with_dict_value_first, None,
                            destination_class)


def add_methods_to_grouped_fields_class():
    for name in (regular_method_names(ScalarDataFields) -
                 regular_method_names(ModelDataFields)):
        method = get_method_from_class(ScalarDataFields, name)
        new_method = prepend_arg_list_with_dict_value(ModelDataFields,
                                                      method)
        setattr(ModelDataFields, name, new_method)


add_methods_to_grouped_fields_class()
