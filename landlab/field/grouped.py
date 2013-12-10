#! /usr/bin/env python
"""
The ModelDataFields class holds at set of ScalarDataFields that are separated
into groups. A typical use for this class would be to define the groups as
being locations on a grid where the values are defined. For instance, the
groups could be *node*, *cell*, *link*, and *face*.

Most of the method functions for ModelDataFields are the same as those for
the ScalarDataFields class but with the first argument being a string that
defines the group name.

Create a sets of data fields defined at *node* and *cell*. Each set can
have a differenct number of values.

>>> fields = ModelDataFields()
>>> fields.new_field_location('node', 12)
>>> fields.new_field_location('cell', 2)

Create some new value arrays for each of the data fields.

>>> fields.ones('node')
array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
>>> fields.zeros('cell')
array([ 0.,  0.])

Create new value arrays and add them to the data fields. Because the data
fields are in different groups (node and cell), they can have the same
name.

>>> fields.add_ones('node', 'planet_surface__elevation')
>>> fields.at_node['planet_surface__elevation']
array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

>>> fields.add_ones('cell', 'planet_surface__elevation')
>>> fields.at_cell['planet_surface__elevation']
array([ 1.,  1.])
"""

import types
import inspect

from landlab.field import ScalarDataFields


class ModelDataFields(object):
    def __init__(self, **kwds):
        #print 'ModelDataFields.__init__'
        self._groups = dict()
        super(ModelDataFields, self).__init__(**kwds)

    @property
    def groups(self):
        return set(self._groups.keys())

    def has_group(self, group):
        return group in self._groups

    def new_field_location(self, group, size):
        if self.has_group(group):
            raise ValueError('ModelDataFields already contains %s' % group)
        else:
            self._groups[group] = ScalarDataFields(size)
            setattr(self, 'at_' + group, self[group])

    def field_values(self, group, field):
        return self[group][field]

    def field_units(self, group, field):
        return self[group].units[field]

    def __getitem__(self, group):
        return self._groups[group]


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
    def func_with_dict_value_first(self, group, *args, **kwds):
        return method(self[group], *args, **kwds)
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
