#! /usr/bin/env python
import types
import inspect

from landlab.field import ScalarDataFields


class ModelDataFields(object):
    """
    The ModelDataFields class holds a set of ScalarDataFields that are separated
    into *groups*. A typical use for this class would be to define the groups as
    being locations on a grid where the values are defined. For instance, the
    groups could be *node*, *cell*, *link*, and *face*.

    Most of the method functions for ModelDataFields are the same as those for
    the ScalarDataFields class but with the first argument being a string that
    defines the group name.

    Attributes
    ----------
    groups

    Examples
    --------

    Create two groups of data fields defined at *node* and *cell*. Each set can
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
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    >>> fields.at_node['planet_surface__elevation']
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    >>> fields.add_ones('cell', 'planet_surface__elevation')
    array([ 1.,  1.])
    >>> fields.at_cell['planet_surface__elevation']
    array([ 1.,  1.])
    """
    def __init__(self, **kwds):
        self._groups = dict()
        super(ModelDataFields, self).__init__(**kwds)

    @property
    def groups(self):
        """Names of all the groups held in the field.

        Returns
        -------
        list
            List of quantity names.
        """
        return set(self._groups.keys())

    def has_group(self, group):
        """Check if the field has *group*.

        Parameters
        ----------
        group: str
            Name of the group.

        Returns
        -------
        boolean
            True if the field contains *group*, otherwise False.
        """
        return group in self._groups

    def new_field_location(self, group, size):
        """Add a new quantity to a field.

        Parameters
        ----------
        group: str
            Name of the new group to add to the field.
        size: int
            Number of elements in the new quantity.

        Raises
        ------
        ValueError
            If the field already contains the group.
        """
        if self.has_group(group):
            raise ValueError('ModelDataFields already contains %s' % group)
        else:
            self._groups[group] = ScalarDataFields(size)
            setattr(self, 'at_' + group, self[group])

    def field_values(self, group, field):
        """Get values of a field.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str
            Name of the field withing *group*.

        Returns
        -------
        array
            The values of the field.

        Raises
        ------
        KeyError
            If either *field* or *group* does not exist.
        """
        return self[group][field]

    def field_units(self, group, field):
        """Get units for a field.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str
            Name of the field withing *group*.

        Returns
        -------
        str
            The units of the field.

        Raises
        ------
        KeyError
            If either *field* or *group* does not exist.

        """
        return self[group].units[field]

    def __getitem__(self, group):
        return self._groups[group]


def _get_method_from_class(clazz, method_name):
    method = getattr(clazz, method_name)
    if inspect.ismethod(method):
        return method
    else:
        raise TypeError('%s is not a bound method' % method_name)


def _is_regular_method_name(name):
    return not (name.startswith('__') and name.endswith('__'))


def _regular_method_names(clazz):
    methods = set()
    for (name, _) in inspect.getmembers(clazz, inspect.ismethod):
        if _is_regular_method_name(name):
            methods.add(name)
    return methods


def _prepend_arg_list_with_dict_value(destination_class, method):
    def func_with_dict_value_first(self, group, *args, **kwds):
        return method(self[group], *args, **kwds)
    func_with_dict_value_first.__doc__ = method.__doc__
    func_with_dict_value_first.__name__ = method.__name__
    return types.MethodType(func_with_dict_value_first, None,
                            destination_class)


def _add_methods_to_grouped_fields_class():
    for name in (_regular_method_names(ScalarDataFields) -
                 _regular_method_names(ModelDataFields)):
        method = _get_method_from_class(ScalarDataFields, name)
        new_method = _prepend_arg_list_with_dict_value(ModelDataFields,
                                                      method)
        setattr(ModelDataFields, name, new_method)


_add_methods_to_grouped_fields_class()
