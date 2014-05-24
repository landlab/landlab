#! /usr/bin/env python
"""Store collections of data fields.
"""
import types
import inspect

from .scalar_data_fields import ScalarDataFields, FieldError

class Error(Exception):
    """Base class for errors in this module."""
    pass


class GroupError(Error, KeyError):
    """Raise this error for a missing group name."""
    def __init__(self, group):
        self._group = group

    def __str__(self):
        return self._group


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

    See Also
    --------
    landlab.field.ScalarDataFields : Data fields within a *group* are
        stored as :class:`landlab.field.ScalarDataFields`.
    landlab.grid.ModelGrid : Inherits from ModelDataFields.

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
        """List of group names.

        Returns
        -------
        set
            Set of quantity names.
        """
        return set(self._groups.keys())

    def has_group(self, group):
        """Check if a group exists.

        Parameters
        ----------
        group: str
            Name of the group.

        Returns
        -------
        boolean
            True if the field contains *group*, otherwise False.

        Examples
        --------
        Check if the field has the groups named *node* or *cell*.

        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 12)
        >>> fields.has_group('node')
        True
        >>> fields.has_group('cell')
        False
        """
        return group in self._groups

    def has_field(self, group, field):
        """Check if a field is in a group.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str
            Name of the field.

        Returns
        -------
        boolean
            ``True`` if the group contains the field, otherwise ``False``.

        Examples
        --------
        Check if the field named ``planet_surface__elevation`` is contained
        in a group.

        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 12)
        >>> _ = fields.add_ones('node', 'planet_surface__elevation')
        >>> fields.has_field('node', 'planet_surface__elevation')
        True
        >>> fields.has_field('cell', 'planet_surface__elevation')
        False
        """
        return group in self._groups

    def keys(self, group):
        """List of field names in a group.

        Returns a list of the field names as a list of strings.

        Parameters
        ----------
        group : str
            Group name.

        Returns
        -------
        list
            List of field names.

        Examples
        --------
        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 4)
        >>> fields.keys('node')
        []
        >>> _ = fields.add_empty('node', 'planet_surface__elevation')
        >>> fields.keys('node')
        ['planet_surface__elevation']
        """
        return self[group].keys()

    def size(self, group):
        """Size of the arrays stored in a group.

        Parameters
        ----------
        group : str
            Group name.

        Returns
        -------
        int
            Array size.

        Examples
        --------
        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 4)
        >>> fields.size('node')
        4
        """
        return self[group].size

    def new_field_location(self, group, size):
        """Add a new quantity to a field.

        Create an empty group into which new fields can be added. The new group
        is created but no memory allocated yet. The dictionary of the new group
        can be through a new *at_* attribute of the class instance.

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

        Examples
        --------
        Create a collection of fields and add two groups, *node* and *cell*,
        to it.

        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 12)
        >>> fields.new_field_location('cell', 2)

        The group names in the collection are retrieved with the *groups*
        attribute as a `set`.

        >>> names = list(fields.groups)
        >>> names.sort()
        >>> names
        ['cell', 'node']

        Access the new (empty) groups with the *at_* attributes.

        >>> fields.at_cell, fields.at_node
        ({}, {})
        """
        if self.has_group(group):
            raise ValueError('ModelDataFields already contains %s' % group)
        else:
            self._groups[group] = ScalarDataFields(size)
            setattr(self, 'at_' + group, self[group])

    def field_values(self, group, field):
        """Get values of a field.

        Given a *group* and a *field*, return a reference to the associated
        data array.

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
        GroupError
            If *group* does not exits
        FieldError
            If *field* does not exits

        Examples
        --------
        Create a group of fields called *node*.

        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 4)

        Add a field, initialized to ones, called *planet_surface__elevation*
        to the *node* group. The *field_values* method returns a reference
        to the field's data.

        >>> _ = fields.add_ones('node', 'planet_surface__elevation')
        >>> fields.field_values('node', 'planet_surface__elevation')
        array([ 1.,  1.,  1.,  1.])

        Raise FieldError if *field* does not exist in *group*.

        >>> fields.field_values('node', 'planet_surface__temperature')
        Traceback (most recent call last):
            ...
        FieldError: planet_surface__temperature

        If *group* does not exists, Raise GroupError.

        >>> fields.field_values('cell', 'planet_surface__elevation')
        Traceback (most recent call last):
            ...
        GroupError: cell
        """
        return self[group][field]

    def field_units(self, group, field):
        """Get units for a field.

        Returns the unit string associated with the data array in *group* and
        *field*.

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
        try:
            return self._groups[group]
        except KeyError:
            raise GroupError(group)


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
