#! /usr/bin/env python
"""
Defines the base component class from which Landlab components inherit.

Base component class methods
++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.core.model_component.Component.from_path
    ~landlab.core.model_component.Component.name
    ~landlab.core.model_component.Component.units
    ~landlab.core.model_component.Component.definitions
    ~landlab.core.model_component.Component.input_var_names
    ~landlab.core.model_component.Component.output_var_names
    ~landlab.core.model_component.Component.optional_var_names
    ~landlab.core.model_component.Component.var_type
    ~landlab.core.model_component.Component.var_units
    ~landlab.core.model_component.Component.var_definition
    ~landlab.core.model_component.Component.var_mapping
    ~landlab.core.model_component.Component.var_loc
    ~landlab.core.model_component.Component.var_help
    ~landlab.core.model_component.Component.initialize_output_fields
    ~landlab.core.model_component.Component.initialize_optional_output_fields
    ~landlab.core.model_component.Component.shape
    ~landlab.core.model_component.Component.grid
    ~landlab.core.model_component.Component.coords
    ~landlab.core.model_component.Component.imshow
"""

from __future__ import print_function

import os
import textwrap
import warnings
import inspect


_VAR_HELP_MESSAGE = """
name: {name}
description:
{desc}
units: {units}
at: {loc}
intent: {intent}
"""


class classproperty(property):
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()


class Component(object):
    """
    Defines the base component class from which Landlab components inherit.

    Base component class methods
    ++++++++++++++++++++++++++++

    .. autosummary::
        :toctree: generated/

        ~landlab.core.model_component.Component.from_path
        ~landlab.core.model_component.Component.name
        ~landlab.core.model_component.Component.units
        ~landlab.core.model_component.Component.definitions
        ~landlab.core.model_component.Component.input_var_names
        ~landlab.core.model_component.Component.output_var_names
        ~landlab.core.model_component.Component.optional_var_names
        ~landlab.core.model_component.Component.var_type
        ~landlab.core.model_component.Component.var_units
        ~landlab.core.model_component.Component.var_definition
        ~landlab.core.model_component.Component.var_mapping
        ~landlab.core.model_component.Component.var_loc
        ~landlab.core.model_component.Component.var_help
        ~landlab.core.model_component.Component.initialize_output_fields
        ~landlab.core.model_component.Component.initialize_optional_output_fields
        ~landlab.core.model_component.Component.shape
        ~landlab.core.model_component.Component.grid
        ~landlab.core.model_component.Component.coords
        ~landlab.core.model_component.Component.imshow
    """
    _input_var_names = set()
    _output_var_names = set()
    _optional_var_names = set()
    _var_units = dict()

    def __init__(self, grid, map_vars=None, **kwds):
        map_vars = map_vars or {}
        self._grid = grid

        for (location, vars) in map_vars.items():
            for (dest, src) in vars.items():
                grid.add_field(location, dest,
                               grid.field_values(location, src))

        for key in kwds:
            component_name = inspect.getmro(self.__class__)[0].__name__
            warnings.warn(
                "Ingnoring unrecognized input parameter, '{param}', for "
                "{name} component".format(name=component_name, param=key))

    @classmethod
    def from_path(cls, grid, path):
        """Create a component from an input file.

        Parameters
        ----------
        grid : ModelGrid
            A landlab grid.
        path : str or file_like
            Path to a parameter file, contents of a parameter file, or
            a file-like object.

        Returns
        -------
        Component
            A newly-created component.
        """
        if os.path.isfile(path):
            with open(path, 'r') as fp:
                params = load_params(fp)
        else:
            params = load_params(path)
        return cls(grid, **params)

    @classproperty
    @classmethod
    def input_var_names(cls):
        """Names of fields that are used by the component.

        Returns
        -------
        tuple of str
            Tuple of field names.
        """
        return tuple(cls._input_var_names)

    @classproperty
    @classmethod
    def output_var_names(self):
        """Names of fields that are provided by the component.

        Returns
        -------
        tuple of str
            Tuple of field names.
        """
        return tuple(self._output_var_names)

    @classproperty
    @classmethod
    def optional_var_names(self):
        """
        Names of fields that are optionally provided by the component, if
        any.

        Returns
        -------
        tuple of str
            Tuple of field names.
        """
        try:
            return tuple(self._optional_var_names)
        except AttributeError:
            return ()

    @classmethod
    def var_type(cls, name):
        """
        Returns the dtype of a field (float, int, bool, str...), if declared.
        Default is float.

        Parameters
        ----------
        name : str
            A field name.

        Returns
        -------
        dtype
            The dtype of the field.
        """
        try:
            return cls._var_type[name]
        except AttributeError:
            return float

    @classproperty
    @classmethod
    def name(self):
        """Name of the component.

        Returns
        -------
        str
            Component name.
        """
        return self._name

    @classproperty
    @classmethod
    def units(self):
        """Get the units for all field values.

        Returns
        -------
        tuple or str
            Units for each field.
        """
        return tuple(self._var_units.items())

    @classmethod
    def var_units(cls, name):
        """Get the units of a particular field.

        Parameters
        ----------
        name : str
            A field name.

        Returns
        -------
        str
            Units for the given field.
        """
        return cls._var_units[name]

    @classproperty
    @classmethod
    def definitions(cls):
        """Get a description of each field.

        Returns
        -------
        tuple of (*name*, *description*)
            A description of each field.
        """
        return tuple(cls._var_doc.items())

    @classmethod
    def var_definition(cls, name):
        """Get a description of a particular field.

        Parameters
        ----------
        name : str
            A field name.

        Returns
        -------
        tuple of (*name*, *description*)
            A description of each field.
        """
        return cls._var_doc[name]

    @classmethod
    def var_help(cls, name):
        """Print a help message for a particular field.

        Parameters
        ----------
        name : str
            A field name.
        """
        desc = os.linesep.join(textwrap.wrap(cls._var_doc[name],
                                             initial_indent='  ',
                                             subsequent_indent='  '))
        units = cls._var_units[name]
        loc = cls._var_mapping[name]

        intent = ''
        if name in cls._input_var_names:
            intent = 'in'
        if name in cls._output_var_names:
            intent += 'out'

        help = _VAR_HELP_MESSAGE.format(name=name, desc=desc, units=units,
                                        loc=loc, intent=intent)

        print(help.strip())

    @classproperty
    @classmethod
    def var_mapping(self):
        """Location where variables are defined.

        Returns
        -------
        tuple of (name, location)
            Tuple of variable name and location ('node', 'link', etc.) pairs.
        """
        return tuple(self._var_mapping.items())

    @classmethod
    def var_loc(cls, name):
        """Location where a particular variable is defined.

        Parameters
        ----------
        name : str
            A field name.

        Returns
        -------
        str
            The location ('node', 'link', etc.) where a variable is defined.
        """
        return cls._var_mapping[name]

    def initialize_output_fields(self):
        """
        Create fields for a component based on its input and output var names.

        This method will create new fields (without overwrite) for any fields
        output by, but not supplied to, the component. New fields are
        initialized to zero. Ignores optional fields, if specified by
        _optional_var_names. New fields are created as arrays of floats, unless
        the component also contains the specifying property _var_type.
        """
        for field_to_set in (set(self.output_var_names) -
                             set(self.input_var_names) -
                             set(self.optional_var_names)):
            grp = self.var_loc(field_to_set)
            type_in = self.var_type(field_to_set)
            init_vals = self.grid.zeros(grp, dtype=type_in)
            units_in = self.var_units(field_to_set)
            self.grid.add_field(grp,
                                field_to_set,
                                init_vals,
                                units=units_in,
                                copy=False,
                                noclobber=True)

    def initialize_optional_output_fields(self):
        """
        Create fields for a component based on its optional field outputs,
        if declared in _optional_var_names.

        This method will create new fields (without overwrite) for any fields
        output by the component as optional. New fields are
        initialized to zero. New fields are created as arrays of floats, unless
        the component also contains the specifying property _var_type.
        """
        for field_to_set in (set(self.optional_var_names) -
                             set(self.input_var_names)):
            self.grid.add_field(self.var_loc(field_to_set),
                                field_to_set,
                                self.grid.zeros(
                                    dtype=self.var_type(field_to_set)),
                                units=self.var_units(field_to_set),
                                noclobber=True)

    @property
    def shape(self):
        """Return the grid shape attached to the component, if defined."""
        return self.grid._shape

    @property
    def grid(self):
        """Return the grid attached to the component."""
        return self._grid

    @property
    def coords(self):
        """Return the coordinates of nodes on grid attached to the component.
        """
        return (self.grid.node_x, self.grid.node_y)

    def imshow(self, name, **kwds):
        """Plot data on the grid attached to the component.
        """
        self._grid.imshow(name, **kwds)
