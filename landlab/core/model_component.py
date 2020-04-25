#! /usr/bin/env python
"""Defines the base component class from which Landlab components inherit.

Base component class methods
++++++++++++++++++++++++++++

.. autosummary::

    ~landlab.core.model_component.Component.name
    ~landlab.core.model_component.Component.from_path
    ~landlab.core.model_component.Component.unit_agnostic
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
"""

import os
import textwrap

import numpy as np

from .. import registry
from ..field import FieldError
from .model_parameter_loader import load_params

_VAR_HELP_MESSAGE = """
name: {name}
description:
{desc}
units: {units}
unit agnostic: {unit_agnostic}
at: {loc}
intent: {intent}
"""


class classproperty(property):
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()


class Component:
    """Base component class from which Landlab components inherit."""

    _info = {}
    _name = None
    _cite_as = ""
    _unit_agnostic = None

    def __new__(cls, *args, **kwds):
        registry.add(cls)
        return object.__new__(cls)

    def __init__(self, grid):
        self._grid = grid
        self._current_time = None
        # ensure that required input fields exist
        for name in self._info.keys():
            at = self._info[name]["mapping"]
            optional = self._info[name]["optional"]
            in_true = "in" in self._info[name]["intent"]
            if (in_true) and (not optional):
                # if required input, verify that it exists.
                if name not in self._grid[at]:
                    raise FieldError(
                        "{component} is missing required input field: {name} at {at}".format(
                            component=self._name, name=name, at=at
                        )
                    )

                # if required input exists, check dtype.
                field = self._grid[at][name]
                dtype = self._info[name]["dtype"]

                if field.dtype != dtype:
                    raise FieldError(
                        "{component} required input variable: {name} at {at} has incorrect dtype. dtype must be {dtype} and is {actual}".format(
                            component=self._name,
                            name=name,
                            at=at,
                            dtype=dtype,
                            actual=field.dtype,
                        )
                    )

            # if optional input exists, check dtype
            if in_true and optional:

                if name in self._grid[at]:
                    field = self._grid[at][name]
                    dtype = self._info[name]["dtype"]

                    if field.dtype != dtype:
                        raise FieldError(
                            "{component} optional input variable: {name} at {at} has incorrect dtype. dtype must be {dtype} and is {actual}".format(
                                component=self._name,
                                name=name,
                                at=at,
                                dtype=dtype,
                                actual=field.dtype,
                            )
                        )

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
            with open(path, "r") as fp:
                params = load_params(fp)
        else:
            params = load_params(path)
        return cls(grid, **params)

    @classproperty
    @classmethod
    def cite_as(cls):
        """Citation information for component.

        Return required software citation, if any. An empty string indicates
        that no citations other than the standard Landlab package citations are
        needed for the component.

        Citations are provided in BibTeX format.

        Returns
        -------
        cite_as
        """
        return cls._cite_as

    @property
    def current_time(self):
        """Current time.

        Some components may keep track of the current time. In this case, the
        ``current_time`` attribute is incremented. Otherwise it is set to None.

        Returns
        -------
        current_time
        """
        return self._current_time

    @current_time.setter
    def current_time(self, new_time):
        if self._current_time is not None:
            assert new_time > self._current_time
        self._current_time = new_time

    @classproperty
    @classmethod
    def input_var_names(cls):
        """Names of fields that are used by the component.

        Returns
        -------
        tuple of str
            Tuple of field names.
        """
        input_var_names = [
            name
            for name in cls._info.keys()
            if (not cls._info[name]["optional"]) and ("in" in cls._info[name]["intent"])
        ]
        return tuple(sorted(input_var_names))

    @classproperty
    @classmethod
    def output_var_names(cls):
        """Names of fields that are provided by the component.

        Returns
        -------
        tuple of str
            Tuple of field names.
        """
        output_var_names = [
            name
            for name in cls._info.keys()
            if (not cls._info[name]["optional"])
            and ("out" in cls._info[name]["intent"])
        ]
        return tuple(sorted(output_var_names))

    @classproperty
    @classmethod
    def optional_var_names(cls):
        """Names of fields that are optionally provided by the component, if
        any.

        Returns
        -------
        tuple of str
            Tuple of field names.
        """
        optional_var_names = [
            name for name in cls._info.keys() if cls._info[name]["optional"]
        ]
        return tuple(sorted(optional_var_names))

    @classmethod
    def var_type(cls, name):
        """Returns the dtype of a field (float, int, bool, str...).

        Parameters
        ----------
        name : str
            A field name.

        Returns
        -------
        dtype
            The dtype of the field.
        """
        return cls._info[name]["dtype"]

    @classproperty
    @classmethod
    def name(cls):
        """Name of the component.

        Returns
        -------
        str
            Component name.
        """
        return cls._name

    @classproperty
    @classmethod
    def unit_agnostic(cls):
        """Whether the component is unit agnostic.

        If True, then the component is unit agnostic. Under this condition a
        user must still provide consistent units across all input arguments,
        keyword arguments, and fields. However, when ``unit_agnostic`` is True
        the units specified can be interpreted as dimensions.

        When False, then the component requires inputs in the specified units.

        Returns
        -------
        bool

        """
        return cls._unit_agnostic

    @classproperty
    @classmethod
    def units(cls):
        """Get the units for all field values.

        Returns
        -------
        tuple or str
            Units for each field.
        """
        return tuple(
            sorted([(name, cls._info[name]["units"]) for name in cls._info.keys()])
        )

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
        return cls._info[name]["units"]

    @classproperty
    @classmethod
    def definitions(cls):
        """Get a description of each field.

        Returns
        -------
        tuple of (*name*, *description*)
            A description of each field.
        """
        return tuple(
            sorted([(name, cls._info[name]["doc"]) for name in cls._info.keys()])
        )

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
        return cls._info[name]["doc"]

    @classmethod
    def var_help(cls, name):
        """Print a help message for a particular field.

        Parameters
        ----------
        name : str
            A field name.
        """
        desc = os.linesep.join(
            textwrap.wrap(
                cls._info[name]["doc"], initial_indent="  ", subsequent_indent="  "
            )
        )
        units = cls._info[name]["units"]
        loc = cls._info[name]["mapping"]
        intent = cls._info[name]["intent"]

        help = _VAR_HELP_MESSAGE.format(
            name=name,
            desc=desc,
            units=units,
            loc=loc,
            intent=intent,
            unit_agnostic=cls._unit_agnostic,
        )

        print(help.strip())

    @classproperty
    @classmethod
    def var_mapping(cls):
        """Location where variables are defined.

        Returns
        -------
        tuple of (name, location)
            Tuple of variable name and location ('node', 'link', etc.) pairs.
        """
        return tuple(
            sorted([(name, cls._info[name]["mapping"]) for name in cls._info.keys()])
        )

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
        return cls._info[name]["mapping"]

    def initialize_output_fields(self, values_per_element=None):
        """Create fields for a component based on its input and output var
        names.

        This method will create new fields (without overwrite) for any fields
        output by, but not supplied to, the component. New fields are
        initialized to zero. Ignores optional fields. New fields are created as
        arrays of floats, unless the component specifies the variable type.

        Parameters
        ----------
        values_per_element: int (optional)
            On occasion, it is necessary to create a field that is of size
            (n_grid_elements, values_per_element) instead of the default size
            (n_grid_elements,). Use this keyword argument to acomplish this
            task.
        """
        for name in self._info.keys():
            at = self._info[name]["mapping"]
            optional = self._info[name]["optional"]
            out_true = "out" in self._info[name]["intent"]
            if (out_true) and (not optional) and (name not in self._grid[at]):

                type_in = self.var_type(name)
                num_elements = self._grid.size(at)

                if values_per_element is None:
                    size = num_elements
                else:
                    size = (num_elements, values_per_element)

                init_vals = np.zeros(size, dtype=type_in)
                units_in = self.var_units(name)

                self.grid.add_field(name, init_vals, at=at, units=units_in, copy=False)

    def initialize_optional_output_fields(self):
        """Create fields for a component based on its optional field outputs,
        if declared in _optional_var_names.

        This method will create new fields (without overwrite) for any
        fields output by the component as optional. New fields are
        initialized to zero. New fields are created as arrays of floats,
        unless the component also contains the specifying property
        _var_type.
        """

        for name in self._info.keys():
            at = self._info[name]["mapping"]
            optional = self._info[name]["optional"]
            out_true = "out" in self._info[name]["intent"]
            if (out_true) and (optional) and (name not in self._grid[at]):

                type_in = self.var_type(name)
                init_vals = self.grid.zeros(at, dtype=type_in)
                units_in = self.var_units(name)

                self.grid.add_field(name, init_vals, at=at, units=units_in, copy=False)

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
        """Return the coordinates of nodes on grid attached to the
        component."""
        return (self.grid.node_x, self.grid.node_y)
