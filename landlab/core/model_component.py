#! /usr/bin/env python
from __future__ import print_function

import os
import textwrap


_VAR_HELP_MESSAGE = """
name: {name}
description:
{desc}
units: {units}
at: {loc}
"""


class classproperty(property):
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()


class Component(object):
    _input_var_names = set()
    _output_var_names = set()
    _var_units = dict()

    def __init__(self, grid, map_vars=None):
        map_vars = map_vars or {}
        self._grid = grid

        for (location, vars) in map_vars.items():
            for (dest, src) in vars.items():
                grid.add_field(location, dest,
                               grid.field_values(location, src))

    @classproperty
    @classmethod
    def input_var_names(cls):
        return tuple(cls._input_var_names)

    @classproperty
    @classmethod
    def output_var_names(self):
        return tuple(self._output_var_names)

    @classproperty
    @classmethod
    def name(self):
        return self._name

    @classproperty
    @classmethod
    def units(self):
        return tuple(self._var_units.items())

    @classproperty
    @classmethod
    def var_units(self):
        return tuple(self._var_units.items())

    @classproperty
    @classmethod
    def var_definitions(self):
        return tuple(self._var_doc.items())

    @classmethod
    def var_help(cls, name):
        desc = os.linesep.join(textwrap.wrap(cls._var_doc[name],
                                             initial_indent='  ',
                                             subsequent_indent='  '))
        units = cls._var_units[name]
        loc = cls._var_mapping[name]

        help = _VAR_HELP_MESSAGE.format(name=name, desc=desc, units=units,
                                        loc=loc)

        print(help.strip())

    @classproperty
    @classmethod
    def var_mapping(self):
        """var_mapping
        This is 'node', 'cell', 'active_link', etc.
        """
        return tuple(self._var_mapping.items())

    @property
    def shape(self):
        return self.grid._shape

    @property
    def grid(self):
        return self._grid

    @property
    def coords(self):
        return (self.grid.node_x, self.grid.node_y)

    def imshow(self, name, **kwds):
        self._grid.imshow(name, **kwds)
