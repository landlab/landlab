#! /usr/bin/env python

import inspect

from landlab.plugin import load_landlab_plugins
from landlab.collections import Palette, Arena


class Error(Exception):
    """
    Base exception for this module
    """

class Framework(object):
    """
    A framework for connecting and running component from The Landlab.
    """
    def __init__(self):
        self._palette = Palette(**load_landlab_plugins())
        self._arena = Arena()

    def instantiate(self, name):
        """
        Instantiate a component from the palette and move it to the arena.

        :name: Name of the component to instantiate
        """
        try:
            self._arena.instantiate(self._palette.get(name), name)
        except KeyError:
            pass

    def remove(self, name):
        """
        Remove a component from the arena.

        :name: Name of the component to remove
        """
        try:
            self._arena.remove(name)
        except KeyError:
            pass

    def list_palette(self):
        """
        Get a list of names of the components in the palette.

        :returns: A list of component names
        """
        return self._palette.list()

    def list_arena(self):
        """
        Get a list of names of the components in the arena.

        :returns: A list of component names
        """
        return self._arena.list()

    def arena_uses(self):
        """
        Get a list of variable names components in the arena use.

        :returns: A list of variable names
        """
        return self._arena.uses()

    def arena_provides(self):
        """
        Get a list of variable names components in the arena provide.

        :returns: A list of variable names
        """
        return self._arena.provides()

    def palette_uses(self):
        """
        Get a list of variable names components in the palette use.

        :returns: A list of variable names
        """
        return self._palette.uses()

    def palette_provides(self):
        """
        Get a list of variable names components in the palette provide.

        :returns: A list of variable names
        """
        return self._palette.provides()

    def __repr__(self):
        return 'Framework(%s)' % ', '.join(self._palette.keys())
