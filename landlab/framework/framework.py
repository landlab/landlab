#! /usr/bin/env python
from landlab import Arena, Palette

from .component import load_landlab_components


class Error(Exception):
    """
    Base exception for this module
    """


class Framework(object):
    """
    A framework for connecting and running component from The Landlab.
    """

    def __init__(self):
        self._palette = Palette(**load_landlab_components())
        self._arena = Arena()

    def instantiate(self, name):
        """
        Instantiate a component called *name* from the palette and move it to
        the arena.
        """
        try:
            self._arena.instantiate(self._palette.get(name), name)
        except KeyError:
            pass

    def remove(self, name):
        """
        Remove a component called *name* from the arena.
        """
        try:
            self._arena.remove(name)
        except KeyError:
            pass

    def list_palette(self):
        """
        Get a list of names of the components in the palette.
        """
        return self._palette.list()

    def list_arena(self):
        """
        Get a list of names of the components in the arena.
        """
        return self._arena.list()

    def arena_uses(self):
        """
        Get a list of variable names that components in the arena use.
        """
        return self._arena.uses()

    def arena_provides(self):
        """
        Get a list of variable names that components in the arena provide.
        """
        return self._arena.provides()

    def palette_uses(self):
        """
        Get a list of variable names that components in the palette use.
        """
        return self._palette.uses()

    def palette_provides(self):
        """
        Get a list of variable names that components in the palette provide.
        """
        return self._palette.provides()

    def __repr__(self):
        return "Framework(%s)" % ", ".join(self._palette.keys())
