#! /usr/bin/env python
"""
Utility functions for loading components for The Landlab.
"""

import os

from landlab.framework.interfaces import BmiBase

_COMPONENT_PATH = [os.path.join(os.path.dirname(__file__), "..", "components")]

try:
    paths = os.environ["LANDLAB_PATH"].split(os.pathsep)
except KeyError:
    pass
else:
    _COMPONENT_PATH = paths + _COMPONENT_PATH


def iscomponent(value, cls):
    """
    Check if *value* is a component for The Landlab. *value* is a component
    if it implements the *cls* or it is an instance of *cls*.

    Returns ``True`` if *value* is a component, otherwise ``False``.
    """
    try:
        return (
            cls in value.__implements__
            or cls.__name__ in value.__implements__
            or isinstance(value, cls)
        )
    except AttributeError:
        return False


def load_components_from_dir(path, cls):
    """
    Look for components for Landlab in *path*. Identify components as being an
    instance of *cls*. Returns a dictionary of discovered component names as
    keys and component classes as values.
    """
    import sys
    import imp

    components = {}

    sys.path.insert(0, path)
    cwd = os.getcwd()

    os.chdir(path)
    for file_name in os.listdir("."):
        if os.path.isfile(file_name) and file_name.endswith(".py"):
            (mod_name, _) = os.path.splitext(file_name)
            mod = imp.load_module(mod_name, *imp.find_module(mod_name))
            for (name, value) in mod.__dict__.items():
                if iscomponent(value, cls):
                    components[name] = value

    os.chdir(cwd)
    sys.path.pop(0)

    return components


def load_components(cls, paths=None):
    """
    Load components from a series of directories.

    Components found earlier in the search path order override those
    discovered later. Use the *paths* keyword to specify a list of paths to
    search for components.

    .. seealso::

        :func:`load_components_from_dir`
    """
    if not paths:
        paths = _COMPONENT_PATH

    components = {}
    for path in paths[::-1]:
        components.update(load_components_from_dir(path, cls))
    return components


def load_landlab_components(paths=None):
    """
    Load components for The Landlab. These are classes that implement BmiBase.
    See :func:`load_components_from_dir` for the meaning of *paths* keyword.

    .. seealso::

        :func:`load_components_from_dir`
    """
    return load_components(BmiBase, paths=paths)
