#! /usr/bin/env python
"""
Utility functions for loading plugins.
"""

import os
from cmt.bmi import BmiBase


_PLUGIN_PATH = [
    os.path.join(os.path.dirname(__file__), '..', 'plugins'),
]

try:
    paths = os.environ['LANDLAB_PATH'].split(os.pathsep)
except KeyError:
    pass
else:
    _PLUGIN_PATH = paths + _PLUGIN_PATH


def isplugin(value, cls):
    """
    Check if *value* is a plugin component for The Landlab. *value* is a plugin
    if it implements the *cls* or it is an instance of *cls*.

    Returns ``True`` if *value* is a plugin, otherwise ``False``.
    """
    try:
        return (cls in value.__implements__ or
                cls.__name__ in value.__implements__ or
                isinstance(value, cls))
    except AttributeError:
        return False


def load_plugins_from_dir(path, cls):
    """
    Look for plugins for The Landlab in *path*. Identify plugins as being an
    instance of *cls*. Returns a dictionary of discovered plugin names as
    keys and plugin classes as values.
    """
    import sys
    import imp

    plugins = {}

    sys.path.insert(0, path)
    cwd = os.getcwd()

    os.chdir(path)
    for name in os.listdir('.'):
        if os.path.isfile(name) and name.endswith('.py'):
            (mod_name, _) = os.path.splitext(name)
            mod = imp.load_module(mod_name, *imp.find_module(mod_name))
            for (name, value) in mod.__dict__.items():
                if isplugin(value, cls):
                    plugins[name] = value

    os.chdir(cwd)
    sys.path.pop(0)

    return plugins


def load_plugins(cls, paths=_PLUGIN_PATH):
    """
    Load plugins from a series of directories. Plugins found earlier in the
    search path order override those discovered later. Use the *paths* keyword
    to specify a list of paths to search for plugins.

    .. seealso::
        load_plugins_from_dir
    """
    plugins = {}
    for path in paths[::-1]:
        plugins.update(load_plugins_from_dir(path, cls))
    return plugins


def load_landlab_plugins(paths=_PLUGIN_PATH):
    """
    Load plugins for The Landlab. These are classes that implement BmiBase. Use
    the *paths* keyword to specify a list of paths to search for plugins.

    .. seealso::
        load_plugins_from_dir
    """
    plugins = {}
    for path in paths[::-1]:
        plugins.update(load_plugins_from_dir(path, cls))
    return plugins

def load_landlab_plugins(paths=_PLUGIN_PATH):
    """
    Load plugins for The Landlab. These are classes that implement BmiBase.

    .. seealso::
        load_plugins_from_dir
    """
    return load_plugins(BmiBase, paths=paths)

    .. seealso::
        load_plugins_from_dir
    """
    return load_plugins(BmiBase, paths=paths)
