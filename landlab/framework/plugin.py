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
    Check if a class is a plugin. A class is a plugin if it implements the
    plugin class or it is an instance of the plugin class.

    :value: The class to check
    :cls: The plugin class
    """
    try:
        return (cls in value.__implements__ or
                cls.__name__ in value.__implements__ or
                isinstance(value, cls))
    except AttributeError:
        return False


def load_plugins_from_dir(path, cls):
    """
    Look for plugins in a given directory. Identify plugins as being an
    instance of the given class, cls.

    :path: Path to directory containing plugins
    :cls: Class that defines the plugin

    :returns: A list of discovered plugins
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
    search path order override those discovered later.

    :cls: Class othat defines the plugin
    :keyword paths: A list of search paths for plugins

    :returns: A list of discovered plugins
    """
    plugins = {}
    for path in paths[::-1]:
        plugins.update(load_plugins_from_dir(path, cls))
    return plugins

def load_landlab_plugins(paths=_PLUGIN_PATH):
    """
    Load plugins for The Landlab. These are classes that implement BmiBase.

    :keyword paths: A list of search paths for plugins

    :returns: A list of discovered plugins
    """
    return load_plugins(BmiBase, paths=paths)
