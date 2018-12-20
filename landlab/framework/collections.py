#! /usr/bin/env python
"""
Collections of components
"""

import inspect

from .decorators import is_implementation

# from landlab.interfaces import BmiBase
from .interfaces import BmiBase


class Error(Exception):
    """
    Base exception for this module
    """


class UnknownComponentError(Error):
    """
    Raise this exception when requesting a component unknown to the framework
    """

    def __init__(self, name):
        self._name = name

    def __str__(self):
        return self._name


class MultipleProvidersError(Error):
    """
    Raise this exception if multiple components provide the same variable
    """

    def __init__(self, var_name):
        self._name = var_name

    def __str__(self):
        return "Multiple component provide variable %s" % self._name


class NoProvidersError(Error):
    """
    Raise this exception if no components provide the a variable
    """

    def __init__(self, var_name):
        self._name = var_name

    def __str__(self):
        return "No components provide variable %s" % self._name


class BadVarNameError(Error):
    """
    Raise this exception if a variable name is not found for a component.
    """

    def __init__(self, var_name):
        self._name = var_name

    def __str__(self):
        return "Component does not contain variable %s" % self._name


def get_var_names(component, intent="input"):
    """
    Get a list of input or output variable names from *component* (a BMI-like
    object). Use the *intent* keyword to specify whether to return input or
    output variable names. *intent* must be one of *input* or *output*.
    """
    assert intent in ["input", "output"]

    func = getattr(component, "get_" + intent + "_var_names")
    try:
        var_names = func()
    except TypeError:
        var_names = getattr(component, "_" + intent + "_var_names")

    return var_names


class Collection(dict):
    """
    A collection of components that implement a BmiBase
    """

    def __init__(self, **kwds):
        super(Collection, self).__init__(**kwds)
        for (_, component) in self.items():
            assert is_implementation(component, BmiBase)

    def list(self):
        """
        Get a list of component names in the collection.
        """
        return self.keys()

    def uses(self):
        """
        Get a list of variable names that components in the collection uses.
        """
        return self._var_names(intent="input")

    def provides(self):
        """
        Get a list of variable names that components in the collection
        provides.
        """
        return self._var_names(intent="output")

    def find_provider(self, var_name):
        """
        Find a component in the collection that provides the CSDMS standard name,
        *var_name*. The returned list contains the names of the components
        providing *var_name*.
        """
        return self._find_var_name(var_name, intent="output")

    def find_user(self, var_name):
        """
        Find components in the collection that use the specified CSDMS standard
        variable name, *var_name*. The returned list contains the names of the
        components using *var_name*.
        """
        return self._find_var_name(var_name, intent="input")

    def find_connections(self):
        """
        Find connections between components that use and provide variables.
        The returned dictionary will have keys that are names of components
        that provide varaibles. Its keys will be a dictionary of variable
        names and list of components that use the variable.
        """
        connections = dict()
        for var_name in self.uses():
            providers = self.find_provider(var_name)
            if len(providers) > 1:
                raise MultipleProvidersError(var_name)
            elif len(providers) == 0:
                raise NoProvidersError(var_name)
            else:
                try:
                    connections[providers[0]].update(
                        {var_name: self.find_user(var_name)}
                    )
                except KeyError:
                    connections[providers[0]] = {var_name: self.find_user(var_name)}
        return connections

    def _find_var_name(self, needle, intent="input"):
        """
        Find a variable either used or provided by components in the
        collection.

        :needle: The variable name to search for
        :keyword intent: One of 'input', or 'output'

        :returns: List of component names
        """
        assert intent in ["input", "output"]

        names = set()
        for (name, component) in self.items():
            try:
                for var_name in get_var_names(component, intent=intent):
                    if var_name == needle:
                        names.add(name)
            except AttributeError:
                pass
        return list(names)

    def _var_names(self, intent="input"):
        """
        Get variable names either used or provided by components in the
        collection.

        :keyword intent: One of 'input', or 'output'

        :returns: List of variable names
        """
        assert intent in ["input", "output"]

        names = set()
        for (name, component) in self.items():
            try:
                for var_name in get_var_names(component, intent=intent):
                    names.add(var_name)
            except AttributeError:
                pass
        return list(names)


class Palette(Collection):
    """
    A collection of component classes that have yet to be instantiated.
    """

    def __init__(self, *args, **kwds):
        super(Palette, self).__init__(*args, **kwds)


class Arena(Collection):
    """
    A collection of component instances.
    """

    def __init__(self):
        super(Arena, self).__init__()
        self._connections = dict()

    def instantiate(self, cls, name):
        """
        Instantiate a component and call it *name*. The component, *cls*, is
        either an instance of a BMI-like object or a class that implements the
        BMI. If *cls* is a class, it is instantiated by calling it's __init__
        method without any arguments.
        """
        if inspect.isclass(cls):
            self[name] = cls()
        else:
            self[name] = cls

        try:
            assert is_implementation(type(self[name]), BmiBase)
        except AssertionError:
            self.pop(name)
            raise TypeError("Class is not an implementation of cmt.bmi.BmiBase")
        else:
            self._connections[name] = dict()

    def connect(self, user_name, provider_name, var_name):
        """
        Connect two components through a variable. *user_name* is the name of the
        component that uses *var_name* and *provider_name* is the name of the
        component that provides *var_name*.

        If the arena doesn't contain either *user_name* or *provider_name*,
        :func:`UnknownComponentError` is raised.
        """
        try:
            user = self[user_name]
        except KeyError:
            raise UnknownComponentError(user_name)
        try:
            provider = self[provider_name]
        except KeyError:
            raise UnknownComponentError(provider_name)

        try:
            user.set_value(var_name, provider.get_value(var_name))
        except BadVarNameError:
            raise
        except AttributeError:
            raise

        self._connections[user_name].update({var_name: provider_name})

    def _component_children(self, name):
        """
        Get a set of component names that represent the components that
        provide data to a component.

        :name: Name of the base component

        :return: A set of component names
        """
        children = set()
        for child in self._connections[name].values():
            children.add(child)
        return children

    def walk(self, root, tree=None):
        """Walk a connected set of components.

        Walk a connected set of components with the component named *root*. If
        the *tree* keyword is given, treat it as a list of components already
        in the tree and add to that list. If a component is already in the tree,
        do not iterate through that component. This function returns a list of
        component names in the order they are visited by the walk.
        """
        if tree is None:
            tree = []
        if root in tree:
            return tree

        tree.append(root)
        for child in self._component_children(root):
            if child in tree:
                continue
            else:
                self.walk(child, tree=tree)
        return tree
