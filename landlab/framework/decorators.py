#! /usr/bin/env python
"""
Decorators for TheLandlab package.
"""

import inspect
import re
import types

import six


def camel_case(text, sep=None):
    """Convert to camel case.

    Convert *text* to camel case. Use the *sep* keyword to specify the word
    separator. The default is to split on whitespace.

    >>> from landlab.framework.decorators import camel_case
    >>> camel_case("eric idle")
    'EricIdle'
    >>> camel_case("terry_gilliam", sep="_")
    'TerryGilliam'
    >>> camel_case("MONTY Python")
    'MONTYPython'
    >>> camel_case("GrahamChapman")
    'GrahamChapman'
    """
    return "".join([word[0].upper() + word[1:] for word in text.split(sep)])


def snake_case(text):
    """Convert camel case to snake case.

    Examples
    --------
    >>> from landlab.framework.decorators import snake_case
    >>> snake_case("EricIdle")
    'eric_idle'
    >>> snake_case("MONTYPython")
    'monty_python'
    """
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", text)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


class Error(Exception):
    """
    Exceptions for this module.
    """

    pass


class InterfaceImplementationError(Error):
    """
    Raise this error if the class does not implement an interface.
    """

    def __init__(self, cls, interface):
        self.cls = cls.__name__
        self.interface = interface.__name__

    def __str__(self):
        return "Class '%s' does not implement interface '%s'" % (
            self.cls,
            self.interface,
        )


def is_implementation(cls, interface):
    """
    Check if *cls* implements *interface*. A class implements the interface
    class if it has the same members, the members have the same type, and
    methods have the same signature.

    Returns ``True`` if *cls* implements *interface*, otherwise ``False``.
    """
    for (name, value) in inspect.getmembers(interface):
        if isinstance(value, types.MethodType):
            try:
                cls_args = inspect.getargspec(getattr(cls, name))
                interface_args = inspect.getargspec(value)
            except AttributeError:
                six.print_("Missing attribute %s" % name)
                return False
            try:
                assert len(cls_args.args) == len(interface_args.args)
            except AssertionError:
                six.print_("Mismatch in number of args for %s" % name)
                return False
        else:
            try:
                assert isinstance(getattr(cls, name), type(getattr(interface, name)))
                # assert(type(getattr(cls, name)) == type(getattr(interface, name)))
            except (AttributeError, AssertionError):
                six.print_("Missing member or type mismatch for %s" % name)
                return False
    return True


class ImplementsOrRaise(object):
    """
    Decorator to indicate if a class implements interfaces. If the class
    does not implement the interface, raise an InterfaceImplementationError.
    If the class does implement the interface, decorate it with a
    __implements__ data mamber that is a tuple of the interfaces it implements.
    """

    def __init__(self, *interfaces):
        self._interfaces = interfaces

    def __call__(self, cls):
        cls.__implements__ = ()
        for interface in self._interfaces:
            if is_implementation(cls, interface):
                cls.__implements__ += (interface.__name__,)
            else:
                raise InterfaceImplementationError(cls, interface)
        return cls


class Implements(object):
    """
    Decorator to indicate if a class implements interfaces. Similar to the
    ImplementsOrRaise decorator except that this decorator silently ignores
    implemention errors. If the class does implement the interface, decorate
    it with a __implements__ data mamber that is a tuple of the interfaces it
    implements. Otherwise, don't do anything.
    """

    def __init__(self, *interfaces):
        self._interfaces = interfaces

    def __call__(self, cls):
        cls.__implements__ = ()
        for interface in self._interfaces:
            if is_implementation(cls, interface):
                cls.__implements__ += (interface.__name__,)
            else:
                pass
        return cls
