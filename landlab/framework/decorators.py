#! /usr/bin/env python
"""Decorators for TheLandlab package."""
import re


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
