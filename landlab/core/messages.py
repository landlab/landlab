#! /usr/bin/env python
"""Print user messages formatted landlab-style.

This module provides functions for printing nicely-formatted
messages to the user.  Messages are formatted in a particular
style so that all of landlab messages will have a similar
look. Anytime landlab prints something for an end-user to see,
this module should be used.

This module also provides convenience functions for print
particular types of messages. Warning and error messages,
for instance.

Examples
--------
>>> from __future__ import print_function

Oftentimes when writing code we may need to print a lengthy
message for the user. This may result in code that looks like
the following.

>>> message = ('Lorem ipsum dolor sit amet, consectetur '
...            'adipiscing elit, sed do eiusmod tempor '
...            'incididunt ut labore et dolore magna aliqua. '
...            'Ut enim ad minim veniam, quis nostrud exercitation '
...            'ullamco laboris nisi ut aliquip ex ea commodo '
...            'consequat.')

Printing this message string would result in one long line that
would, most likely, extend beyond the user's terminal and be
difficult to read. One solution would be to join the lines
by line separators but then that would result in a bunch of really
short lines.

To help with this, landlab provides a set of functions with the
most basic being `format_message`.

>>> from landlab.core.messages import format_message
>>> print(format_message(message))
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad
minim veniam, quis nostrud exercitation ullamco laboris nisi ut
aliquip ex ea commodo consequat.

landlab also provides functions for printing warning and error
messages.


>>> from landlab.core.messages import warning_message
>>> message = ('Lorem ipsum dolor sit amet, consectetur\\n'
...            'adipiscing elit, sed do eiusmod tempor\\n'
...            'incididunt ut labore et dolore magna aliqua.')
>>> print(warning_message(message))
WARNING
=======
<BLANKLINE>
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
eiusmod tempor incididunt ut labore et dolore magna aliqua.

>>> from landlab.core.messages import error_message
>>> print(error_message(message))
ERROR
=====
<BLANKLINE>
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
eiusmod tempor incididunt ut labore et dolore magna aliqua.

Another common design pattern used in landlab is to allow the
user, usually through a keyword, to control what happens
if a particular assertion fails. For instance, the user
may want the code to raise an error, or print a warning
message, or do nothing at all. The `assert_or_print` function
should be used in these cases.

>>> from landlab.core.messages import assert_or_print

>>> dt = 1e6
>>> assert_or_print(dt < 1, 'Unstable time step!', onerror='pass')
>>> assert_or_print(dt < 1, 'Unstable time step!', onerror='warn')
...     #doctest: +SKIP
WARNING
=======

Unstable time step!

>>> assert_or_print(dt < 1, 'Unstable time step!', onerror='raise')
...     #doctest: +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
AssertionError
"""
from __future__ import print_function

import os
import re
import sys
import textwrap

import six


def indent_and_wrap(content, indent=""):
    """Indent and wrap some text

    Lines are first dedented to remove common leading whitespace,
    then indented according to the value of *indent*, and then
    wrapped at 70 characters (indenting if necessary with subsequent
    indent being twice *indent*).

    Note that when looking for common whitespace, the first line is
    ignored.

    Parameters
    ----------
    content : str
        The content to wrap.

    Returns
    -------
    str
        The content properly wrapped and indented.

    Examples
    --------
    >>> from __future__ import print_function
    >>> from landlab.core.messages import indent_and_wrap
    >>> content = '''@book{knuth1998art,
    ...     title={The art of computer programming: sorting and searching},
    ...     author={Knuth, Donald Ervin},
    ...     volume={3},
    ...     year={1998},
    ...     publisher={Pearson Education}
    ...     }'''
    >>> print(indent_and_wrap(content))
    @book{knuth1998art,
    title={The art of computer programming: sorting and searching},
    author={Knuth, Donald Ervin},
    volume={3},
    year={1998},
    publisher={Pearson Education}
    }
    """
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=2 * indent)
    lines = content.splitlines()
    first_line, the_rest = [lines[0].strip()], lines[1:]
    if the_rest:
        the_rest = textwrap.dedent(os.linesep.join(the_rest)).splitlines()
    if first_line[0]:
        lines = first_line + the_rest
    else:
        lines = the_rest
    return os.linesep.join([os.linesep.join(wrapper.wrap(line)) for line in lines])


def split_paragraphs(msg, linesep=os.linesep):
    """Split text into paragraphs.

    Split a block of text into paragraphs. A paragraph is
    defined as adjacent new-line characters (possibly separated
    by some whitespace).

    Parameters
    ----------
    msg : str
        Text to split into paragraphs.
    linesep : str, optional
        Line separator used in the message string.

    Returns
    -------
    list of str
        List of paragraphs.

    Examples
    --------
    >>> from landlab.core.messages import split_paragraphs
    >>> text = '''
    ... Pharetra pharetra massa massa ultricies mi quis hendrerit.
    ...
    ... Dictumst vestibulum rhoncus est pellentesque.
    ... '''
    >>> split_paragraphs(text, linesep='\\n') #doctest: +NORMALIZE_WHITESPACE
    ['Pharetra pharetra massa massa ultricies mi quis hendrerit.',
     'Dictumst vestibulum rhoncus est pellentesque.']

    >>> text = '''
    ... Pharetra pharetra massa massa ultricies mi quis hendrerit.
    ... Dictumst vestibulum rhoncus est pellentesque.
    ... '''
    >>> len(split_paragraphs(text, linesep='\\n'))
    1
    """
    pattern = linesep + r"\s*" + linesep
    parsep = linesep * 2
    return re.sub(pattern, parsep, msg.strip()).split(parsep)


def format_message(msg, header=None, footer=None, linesep=os.linesep):
    """Format a message, landlab-style.

    Create a nicely formatted message that splits paragraphs,
    dedents paragraphs, and wraps text at 70 characters.
    Optionally, add a header and footer to the resulting message.


    Parameters
    ----------
    msg : str
        The message to be formatted.
    header : str or list of str, optional
        String to add before *msg*.
    footer : str or list of str, optional
        String to add after *msg*.

    Returns
    -------
    str
        The formatted message.

    Examples
    --------
    >>> from __future__ import print_function
    >>> from landlab.core.messages import format_message
    >>> text = '''
    ... Lorem ipsum dolor sit amet, consectetur
    ... adipiscing elit, sed do eiusmod tempor
    ... incididunt ut labore et dolore magna aliqua.
    ...
    ... Pharetra pharetra massa massa ultricies mi
    ... quis hendrerit.
    ...
    ... Dictumst vestibulum rhoncus est pellentesque.
    ... Sed viverra tellus in hac habitasse platea
    ... dictumst vestibulum rhoncus.'''
    >>> print(format_message(text))
    Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
    eiusmod tempor incididunt ut labore et dolore magna aliqua.
    <BLANKLINE>
    Pharetra pharetra massa massa ultricies mi quis hendrerit.
    <BLANKLINE>
    Dictumst vestibulum rhoncus est pellentesque. Sed viverra tellus in
    hac habitasse platea dictumst vestibulum rhoncus.
    """
    if isinstance(header, six.string_types):
        header = [header]
    header = header or []

    if isinstance(footer, six.string_types):
        footer = [footer]
    footer = footer or []

    paragraphs = header
    if msg is not None:
        for paragraph in split_paragraphs(msg.strip(), linesep=linesep):
            paragraphs.append(
                os.linesep.join(textwrap.wrap(textwrap.dedent(paragraph)))
            )
    paragraphs += footer

    return (os.linesep * 2).join(paragraphs)


def deprecation_message(msg=None, **kwds):
    """Create a deprecation message, landlab-style.

    Parameters
    ----------
    msg : str, optional
        Warning message.

    Returns
    -------
    str
        The formatted warning message.

    Examples
    --------
    >>> from __future__ import print_function
    >>> from landlab.core.messages import deprecation_message
    >>> print(deprecation_message("Dictumst vestibulum rhoncus est pellentesque."))
    DEPRECATION WARNING
    ===================
    <BLANKLINE>
    Dictumst vestibulum rhoncus est pellentesque.

    >>> print(
    ...     deprecation_message(
    ...         "Dictumst vestibulum rhoncus est pellentesque.",
    ...         use="Lorem ipsum dolor sit amet",
    ...     )
    ... )
    DEPRECATION WARNING
    ===================
    <BLANKLINE>
    Dictumst vestibulum rhoncus est pellentesque.
    <BLANKLINE>
    Example
    -------
    Lorem ipsum dolor sit amet
    """
    use = kwds.pop("use", None)
    if use:
        footer = os.linesep.join(["Example", "-------", use])
    else:
        footer = None
    header = "Deprecation warning".upper()
    return format_message(
        msg, header=os.linesep.join([header, "=" * len(header)]), footer=footer, **kwds
    )


def warning_message(msg=None, **kwds):
    """Create a warning message, landlab-style.

    Parameters
    ----------
    msg : str, optional
        Warning message.

    Returns
    -------
    str
        The formatted warning message.

    Examples
    --------
    >>> from __future__ import print_function
    >>> from landlab.core.messages import warning_message
    >>> print(warning_message('Dictumst vestibulum rhoncus est pellentesque.'))
    WARNING
    =======
    <BLANKLINE>
    Dictumst vestibulum rhoncus est pellentesque.
    """
    header = "Warning".upper()
    return format_message(
        msg, header=os.linesep.join([header, "=" * len(header)]), **kwds
    )


def error_message(msg=None, **kwds):
    """Create an error  message, landlab-style.

    Parameters
    ----------
    msg : str, optional
        Warning message.

    Returns
    -------
    str
        The formatted warning message.

    Examples
    --------
    >>> from __future__ import print_function
    >>> from landlab.core.messages import error_message
    >>> print(error_message('Dictumst vestibulum rhoncus est pellentesque.'))
    ERROR
    =====
    <BLANKLINE>
    Dictumst vestibulum rhoncus est pellentesque.
    """
    header = "Error".upper()
    return format_message(
        msg, header=os.linesep.join([header, "=" * len(header)]), **kwds
    )


def assert_or_print(cond, msg=None, onerror="raise", file=sys.stdout):
    """Make an assertion printing a message if it fails.

    Specify an action to take if an assertion fails, depending on
    the values of *onerror*. *onerror* must be one of:
        *  "pass": do nothing if the assertion passes or fails.
        *  "warn": print a warning message if the assertion fails.
        *  "error": print an error message and raise an `AssertionError`
           on failure.

    Parameters
    ----------
    cond : expression
        An expression to test.
    msg : str, optional
        Message to print if the condition is not met.
    onerror: {'pass', 'warn', 'raise'}, optional
        What to do if the condition evaluates to `False`.
    file : file_like, optional
        File-like object where the message is printed.

    Examples
    --------
    >>> from landlab.core.messages import assert_or_print

    >>> assert_or_print(True, 'Lorem ipsum', onerror='pass')
    >>> assert_or_print(False, 'Lorem ipsum', onerror='pass')

    >>> assert_or_print(True, 'Lorem ipsum', onerror='warn')
    >>> assert_or_print(False, 'Lorem ipsum', onerror='warn')

    >>> assert_or_print(True, 'Lorem ipsum', onerror='raise')
    >>> assert_or_print(False, 'Lorem ipsum', onerror='raise')
    ...     #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    AssertionError
    """
    if onerror not in ("pass", "warn", "raise"):
        raise ValueError("onerror must be one of 'pass', 'warn', or 'raise'")

    try:
        assert cond
    except AssertionError:
        if onerror == "warn":
            print(warning_message(msg), file=file, end="")
        elif onerror == "raise":
            print(error_message(msg), file=file, end="")
            raise
