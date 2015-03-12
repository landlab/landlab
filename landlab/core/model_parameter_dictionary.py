#! /usr/env/python
"""
The Model Parameter Dictionary is a tool for numerical modelers to
easily read and access model parameters from a simple formatted
input (text) file. Each parameter has a KEY, which identifies the
parameter, and a VALUE, which can be a number or a string. A
ModelParameterDictionary object reads model parameters from an input
file to a Dictionary, and provides functions for the user to look up
particular parameters by key name.

The format of the input file looks like::

    >>> from StringIO import StringIO
    >>> param_file = StringIO('''
    ... PI: the text "PI" is an example of a KEY
    ... 3.1416
    ... AVOGADROS_NUMBER: this is another
    ... 6.022e23
    ... FAVORITE_FRUIT: yet another
    ... mangoes
    ... NUMBER_OF_MANGO_WALKS: this one is an integer
    ... 4
    ... ALSO_LIKES_APPLES: this is a boolean
    ... true
    ... ''')

Example code that reads these parameters from a file called
*myinputs.txt*:

    >>> from landlab import ModelParameterDictionary
    >>> my_param_dict = ModelParameterDictionary()
    >>> my_param_dict.read_from_file(param_file)
    >>> pi = my_param_dict.read_float('PI')
    >>> avogado = my_param_dict.read_float('AVOGADROS_NUMBER')
    >>> fruit = my_param_dict.read_string('FAVORITE_FRUIT')
    >>> nmang = my_param_dict.read_int('NUMBER_OF_MANGO_WALKS')
    >>> apples_ok = my_param_dict.read_bool('ALSO_LIKES_APPLES')

As in Python, hash marks (#) denote comments. The rules are that each
key must have one and only one parameter value, and each value must
appear on a separate line immediately below the key line.

Also available are functions to read input parameters from the
command line (e.g., read_float_cmdline( 'PI' ) )
"""

# Licensing information:
#
# model_parameter_dictionary.py: reads formatted input from a text file
#   for use in specifying parameter values in numerical models.
#
# Copyright (C) 2011 Gregory E. Tucker
#
# Developer can be contacted at:
#   Cooperative Institute for Research in Environmental Sciences (CIRES)
#   University of Colorado Boulder
#   Campus Box 399
#   Boulder, CO 80309 USA
#   Phone: (+1) 303-492-6985
#   Email: gtucker@colorado.edu
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA.


import warnings
import types
 
import scipy.io
import numpy as np


_VALID_TRUE_VALUES = set(['TRUE', '1', 1])
_VALID_FALSE_VALUES = set(['FALSE', '0', 0])
_VALID_BOOLEAN_VALUES = _VALID_TRUE_VALUES | _VALID_FALSE_VALUES


class Error(Exception):
    """
    Base class for exceptions raised from this module.
    """
    pass


class MissingKeyError(Error):
    """
    Raise this error if the parameter dictionary file does not contain a
    requested *key*.
    """
    def __init__(self, key):
        self._key = key

    def __str__(self):
        return self._key


class ParameterValueError(Error):
    """
    Raise this error if a parameter value given by *key* is not of the
    expected type.
    """
    def __init__(self, key, val, expected_type):
        self._key = key
        self._val = val
        self._type = expected_type

    def __str__(self):
        return '%s: %s is not of type %s' % (self._key, self._val, self._type)


def _to_bool(string):
    upper = string.upper()
    if upper in _VALID_TRUE_VALUES:
        return True
    elif upper in _VALID_FALSE_VALUES:
        return False
    else:
        raise ValueError('invalid literal for _to_bool()')


_CONVERT_FROM_STR = {
    'float': float,
    'int': int,
    'str': str,
    'bool': _to_bool,
}


_VALID_VALUE_TYPES = set(_CONVERT_FROM_STR)

class ModelParameterDictionary(dict):
    """
    Reads model parameters from an input file to a dictionary
    and provides functions for the user to look up particular parameters
    by key name.

    If the keyword *auto_type* is True, then guess at the type for each value.
    Use *from_file* to read in a parameter file from a file-like object or a
    file with the given file name.

    Example
    =======

    Create a file-like object that contains a model parameter dictionary.

    >>> from StringIO import StringIO
    >>> test_file = StringIO('''
    ... INT_VAL:
    ... 1
    ... DBL_VAL:
    ... 1.2
    ... BOOL_VAL:
    ... true
    ... INT_ARRAY: 
    ... 1,2,3
    ... DBL_ARRAY: 
    ... 1.,2.,3.
    ... STR_VAL:
    ... landlab is awesome!
    ... ''')

    Create a ModelParameterDictionary, fill it with values from the
    parameter dictionary, and try to convert each value string to its
    intended type.

    >>> from landlab import ModelParameterDictionary
    >>> params = ModelParameterDictionary(auto_type=True, from_file=test_file)

    The returned ModelParameterDictionary can now be used just like a
    regular Python dictionary to get items, keys, etc.

    >>> print sorted(params.keys())
    ['BOOL_VAL', 'DBL_ARRAY', 'DBL_VAL', 'INT_ARRAY', 'INT_VAL', 'STR_VAL']

    >>> print params['INT_VAL']
    1
    >>> print params['DBL_VAL']
    1.2
    >>> print params['BOOL_VAL']
    True
    >>> print params['STR_VAL']
    landlab is awesome!

    Lines containing commas are converted to numpy arrays. The type of the
    array is determined by the values.

    >>> isinstance(params['DBL_ARRAY'], np.ndarray)
    True
    >>> print params['INT_ARRAY']
    [1 2 3]
    >>> print params['DBL_ARRAY']
    [ 1.  2.  3.]
    """
    def __init__(self, from_file=None, auto_type=False):
        super(ModelParameterDictionary, self).__init__()

        self._auto_type = auto_type
        if from_file is not None:
            self.read_from_file(from_file)

    def read_from_file(self, param_file):
        """
        Read and parse parameter dictionary information from a file or
        file-like object in *param_file*.

        The format of the parameter file should be like::

            # A comment line
            SOME_KEY: this the key for some parameter
            1.234

        In other words, the rules are:
            - Comments are preceded by hash characters
            - Each parameter has two consecutive lines, one for the
                key and one for the value
            - The key must be followed by a space, colon, or eol
            - The parameter can be numeric or text
        """
        if isinstance(param_file, types.StringTypes):
            try:
                with open(param_file, 'r') as opened_file:
                    self._read_from_file_like(opened_file)
            except IOError:
                raise
        else:
            self._read_from_file_like(param_file)

    @staticmethod
    def _get_stripped_lines(param_file):
        """
        Strip lines from the iterable, *param_file*. Ignore lines, that upon
        being stripped, are either blank or only contain a comment.
        Comments are lines that contain a hash preceeded only by whitespace.

        Returns a list of the stripped lines.
        """
        stripped_line_list = []
        for line in param_file:
            line = line.strip()   # strip leading spaces
            if len(line) > 0 and line[0] != '#':
                stripped_line_list.append(line)
        return stripped_line_list

    def _read_from_file_like(self, param_file):
        """
        Read parameters from the file-like object, *param_file*. In fact, 
        *param_file* really only needs to be an iterable of strings.
        """
        stripped_line_list = self._get_stripped_lines(param_file)

        iskey = True
        for line in stripped_line_list:
            if iskey:
                # Strip out everything after the first space or colon
                first_colon = line.find(':')
                if first_colon == -1:
                    first_colon = len(line)
                first_space = line.find(' ')
                if first_space == -1:
                    first_space = len(line)
                last_key_char = min(first_colon, first_space)
                last_key = line[0:last_key_char]  # remember the key
                iskey = False
            else:
                if self._auto_type:
                    self[last_key] = self._auto_type_value(line)
                else:
                    self[last_key] = line
                iskey = True

    @staticmethod
    def _auto_type_value(line):
        """
        Read a value from a string and try to guess its type. If the line
        contains any commas, try to convert it to a numpy array of ints and,
        if that doesn't work, an array of floats, otherwise it's just a
        string. 

        If there are no commas, the order of types is bool, int, float, and
        str.
        """
        if ',' in line:
            try:
                return np.array(line.split(','), np.int)
            except ValueError:
                try:
                    return np.array(line.split(','), np.float)
                except ValueError:
                    return line
        else:
            if line.upper() in _VALID_BOOLEAN_VALUES:
                return line.upper() in _VALID_TRUE_VALUES
            else:
                try:
                    return int(line)
                except ValueError:
                    try:
                        return float(line)
                    except ValueError:
                        return line

    def params(self):
        """
        Return a list of all the parameters names in the parameter dictionary.
        """
        return self.keys()

    def get(self, key, *args, **kwds):
        """get(key, [default], ptype=str)

        Get a value from a model parameter dictionary. Use the *ptype*
        keyword to convert the value to a given type. *ptype* is a function
        that converts the retreived value to the desired value. If a second
        argument after *key* is provided, use it as a default in case *key*
        is not contained in the ModelParameterDictionary.

        >>> from StringIO import StringIO
        >>> from landlab import ModelParameterDictionary
        >>> params = ModelParameterDictionary(StringIO(
        ... '''
        ... MY_INT:
        ... 1
        ... '''))
        >>> params.get('MY_INT')
        '1'
        >>> params.get('MY_INT', ptype=int)
        1
        >>> params.get('MY_MISSING_INT', 2, ptype=int)
        2

        Be careful when dealing with booleans. If you want to be returned
        a boolean value, *DO NOT* set the *ptype* keyword to the builtin
        *bool*. This will not work as the Python *bool* function does not
        convert strings to booleans as you might expect. For example,

        >>> bool('True')
        True
        >>> bool('False')
        True

        If you would like to get a boolean, use *ptype='bool'*.

        .. note:: Use *ptype='bool'* not *ptype=bool*.
            If you use *bool* to convert a string the returned boolean will
            be `True` for *any* non-empty string. This is just how the
            Python built-in *bool* works,

            >>> bool('0')
            True
            >>> bool('1')
            True
            >>> bool('')
            False

        >>> from StringIO import StringIO
        >>> params = ModelParameterDictionary(StringIO(
        ... '''
        ... MY_BOOL:
        ... false
        ... '''))
        >>> params.get('MY_BOOL')
        'false'
        >>> params.get('MY_BOOL', ptype='bool')
        False
        """
        ptype = kwds.pop('ptype', str)
        assert(len(kwds) == 0)

        if ptype is bool:
            warnings.warn(
                "Using bool function to convert value. This probably isn't "
                "what you want. Use ptype='bool' instead.",
                RuntimeWarning)

        value = super(ModelParameterDictionary, self).get(key, *args)
        if value is None:
            raise MissingKeyError(key)

        if isinstance(ptype, str):
            try:
                converter = _CONVERT_FROM_STR[ptype]
            except KeyError:
                converter = str
        else:
            converter = ptype

        try:
            typed_value = converter(value)
        except ValueError:
            raise ParameterValueError(key, value, ptype)
        else:
            return typed_value


    def read_int(self, key, *args):
        """
        Locate *key* in the input file and return it as an integer.

        >>> from StringIO import StringIO
        >>> from landlab import ModelParameterDictionary
        >>> params = ModelParameterDictionary(StringIO(
        ... '''
        ... MY_INT:
        ... 1
        ... '''))
        >>> params.read_int('MY_INT')
        1

        Raise an error if *key* isn't in the dictionary or if its value is
        not an integer.
        """
        return self.get(key, *args, ptype=int)

    def read_float(self, key):
        """
        Locate *key* in the input file and return it as a float.

        >>> from StringIO import StringIO
        >>> from landlab import ModelParameterDictionary
        >>> params = ModelParameterDictionary(StringIO(
        ... '''
        ... MY_FLOAT:
        ... 3.14
        ... '''))
        >>> print round(params.read_float('MY_FLOAT'), 6)
        3.14

        An error is generated if *key* isn't in the dictionary or
        if its value is not a number.
        """
        try:
            my_float = float(self[key])
        except KeyError:
            raise MissingKeyError(key)
        except ValueError:
            raise ParameterValueError(key, self[key], 'float')
        else:
            return my_float

    def read_string(self, key):
        """
        Locate *key* in the input file and return it as a string.

        >>> from StringIO import StringIO
        >>> from landlab import ModelParameterDictionary
        >>> params = ModelParameterDictionary(StringIO(
        ... '''
        ... MY_STRING:
        ... landlab
        ... '''))
        >>> params.read_string('MY_STRING')
        'landlab'

        An error is generated if *key* isn't in the dictionary.
        """
        try:
            my_value = self[key]
        except KeyError:
            raise MissingKeyError(key)
        return str(my_value)

    def read_bool(self, key):
        """
        >>> from StringIO import StringIO
        >>> from landlab import ModelParameterDictionary
        >>> params = ModelParameterDictionary(StringIO(
        ... '''
        ... MY_BOOL:
        ... true
        ... '''))
        >>> params.read_bool('MY_BOOL')
        True

        An error is generated if MY_BOOL isn't 0, 1, True or False
        """
        try:
            my_value = self[key]
        except KeyError:
            raise MissingKeyError(key)

        if my_value.upper() in _VALID_TRUE_VALUES:
            return True
        elif my_value.upper() in _VALID_FALSE_VALUES:
            return False
        else:
            raise ParameterValueError(key, my_value, 'boolean')

    def read_int_cmdline(self, key):
        """
        Read an integer from the command line and use it as a value for
        *key* in the dictonary.

        Usage: i = read_int_cmdline('MY_INT')

        An error is generated if *key* is not an integer.
        """
        my_value = input(key + ': ')
        self[key] = my_value
        if not isinstance(my_value, int):
            raise ParameterValueError(key, my_value, 'int')
        return my_value

    def read_float_cmdline(self, key):
        """
        Read a float from the command line and use it as a value for
        *key* in the dictonary.

        Usage: f = read_float_cmdline('MY_FLOAT')

        An error is generated if *key* is not a float.
        """
        my_value = input(key + ': ')
        self[key] = my_value
        try:
            my_float = float(my_value)
        except ValueError:
            raise ParameterValueError(key, my_value, 'float')
        else:
            return my_float

    def read_string_cmdline(self, key):
        """
        Read a string from the command line and use it as a value for
        *key* in the dictonary.

        """
        my_str = raw_input(key + ': ')
        self[key] = my_str
        return my_str

    def read_bool_cmdline(self, key):
        """
        Read a boolean from the command line and use it as a value for
        *key* in the dictonary.

        Usage: f = read_bool_cmdline('MY_BOOL')

        An error is generated if *key* is not a boolean.
        """
        my_value = raw_input(key + ': ')
        if my_value.upper() in _VALID_TRUE_VALUES:
            return True
        elif my_value.upper() in _VALID_FALSE_VALUES:
            return False
        else:
            raise ParameterValueError(key, my_value, 'boolean')
