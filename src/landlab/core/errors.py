class Error(Exception):
    """Base class for exceptions raised from this module."""

    pass


class MissingKeyError(Error):
    """Error to indicate a missing parameter key.

    Raise this error if the parameter dictionary file does not contain a
    requested *key*.
    """

    def __init__(self, key):
        self._key = key

    def __str__(self):
        return self._key


class ParameterValueError(Error):
    """Error to indicate a bad parameter values.

    Raise this error if a parameter value given by *key* is not of the
    expected type.
    """

    def __init__(self, key, val, expected_type):
        self._key = key
        self._val = val
        self._type = expected_type

    def __str__(self):
        return f"{self._key}: {self._val} is not of type {self._type}"
