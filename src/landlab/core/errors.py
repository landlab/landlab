class Error(Exception):
    """Base class for exceptions raised from this module."""


class MissingKeyError(Error):
    """Error to indicate a missing parameter key.

    Raise this error if the parameter dictionary file does not contain a
    requested *key*.
    """


class ParameterValueError(Error):
    """Error to indicate a bad parameter values.

    Raise this error if a parameter value given by *key* is not of the
    expected type.
    """

    def __str__(self) -> str:
        key, val, expected_type = self.args
        return f"{key!r}: {val!r} is not of type {expected_type!r}"


class ValidationError(Error):
    """Error to indicate that a validation has failed."""
