import warnings
from functools import wraps


def make_return_array_immutable(func):
    @wraps(func)
    def _wrapped(self, *args, **kwds):
        array = func(self, *args, **kwds)
        immutable_array = array.view()
        immutable_array.flags.writeable = False
        return immutable_array
    return _wrapped


def deprecated(func):
    """Mark a function as deprecated
    """
    @wraps(func)
    def _wrapped(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    _wrapped.__dict__.update(func.__dict__)

    return _wrapped
