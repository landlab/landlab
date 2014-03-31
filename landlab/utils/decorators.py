import warnings


def make_return_array_immutable(func):
    def _wrapped(self, *args, **kwds):
        array = func(self, *args, **kwds)
        immutable_array = array.view()
        immutable_array.flags.writeable = False
        return immutable_array
    return _wrapped


def deprecated(func):
    """Mark a function as deprecated
    """
    def _wrapped(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    _wrapped.__name__ = func.__name__
    _wrapped.__doc__ = func.__doc__
    _wrapped.__dict__.update(func.__dict__)

    return _wrapped
