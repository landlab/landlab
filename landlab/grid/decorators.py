"""This module defines decorators used with ModelGrid objects."""
import numpy as np


class override_array_setitem_and_reset(object):

    """Decorator that calls a grid method after setting array values.

    This decorator wraps :any:`ModelGrid` methods that return a numpy array
    so that it returns a wrapped array that overrides the numpy array
    `__setitem__`, `__setslice__`, and `itemset` methods. The wrapped methods
    set values in the array but then also call a grid method that resets some
    state variables of the grid.

    Parameters
    ----------
    reset : str
        The name of the grid method to call after setting values. The
        corresponding method must take no arguments.
    """

    def __init__(self, reset):
        """Initialize the decorator with an argument.

        Parameters
        ----------
        reset : str
            The name of the grid method to call after setting values. The
            corresponding method must take no arguments.
        """
        self._reset = reset

    def __call__(self, func):
        """Get a wrapped version of the method.

        Parameters
        ----------
        func : function
            The grid method to wrap.

        Returns
        -------
        function
            The wrapped function.
        """
        reset = self._reset
        def _wrapped(grid):
            """Embed a grid into a numpy array and override set methods."""
            class array(np.ndarray):

                """Override numpy setters and reset grid topology."""

                def __new__(cls, arr):
                    """Instantiate the class with a view of the base array."""
                    obj = np.asarray(arr).view(cls)
                    obj.grid = grid
                    return obj

                def __array_finalize__(self, obj):
                    if obj is None:
                        return

                def itemset(self, ind, value):
                    """Set value of array, then call reset function."""
                    np.ndarray.itemset(self, ind, value)
                    getattr(self.grid, reset)()

                def __setitem__(self, ind, value):
                    """Set value of array, then call reset function."""
                    np.ndarray.__setitem__(self, ind, value)
                    getattr(self.grid, reset)()

                def __setslice__(self, start, stop, value):
                    """Set values of array, then call reset function."""
                    np.ndarray.__setslice__(self, start, stop, value)
                    getattr(self.grid, reset)()

            return array(func(grid))

        _wrapped.__name__ = func.__name__
        _wrapped.__doc__ = func.__doc__

        return _wrapped


