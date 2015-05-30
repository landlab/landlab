from functools import wraps

import six


def track_this_method(func):
    @wraps(func)
    def _wrapped(self, *args, **kwds):
        if self._DEBUG_TRACK_METHODS:
            six.print_('Entering: %s.%s' % (self.__class__.__name__,
                                            func.__name__))
        ans = func(self, *args, **kwds)
        if self._DEBUG_TRACK_METHODS:
            six.print_('Exiting: %s.%s' % (self.__class__.__name__,
                                           func.__name__))
        return ans
    return _wrapped
