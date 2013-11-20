
def make_return_array_immutable(func):
    def _wrapped(self, *args, **kwds):
        array = func(self, *args, **kwds)
        immutable_array = array.view()
        immutable_array.flags.writeable = False
        return immutable_array
    return _wrapped
