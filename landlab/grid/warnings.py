import os

from ..core.messages import deprecation_message


class DeprecatedSignature(DeprecationWarning):

    msg = "You are using a deprecated calling signature."

    def __init__(self, name, old=None, new=None):
        self._name = name
        self._old = old
        self._new = new

        if old:
            self._old = self._construct_call(name, self._old[0], self._old[1])
        if new:
            self._new = self._construct_call(name, self._new[0], self._new[1])

    @staticmethod
    def _construct_call(name, args, kwds):
        signature = ", ".join(
            [repr(arg) for arg in args]
            + ["{k}={v}".format(k=k, v=repr(v)) for k, v in kwds.items()]
        )
        return "{name}({signature})".format(name=name, signature=signature)

    def __str__(self):
        if self._new:
            use = ">>> grid = {call}".format(call=self._new)
        else:
            use = None

        return os.linesep + deprecation_message(self.msg, use=use)
