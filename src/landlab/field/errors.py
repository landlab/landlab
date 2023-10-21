class Error(Exception):

    """Base class for errors in this module."""

    pass


class FieldError(Error, KeyError):

    """Raise this error for a missing field name."""

    def __init__(self, field):
        self._field = field

    def __str__(self):
        return self._field


class GroupError(Error, KeyError):

    """Raise this error for a missing group name."""

    def __init__(self, group):
        self._group = group

    def __str__(self):
        return self._group


class GroupSizeError(Error, KeyError):

    """Raise this error if a group has changed sizes."""

    def __init__(self, group, old_size, new_size):
        self._group = group
        self._old_size = old_size
        self._new_size = new_size

    def __str__(self):
        return (
            "number of {group} elements has changed. "
            "(was = {was}, now={now})".format(
                group=self._group, was=self._old_size, now=self._new_size
            )
        )
