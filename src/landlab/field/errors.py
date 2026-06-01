class Error(Exception):
    """Base class for errors in this module."""


class FieldError(Error, KeyError):
    """Raise this error for a missing field name."""


class GroupError(Error, KeyError):
    """Raise this error for a missing group name."""


class GroupSizeError(Error, KeyError):
    """Raise this error if a group has changed sizes."""

    def __init__(self, group, old_size, new_size):
        super().__init__(group, old_size, new_size)

    def __str__(self):
        group, old_size, new_size = self.args
        return (
            f"number of {group} elements has changed."
            f" (was = {old_size}, now={new_size})"
        )
