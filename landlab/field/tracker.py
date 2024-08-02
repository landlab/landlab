from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable
from collections.abc import Generator
from collections.abc import Iterable
from contextlib import contextmanager

import numpy as np
from numpy.typing import NDArray

from landlab.field.errors import BadOperationError
from landlab.field.errors import FieldError
from landlab.field.errors import MissingTrackedFieldError
from landlab.field.errors import NothingToTrackError
from landlab.grid.base import ModelGrid


@contextmanager
def open_tracker(
    grid: ModelGrid,
    include: str | Iterable[str] = "*",
    exclude: str | Iterable[str] | None = None,
) -> Generator[FieldTracker, None, None]:
    """Track fields over a block of code.

    Parameters
    ----------
    include : str, or iterable of str, optional
        Glob-style pattern for field names to include.
    exclude : str, or iterable of str, optional
        Glob-style pattern for field names to exclude.

    Returns
    -------
    tuple of str
        The canonical names of the fields being tracked.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.field.tracker import open_tracker

    >>> grid = RasterModelGrid((3, 4))
    >>> _ = grid.add_ones("discharge", at="cell")

    >>> with open_tracker(grid) as tracker:
    ...     grid.at_cell["discharge"] += 4.0
    ...

    >>> tracker.fields
    (('at_cell:discharge', [array([1., 1.]), array([5., 5.])]),)
    """
    tracker = FieldTracker(grid)
    tracker.open(include=include, exclude=exclude)
    try:
        yield tracker
    finally:
        tracker.close()


class FieldTracker:
    """Track a grid's fields.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab ModelGrid

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.field.tracker import FieldTracker
    >>> from landlab.field.tracker import open_tracker

    >>> grid = RasterModelGrid((3, 4))
    >>> _ = grid.add_ones("z", at="node")

    ``open_tracker`` is a context manager that can be used to
    track fields. The following tracks all *at-node* fields.

    >>> with open_tracker(grid, "at_node*") as tracker:
    ...     grid.at_node["z"] += 2.0
    ...
    >>> tracker.fields
    (('at_node:z',
      [array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
       array([3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3.])]),)

    If you prefer, you can also use the ``open`` and ``close`` methods
    for tracking.

    >>> tracker = FieldTracker(grid)
    >>> tracker.open()
    ('at_node:z',)
    >>> grid.at_node["z"] += 6.0
    >>> tracker.close()
    >>> tracker.fields
    (('at_node:z',
      [array([3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3.]),
       array([9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9.])]),)

    Use the ``checkpoint`` method if you would like to save some intermediate
    values.

    >>> _ = grid.add_ones("temperature", at="node")
    >>> with open_tracker(grid, "*temp*") as tracker:
    ...     grid.at_node["temperature"] += 4.0
    ...     tracker.checkpoint()
    ...     grid.at_node["temperature"] += 1.0
    ...

    >>> tracker.fields
    (('at_node:temperature',
      [array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
       array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.]),
       array([6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6.])]),)

    You can reduce the tracked fields to a single array using the ``reduce``
    method. The `op` keyword is a function that accepts two array parameters
    and returns a new array.

    >>> tracker.reduce(op=np.subtract)
    (('at_node:temperature', array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.])),)
    """

    def __init__(self, grid: ModelGrid) -> None:
        self._grid = grid
        self._fields: dict[str, list[NDArray]] = defaultdict(list)
        self._tracking: tuple[str, ...] | None = None

    @property
    def tracking(self) -> tuple[str, ...] | None:
        """Return the name of the fields being tracked."""
        return self._tracking

    @property
    def fields(self) -> tuple[tuple[str, list[NDArray]], ...]:
        """Return the tracked field values."""
        return tuple(sorted(self._fields.items()))

    def open(
        self,
        include: str | Iterable[str] = "*",
        exclude: str | Iterable[str] | None = None,
    ) -> tuple[str, ...]:
        """Start tracking a set of fields.

        Parameters
        ----------
        include : str, or iterable of str, optional
            Glob-style pattern for field names to include.
        exclude : str, or iterable of str, optional
            Glob-style pattern for field names to exclude.

        Returns
        -------
        tuple of str
            The canonical names of the fields being tracked.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.field.tracker import FieldTracker

        >>> grid = RasterModelGrid((3, 4))
        >>> _ = grid.add_ones("elevation", at="node")
        >>> _ = grid.add_ones("temperature", at="node")
        >>> _ = grid.add_ones("discharge", at="link")

        >>> tracker = FieldTracker(grid)
        >>> tracker.open(include=("at_link*", "at_node*"), exclude="*temperature*")
        ('at_link:discharge', 'at_node:elevation')
        >>> tracker.fields
        (('at_link:discharge',
          [array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])]),
         ('at_node:elevation',
          [array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])]))
        """
        names = self._grid.fields(include=include, exclude=exclude)
        if not names:
            raise NothingToTrackError()

        self._fields.clear()

        self._tracking = tuple(sorted(names))
        self.checkpoint()

        return self._tracking

    def close(self) -> None:
        """Stop tracking."""
        if self.tracking is None:
            return

        try:
            self.checkpoint()
        finally:
            self._tracking = None

    def checkpoint(self) -> None:
        """Add current fields values to the tracker.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.field.tracker import open_tracker

        >>> grid = RasterModelGrid((3, 4))
        >>> _ = grid.add_ones("elevation", at="node")
        >>> _ = grid.add_ones("discharge", at="link")

        >>> with open_tracker(grid) as tracker:
        ...     grid.at_node["elevation"] += 2.0
        ...     tracker.checkpoint()
        ...     grid.at_node["elevation"] += 3.0
        ...     grid.at_link["discharge"] += 4.0
        ...
        >>> tracker.fields
        (('at_link:discharge',
          [array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
           array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
           array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.])]),
         ('at_node:elevation',
          [array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
           array([3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3.]),
           array([6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6.])]))
        """
        if self.tracking is None:
            raise BadOperationError(
                "checkpoint", reason="The tracker has not been started."
            )

        fields = FieldTracker._copy_fields(self._grid, self.tracking)

        for name, values in fields:
            self._fields[name].append(values)

    def reduce(
        self, op: Callable[[NDArray, NDArray], NDArray] | None = None
    ) -> tuple[tuple[str, NDArray], ...]:
        """Reduced the tracked fields to a single array.

        Parameters
        ----------
        op : func, optional
            A function used to reduce the saved field values. The function
            should take two arguments, both arrays, and return a new array.

        Returns
        -------
        tuple of (str, ndarray)
            The reduced fields.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.field.tracker import open_tracker

        >>> grid = RasterModelGrid((3, 4))
        >>> _ = grid.add_ones("elevation", at="node")
        >>> _ = grid.add_ones("discharge", at="cell")

        >>> with open_tracker(grid) as tracker:
        ...     grid.at_node["elevation"] += 2.0
        ...     tracker.checkpoint()
        ...     grid.at_node["elevation"] += 3.0
        ...     grid.at_cell["discharge"] += 4.0
        ...

        >>> tracker.reduce(np.subtract)
        (('at_cell:discharge',
          array([4., 4.])),
         ('at_node:elevation',
          array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.])))

        >>> def foo(x, y):
        ...     return (y - x) ** 2
        ...
        >>> tracker.reduce(foo)
        (('at_cell:discharge',
          array([16., 16.])),
         ('at_node:elevation',
          array([25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25.])))
        """
        return FieldTracker._reduce_fields(self._fields, op=op)

    @staticmethod
    def _copy_fields(grid, names: Iterable[str]) -> list[tuple[str, NDArray]]:
        """Copy field values from a grid."""
        fields = []
        for full_name in names:
            at, name = full_name.split(":")
            try:
                fields.append((full_name, getattr(grid, at)[name].copy()))
            except FieldError:
                raise MissingTrackedFieldError(full_name) from None

        return fields

    @staticmethod
    def _reduce_fields(
        fields: dict[str, list[NDArray]],
        op: Callable[[NDArray, NDArray], NDArray] | None = None,
    ):
        """Reduce field values to a single array."""
        op = np.subtract if op is None else op
        reduced_values = []
        for name, values in fields.items():
            reduced_values.append((name, op(values[-1], values[0])))
        return tuple(sorted(reduced_values))
