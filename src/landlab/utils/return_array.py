#! /usr/bin/env python
"""Return array with same shape as grid elements."""

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab.field.errors import FieldError
from landlab.field.graph_field import GraphFields
from landlab.utils.decorators import use_field_name_array_or_value

type FieldLike = str | ArrayLike


def validate_field(
    value: FieldLike,
    *,
    grid: GraphFields | None = None,
    at: str | None = None,
) -> str | NDArray:
    """Validate a field name or array-like value.

    Check that *value* is either the name of an existing field, or an
    array-like object that could represent values located on a grid.
    Unlike :func:`resolve_field`, a field name is *not* resolved to its
    values, so the returned value can be stored and later passed to
    :func:`resolve_field` to pick up the field's current
    values at the time it's needed. *value* is converted to a *numpy*
    array, but, unlike :func:`resolve_field`, it is not broadcast.

    Parameters
    ----------
    value : str or array_like
        A field name, or an array of values.
    grid : GraphFields, optional
        Grid used to validate *value* against. If not given, *value* is
        not checked against a grid (and, if a field name, is returned
        unchecked).
    at : str, optional
        Name of the group (e.g. "node", "link") that *value* is defined
        on. Must be given if, and only if, *grid* is given.

    Returns
    -------
    str or ndarray
        *value* unchanged, if it is a field name, otherwise *value* as an
        array (not copied, if *value* is already an ndarray).

    Raises
    ------
    ValueError
        If only one of *grid* and *at* is given, or if *value* is
        array-like (and not a scalar) but its first dimension does not
        match the number of *at* elements of *grid*.
    landlab.field.errors.FieldError
        If *value* is a field name that does not exist in *grid* at *at*.

    Examples
    --------
    >>> from landlab.utils.return_array import validate_field
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 4))
    >>> _ = grid.add_ones("foo", at="node")

    >>> validate_field("foo", grid=grid, at="node")
    'foo'
    >>> validate_field(range(12), grid=grid, at="node")
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11])
    >>> validate_field(42, grid=grid, at="node")
    array(42)
    """
    if (grid is None) != (at is None):
        raise ValueError("grid and at must both be given, or neither")

    if isinstance(value, str):
        if grid is not None and value not in grid[at]:
            raise FieldError(f"{value!r} is not a field at {at!r}")
        return value

    value = np.asarray(value)

    if at is not None:
        size = grid.size(at)
        if value.ndim > 0 and value.shape[0] != size:
            raise ValueError(
                f"first dimension of value, {value.shape}, does not match the"
                f" number of {at!r} elements of grid, {size}"
            )

    return value


def resolve_field(
    field: FieldLike,
    *,
    grid: GraphFields,
    at: str,
) -> NDArray:
    """Resolve a field name or array to an array of values.

    If *field* is a field name, look up and return its current values on
    *grid*. If *field* is 0-dimensional, broadcast it to the number of
    *at* elements of *grid*. Otherwise, return *field* unchanged. Call
    this each time you need the values of a *field* previously validated with
    :func:`validate_field`, so that a field name always resolves to that
    field's current values, even if they have changed since *field* was
    validated.

    Parameters
    ----------
    field : str or array_like
        A field name, or an array of values, as returned by
        :func:`validate_field`.
    grid : GraphFields
        Grid to look up *field* on, or broadcast it against, as needed.
    at : str
        Name of the group (e.g. "node", "link") that *field* is defined
        on.

    Returns
    -------
    ndarray
        The current values of *field*. If *field* is 0-dimensional, this
        is a read-only, broadcast view rather than a newly allocated array.
        If *field* is already an *numpy* array, it is returned unchanged
        (not copied). Otherwise, a new array is created.

    Raises
    ------
    landlab.field.errors.GroupError
        If *at* is not a group on *grid*.
    landlab.field.errors.FieldError
        If *field* is a field name that does not exist in *grid* at *at*.

    Examples
    --------
    >>> from landlab.utils.return_array import resolve_field
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 4))
    >>> foo = grid.add_ones("foo", at="node")

    >>> resolve_field("foo", grid=grid, at="node") is foo
    True
    >>> resolve_field(range(12), grid=grid, at="node")
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11])
    >>> resolve_field(42, grid=grid, at="node")
    array([42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42])
    """
    if isinstance(field, str):
        return grid.field_values(field, at=at)

    field = np.asarray(field)

    return np.broadcast_to(field, (grid.size(at),)) if field.ndim == 0 else field


@use_field_name_array_or_value("node")
def return_array_at_node(grid, value):
    """Function to return an array stored at node or of shape `(n_nodes,)`.

    This function exists to take advantage of the use_field_name_array_or_value
    decorator which permits providing the surface as a field name or array.

    Parameters
    ----------
    grid : ModelGrid
    value : field name, ndarray of shape `(n_nodes, )`, or single value.

    Returns
    -------
    array : ndarray of shape `(n_nodes, )`
    """
    return value


@use_field_name_array_or_value("link")
def return_array_at_link(grid, value):
    """Function to return an array stored at node or of shape `(n_nodes,)`.

    This function exists to take advantage of the use_field_name_array_or_value
    decorator which permits providing the surface as a field name or array.

    Parameters
    ----------
    grid : ModelGrid
    value : field name, ndarray of shape `(n_nodes, )`, or single value.

    Returns
    -------
    array : ndarray of shape `(n_nodes, )`
    """
    return value
