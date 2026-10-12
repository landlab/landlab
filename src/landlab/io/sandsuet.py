from collections.abc import Iterable
from typing import BinaryIO

import xarray as xr
from requireit import require_instance

from landlab.grid.raster import RasterModelGrid


def dump(
    grid: RasterModelGrid,
    *,
    stream: BinaryIO | None = None,
    include: str = "*",
    exclude: str | None = None,
) -> memoryview | None:
    """Dump grid fields to sandsuet format.

    Parameters
    ----------
    grid : RasterModelGrid
        The grid whose fields will be serialized.
    stream : BinaryIO, optional
        A binary stream to write to. If not provided, the serialized
        dataset is returned as a ``memoryview``.
    include : str, optional
        Glob-style pattern of fields to include (default is ``"*"``).
    exclude : str, optional
        Glob-style pattern of fields to exclude.

    Returns
    -------
    memoryview or None
        The serialized dataset as a ``memoryview`` if no stream is
        provided, otherwise ``None``.

    Raises
    ------
    ValueError
        If the selected fields are not all at the same, supported
        location (``node``, ``patch``, ``corner``, or ``cell``).

    Examples
    --------
    >>> import io
    >>> import xarray as xr
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 4), xy_spacing=10.0)
    >>> _ = grid.add_ones("topographic__elevation", at="node")
    >>> data = dump(grid)
    >>> ds = xr.open_dataset(data)
    >>> ds.attrs["sandsuet_version"]
    '1.0.0'
    >>> list(ds.data_vars)
    ['topographic__elevation']
    """
    grid = require_instance(grid, RasterModelGrid)

    fields = grid.fields(include=include, exclude=exclude)
    at = _validate_field_location(fields)

    if at is None:
        return None

    shape = grid.shape
    if at in ("patch", "corner"):
        shape = (shape[0] - 1, shape[1] - 1)
    elif at == "cell":
        shape = (shape[0] - 2, shape[1] - 2)

    xy = getattr(grid, f"xy_of_{at}")
    x_of_col = xy[:, 0].reshape(shape)[0, :]
    y_of_row = xy[:, 1].reshape(shape)[:, 0]

    coords = {
        "y": xr.DataArray(y_of_row, dims=("y",), name="y"),
        "x": xr.DataArray(x_of_col, dims=("x",), name="x"),
    }

    data_vars = {}
    for full_name in fields:
        _, _, field_name = full_name.partition(":")
        data_vars[field_name] = xr.DataArray(
            grid.field_values(field_name, at=at).reshape(shape),
            dims=("y", "x"),
            name=field_name,
        )

    dataset = xr.Dataset(
        data_vars=data_vars, coords=coords, attrs={"sandsuet_version": "1.0.0"}
    )

    if stream is None:
        return dataset.to_netcdf()
    else:
        dataset.to_netcdf(stream)
        return None


def _validate_field_location(names: Iterable[str]) -> str | None:
    valid_locations = ("node", "corner", "patch", "cell")

    at, bad_names = set(), set()
    for name in names:
        location, _, _ = name.removeprefix("at_").partition(":")
        if location not in valid_locations:
            bad_names.add(name)
        else:
            at.add(location)

    if bad_names or len(at) > 1:
        valid_str = ", ".join(repr(loc) for loc in sorted(valid_locations))
        raise ValueError(f"fields must be at one, and only one, of {valid_str}")

    return at.pop() if len(at) == 1 else None
