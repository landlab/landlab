#! /usr/bin/env python
"""Read/write data from an ESRI ASCII file into a RasterModelGrid.

ESRI ASCII functions
++++++++++++++++++++

.. autosummary::

    ~dump
    ~lazy_load
    ~lazy_loads
    ~load
    ~loads
    ~parse
"""
from __future__ import annotations

import io
from collections import OrderedDict
from collections.abc import Callable
from typing import Literal
from typing import NamedTuple
from typing import TextIO
from typing import overload

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab.grid.raster import RasterModelGrid
from landlab.layers.eventlayers import _valid_keywords_or_raise


class EsriAsciiError(Exception):
    pass


class BadHeaderError(EsriAsciiError):
    def __init__(self, msg: str) -> None:
        self._msg = msg

    def __str__(self) -> str:
        return self._msg


def dump(
    grid: RasterModelGrid,
    stream: TextIO | None = None,
    at: Literal["node", "patch", "corner", "cell"] = "node",
    name: str | None = None,
    nodata_value: int | float = -9999,
) -> str | None:
    """Dump a grid field to ESRII ASCII format.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid.
    stream : file-like, optional
        A ``file_like`` object to write to. If ``None``, return
        a string containing the serialized grid.
    at : str, optional
        Where the field to be written is located on the grid.
    name : str, optional
        The name of the field to be written. If ``None``, only the header
        information will be written.
    nodata_value : float or int, optional
        The value to indicate there is no-data at this point.

    Returns
    -------
    str or None
        The grid field in ESRI ASCII format or ``None`` if a ``file_like``
        object was provided.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.io import esri_ascii

    >>> grid = RasterModelGrid((3, 4), xy_spacing=2.0)
    >>> print(esri_ascii.dump(grid, at="node"))
    NROWS 3
    NCOLS 4
    CELLSIZE 2.0
    XLLCENTER 0.0
    YLLCENTER 0.0
    NODATA_VALUE -9999

    >>> print(esri_ascii.dump(grid, at="cell", nodata_value=0))
    NROWS 1
    NCOLS 2
    CELLSIZE 2.0
    XLLCENTER 2.0
    YLLCENTER 2.0
    NODATA_VALUE 0
    """
    grid = _validate_grid(grid)

    shape = np.asarray(grid.shape)
    if at in ("corner", "patch"):
        shape -= 1
    elif at == "cell":
        shape -= 2

    shift = _get_lower_left_shift(at=at, ref="center")

    header = _Header(
        nrows=shape[0],
        ncols=shape[1],
        xllcenter=grid.xy_of_lower_left[0] - grid.dx * shift,
        yllcenter=grid.xy_of_lower_left[1] - grid.dx * shift,
        cellsize=grid.dx,
        nodata_value=nodata_value,
    )
    lines = [str(header)]

    if name:
        data = getattr(grid, f"at_{at}")[name]
        kwds = {"comments": ""}
        if data.dtype.kind in ("i", "u"):
            kwds["fmt"] = "%d"
        with io.StringIO() as fp:
            np.savetxt(fp, np.flipud(data.reshape(shape)), **kwds)
            lines.append(fp.getvalue())

    content = "".join(lines)
    if stream is None:
        return content
    else:
        stream.write(content)
        return None


def load(
    stream: TextIO,
    at: str = "node",
    name: str | None = None,
    out: RasterModelGrid | None = None,
) -> RasterModelGrid:
    """Parse a RasterModelGrid from a ESRI ASCII format stream.

    Parameters
    ----------
    stream : file_like
        A text stream in ESRI ASCII format.
    at : {'node', 'patch', 'corner', 'cell'}, optional
        Location on the grid where data are placed.
    name : str, optional
        Name of the newly-created field. If `name` is
        not provided, the grid will be created but the
        data not added as a field.
    out : RasterModelGrid, optional
        Place the data from the file onto an existing grid. If not
        provided, create a new grid.

    Returns
    -------
    :
        A newly-created ``RasterModelGrid`` with, optionaly, the data added
        as a field (if `name` was provided).
    """
    return loads(stream.read(), at=at, name=name, out=out)


def loads(
    s: str,
    at: str = "node",
    name: str | None = None,
    out: RasterModelGrid | None = None,
) -> RasterModelGrid:
    """Parse a RasterModelGrid from a ESRI ASCII formatted string.

    Parameters
    ----------
    s : str
        A string in ESRI ASCII format.
    at : {'node', 'patch', 'corner', 'cell'}, optional
        Location on the grid where data are placed.
    name : str, optional
        Name of the newly-created field. If `name` is
        not provided, the grid will be created but the
        data not added to the grid as a field.
    out : RasterModelGrid, optional
        Place the data from the file onto an existing grid. If not
        provided, create a new grid.

    Returns
    -------
    RasterModelGrid
        A newly-created ``RasterModelGrid`` with, optionaly, the data added
        as a field (if `name` was provided).

    Examples
    --------
    >>> from landlab.io import esri_ascii

    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCORNER -2.0
    ... YLLCORNER 4.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... '''.lstrip()
    >>> grid = esri_ascii.loads(contents, at="cell")
    >>> grid.shape
    (3, 4)
    >>> grid.xy_of_lower_left
    (-3.0, 3.0)
    >>> grid.spacing
    (2.0, 2.0)

    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCENTER -2.0
    ... YLLCENTER 4.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... 10 11
    ... '''.lstrip()
    >>> grid = esri_ascii.loads(contents, at="cell", name="foo")
    >>> grid.shape
    (3, 4)
    >>> grid.xy_of_lower_left
    (-4.0, 2.0)
    >>> grid.spacing
    (2.0, 2.0)
    >>> grid.at_cell["foo"]
    array([10., 11.])

    >>> contents = '''
    ... NROWS 3
    ... NCOLS 4
    ... XLLCENTER 0.0
    ... YLLCENTER 0.0
    ... CELLSIZE 1.0
    ... NODATA_VALUE -9999
    ... 1 2 3 4
    ... 5 6 7 8
    ... 9 10 11 12
    ... '''.lstrip()
    >>> grid = esri_ascii.loads(contents, at="node", name="foo")
    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCENTER 1.0
    ... YLLCENTER 1.0
    ... CELLSIZE 1.0
    ... NODATA_VALUE -9999
    ... 10 20
    ... '''.lstrip()
    >>> grid = esri_ascii.loads(contents, at="cell", name="foo", out=grid)
    >>> grid.at_node["foo"].reshape(grid.shape)
    array([[ 9., 10., 11., 12.],
           [ 5.,  6.,  7.,  8.],
           [ 1.,  2.,  3.,  4.]])
    >>> grid.at_cell["foo"]
    array([10., 20.])
    """
    if name is None:
        info = parse(s, with_data=False)
    else:
        info, data = parse(s, with_data=True)
    header = _Header(**info)

    shape = np.asarray(header.shape)
    if at in ("corner", "patch"):
        shape += 1
    elif at == "cell":
        shape += 2

    shift = header.cell_size * _get_lower_left_shift(at=at, ref=header.lower_left_ref)
    xy_of_lower_left = np.asarray(header.lower_left) + shift

    if out is None:
        grid = RasterModelGrid(
            shape, xy_spacing=header.cell_size, xy_of_lower_left=xy_of_lower_left
        )
    else:
        grid = _validate_grid(
            out, shape, xy_spacing=header.cell_size, xy_of_lower_left=xy_of_lower_left
        )

    if name is not None:
        getattr(grid, f"at_{at}")[name] = data

    return grid


def _header_to_grid_args(info: dict[str, int | float], at: str):
    header = _Header(**info)

    shape = np.asarray(header.shape)
    if at in ("corner", "patch"):
        shape += 1
    elif at == "cell":
        shape += 2

    shift = header.cell_size * _get_lower_left_shift(at=at, ref=header.lower_left_ref)
    # xy_of_lower_left = np.asarray(header.lower_left) + shift

    return {
        "shape": tuple(shape),
        "xy_spacing": (header.cell_size, header.cell_size),
        "xy_of_lower_left": tuple(np.asarray(header.lower_left) + shift),
    }


class RasterGridArgs(NamedTuple):
    shape: tuple[int, int]
    xy_spacing: tuple[float, float]
    xy_of_lower_left: tuple[float, float]


@overload
def lazy_load(stream: TextIO, at: str = "node", name: None = ...) -> RasterGridArgs: ...


@overload
def lazy_load(
    stream: TextIO, at: str = "node", name: str = ...
) -> tuple[RasterGridArgs, NDArray]: ...


def lazy_load(
    stream: TextIO,
    at: str = "node",
    name: str | None = None,
) -> RasterGridArgs | tuple[RasterGridArgs, NDArray]:
    """Parse a RasterModelGrid from a ESRI ASCII format stream.

    Parameters
    ----------
    stream : file_like
        A text stream in ESRI ASCII format.
    at : {'node', 'patch', 'corner', 'cell'}, optional
        Location on the grid where data are placed.
    name : str, optional
        Name of the newly-created field. If `name` is
        not provided, the grid will be created but the
        data not added as a field.

    Returns
    -------
    RasterGridArgs
        The header metadata
    NDArray, optional
        The data as a numpy array.
    """
    return lazy_loads(stream.read(), at=at, name=name)


@overload
def lazy_loads(s: str, at: str = "node", name: None = ...) -> RasterGridArgs: ...


@overload
def lazy_loads(
    s: str, at: str = "node", name: str = ...
) -> tuple[RasterGridArgs, NDArray]: ...


def lazy_loads(
    s: str,
    at: str = "node",
    name: str | None = None,
    out: RasterModelGrid | None = None,
) -> RasterGridArgs | tuple[RasterGridArgs, NDArray]:
    """Parse a ESRI ASCII formatted string.

    Parameters
    ----------
    s : str
        A string in ESRI ASCII format.
    at : {'node', 'patch', 'corner', 'cell'}, optional
        Location on the grid where data are placed.
    name : str, optional
        Name of the newly-created field. If `name` is
        not provided, the grid will be created but the
        data not added to the grid as a field.
    out : RasterModelGrid, optional
        Place the data from the file onto an existing grid. If not
        provided, create a new grid.

    Returns
    -------
    RasterGridArgs
        The header metadata
    NDArray, optional
        The data as a numpy array.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.io import esri_ascii

    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCORNER -2.0
    ... YLLCORNER 4.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... 10 20
    ... '''.lstrip()
    >>> args, data = esri_ascii.lazy_loads(contents, at="cell", name="foo")
    >>> args
    RasterGridArgs(shape=(3, 4), xy_spacing=(2.0, 2.0), xy_of_lower_left=(-3.0, 3.0))
    >>> data
    array([10., 20.])

    >>> grid = RasterModelGrid(*args)
    >>> grid.at_cell["foo"] = data

    >>> contents = '''
    ... NROWS 3
    ... NCOLS 4
    ... XLLCORNER -3.0
    ... YLLCORNER 3.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... 1 2 3 4
    ... 5 6 7 8
    ... 9 10 11 12
    ... '''.lstrip()
    >>> args, data = esri_ascii.lazy_loads(contents, at="node", name="foo", out=grid)
    >>> args
    RasterGridArgs(shape=(3, 4), xy_spacing=(2.0, 2.0), xy_of_lower_left=(-3.0, 3.0))
    >>> data
    array([ 9., 10., 11., 12.,  5.,  6.,  7.,  8.,  1.,  2.,  3.,  4.])

    >>> grid.at_cell["foo"]
    array([10., 20.])
    >>> grid.at_node["foo"]
    array([ 9., 10., 11., 12.,  5.,  6.,  7.,  8.,  1.,  2.,  3.,  4.])
    """
    if name is None:
        info = parse(s, with_data=False)
    else:
        info, data = parse(s, with_data=True)

    args = RasterGridArgs(**_header_to_grid_args(info, at=at))

    if name is not None:
        if out is not None:
            grid = _validate_grid(out, *args)
            getattr(grid, f"at_{at}")[name] = data
        return args, data
    else:
        return args


@overload
def parse(s: str, with_data: Literal[False] = ...) -> dict[str, int | float]: ...


@overload
def parse(
    s: str, with_data: Literal[True] = ...
) -> tuple[dict[str, int | float], NDArray]: ...


def parse(
    s: str, with_data: bool = False
) -> dict[str, int | float] | tuple[dict[str, int | float], NDArray]:
    """Parse a ESRI ASCII formatted string.

    Parameters
    ----------
    s : str
        String to parse.
    with_data : bool, optional
        Return the data as a numpy array, otherwise, just
        return the header.

    Returns
    -------
    dict or tuple of (dict, ndarray)
        The header metadata and, optionally, the data.

    Examples
    --------
    >>> from landlab.io import esri_ascii

    >>> contents = '''
    ... NROWS 2
    ... NCOLS 3
    ... XLLCORNER -2.0
    ... YLLCORNER 4.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... 10 20 30
    ... 40 50 60
    ... '''.lstrip()
    >>> list(esri_ascii.parse(contents).items())
    [('nrows', 2),
     ('ncols', 3),
     ('xllcorner', -2.0),
     ('yllcorner', 4.0),
     ('cellsize', 2.0),
     ('nodata_value', -9999.0)]

    >>> info, data = esri_ascii.parse(contents, with_data=True)
    >>> data
    array([40., 50., 60., 10., 20., 30.])
    """
    lines = s.splitlines()

    start_of_data: int | None = 0
    for lineno, line in enumerate(lines):
        if not _Header.is_header_line(line):
            start_of_data = lineno
            break
    else:
        start_of_data = None

    if start_of_data is None:
        header, body = lines, []
    else:
        header, body = lines[:start_of_data], lines[start_of_data:]

    if not header:
        raise EsriAsciiError("missing header")
    if with_data and not body:
        raise EsriAsciiError("missing data")

    info = OrderedDict(_Header.parse_header_line(line) for line in header)

    validated_info = _Header(**info)

    if with_data:
        data = np.loadtxt([" ".join(body)]).reshape(validated_info.shape)

        return info, np.flipud(data).reshape((-1,))
    else:
        return info


def _validate_grid(
    grid: RasterModelGrid,
    shape: ArrayLike | None = None,
    xy_spacing: float | tuple[float, float] | None = None,
    xy_of_lower_left: tuple[float, float] | None = None,
) -> RasterModelGrid:
    if not isinstance(grid, RasterModelGrid):
        raise EsriAsciiError(
            "Not a RasterModelGrid. Only RasterModelGrids can be expressed in"
            " ESRI ASCII format."
        )
    if grid.dx != grid.dy:
        raise EsriAsciiError(f"x and y spacing must be equal ({grid.dx} != {grid.dy}).")

    if shape is not None and not np.all(np.equal(grid.shape, shape)):
        raise EsriAsciiError(f"Grid shape mismatch ({grid.shape} != {shape}).")
    if xy_spacing is not None and not np.all(np.equal(grid.spacing, xy_spacing)):
        raise EsriAsciiError(f"Grid spacing mismatch ({grid.dx} != {xy_spacing}).")
    if xy_of_lower_left is not None and not np.all(
        np.equal(grid.xy_of_lower_left, xy_of_lower_left)
    ):
        raise EsriAsciiError(
            f"Grid lower-left mismatch ({grid.xy_of_lower_left} != {xy_of_lower_left})."
        )

    return grid


def _get_lower_left_shift(at="node", ref="center"):
    """Calculate the shift from the lower left of a grid to ESRI ASCII.

    Parameters
    ----------
    at : {'node', 'patch', 'corner', 'cell'}, optional
        Grid location where data are placed.
    ref : {'center', 'corner'}, optional
        Reference location with respect to the ESRI ASCII cell.

    Returns
    -------
    float
        The unit shift between the lower left of an ESRI ASCII field
        and a ``RasterModelGrid``.
    """
    if ref not in ("center", "corner"):
        raise ValueError(
            f"Unrecognized reference location ({ref!r} not one of 'center', 'corner')"
        )

    if at == "node":
        return 0.0 if ref == "center" else 0.5
    if at == "patch":
        return -0.5 if ref == "center" else 0.0
    if at == "corner":
        return -0.5 if ref == "center" else 0.0
    if at == "cell":
        return -1.0 if ref == "center" else -0.5

    raise ValueError(
        f"Unrecognized grid location ({at!r} not one of"
        " 'cell', 'corner', 'node', 'patch')"
    )


class _Header:
    required_keys = frozenset(("ncols", "nrows", "cellsize"))
    optional_keys = frozenset(
        ("xllcorner", "xllcenter", "yllcorner", "yllcenter", "nodata_value")
    )

    def __init__(self, **kwds):
        _valid_keywords_or_raise(
            kwds, required=_Header.required_keys, optional=_Header.optional_keys
        )

        self._shape = self._validate_shape(kwds["nrows"], kwds["ncols"])
        self._cell_size = self._validate_cell_size(kwds["cellsize"])
        self._nodata_value = kwds.get("nodata_value", -9999)

        lower_left = {
            k: kwds[k]
            for k in ("xllcorner", "xllcenter", "yllcorner", "yllcenter")
            if k in kwds
        }

        self._lower_left, self._lower_left_ref = self._validate_lower_left(**lower_left)

    @staticmethod
    def parse_header_line(line: str) -> tuple[str, int | float]:
        """Parse a header line into a key/value pair."""
        try:
            key, value = line.split(maxsplit=1)
        except ValueError:
            raise BadHeaderError(
                f"header line must contain a key/value pair ({line!r})"
            ) from None
        else:
            key = key.lower()

        convert: dict[str, Callable[[str], int | float]] = {"nrows": int, "ncols": int}

        try:
            return key, convert.get(key, float)(value)
        except ValueError:
            raise BadHeaderError(
                f"unable to convert header value to a number ({line!r})"
            ) from None

    @staticmethod
    def is_header_line(line: str) -> bool:
        """Check if a string is a possible header line."""
        try:
            key, value = line.split(maxsplit=1)
        except ValueError:
            return False
        else:
            return key.lower() in (_Header.required_keys | _Header.optional_keys)

    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @staticmethod
    def _validate_shape(n_rows: int, n_cols: int) -> tuple[int, int]:
        if not isinstance(n_rows, (int, np.integer)) or n_rows <= 0:
            raise BadHeaderError(f"n_rows must be a positive integer ({n_rows!r})")
        if not isinstance(n_cols, (int, np.integer)) or n_cols <= 0:
            raise BadHeaderError(f"n_cols must be a positive integer ({n_cols!r})")

        return n_rows, n_cols

    @property
    def lower_left(self) -> tuple[float, float]:
        return self._lower_left

    @property
    def lower_left_ref(self) -> str:
        return self._lower_left_ref

    @staticmethod
    def _validate_lower_left(cell_size=0.0, **kwds) -> tuple[tuple[float, float], str]:
        if set(kwds) == {"xllcorner", "yllcorner"}:
            lower_left = kwds["xllcorner"], kwds["yllcorner"]
            ref = "corner"
        elif set(kwds) == {"xllcenter", "yllcenter"}:
            lower_left = kwds["xllcenter"], kwds["yllcenter"]
            ref = "center"
        else:
            raise BadHeaderError(
                "header must contain one, and only one, of the pairs"
                " ('xllcenter', 'yllcenter') or ('xllcorner', 'yllcorner')"
                f" (got {', '.join(repr(k) for k in sorted(set(kwds)))})"
            )
        return lower_left, ref

    @property
    def cell_size(self) -> float:
        return self._cell_size

    @staticmethod
    def _validate_cell_size(cell_size: float) -> float:
        if cell_size <= 0.0:
            raise BadHeaderError(f"cell size must be greater than zero ({cell_size})")
        return float(cell_size)

    @property
    def nodata(self) -> int | float:
        return self._nodata_value

    def __str__(self) -> str:
        lines = (
            f"NROWS {self.shape[0]}",
            f"NCOLS {self.shape[1]}",
            f"CELLSIZE {self.cell_size}",
            f"XLL{self.lower_left_ref.upper()} {self.lower_left[0]}",
            f"YLL{self.lower_left_ref.upper()} {self.lower_left[1]}",
            f"NODATA_VALUE {self.nodata}",
        )
        return "\n".join(lines) + "\n"

    def __repr__(self) -> str:
        kwds = {
            "nrows": self.shape[0],
            "ncols": self.shape[1],
            "cellsize": self.cell_size,
            f"xll{self.lower_left_ref}": self.lower_left[0],
            f"yll{self.lower_left_ref}": self.lower_left[1],
            "nodata_value": self.nodata,
        }
        args = [f"{k}={v!r}" for k, v in kwds.items()]
        return f"_Header({', '.join(args)})"
