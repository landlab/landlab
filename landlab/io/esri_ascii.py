#! /usr/bin/env python
"""Read/write data from an ESRI ASCII file into a RasterModelGrid.

ESRI ASCII functions
++++++++++++++++++++

.. autosummary::

    ~landlab.io.esri_ascii.read_asc_header
    ~landlab.io.esri_ascii.read_esri_ascii
    ~landlab.io.esri_ascii.write_esri_ascii
"""
from __future__ import annotations

import io
import os
import pathlib
import re
from collections import OrderedDict
from collections.abc import Callable
from collections.abc import Generator
from collections.abc import Iterable
from typing import Literal
from typing import NamedTuple
from typing import TextIO
from typing import overload

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab.grid.raster import RasterModelGrid
from landlab.layers.eventlayers import _valid_keywords_or_raise
from landlab.utils.add_halo import add_halo


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
    at: str = "node",
    name: str | None = None,
) -> str | None:
    """Dump a grid field to ESRII ASCII format.

    Parameters
    ----------
    grid : RasterModelGrid
        A Landlab grid.
    stream : file_like or None, optional
        A ``file_like`` object to write to. If ``None``, return
        a string containing the serialized grid.
    at : {"node", "patch", "corner", "cell"}, optional
        Where the field to be written is located on the grid.
    name : str or None, optional
        The name of the field to be written. If ``None``, only the header
        information will be written.

    Returns
    -------
    str or None
        The grid field in ESRI ASCII format or ``None`` if a ``file_like``
        object was provided.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.esri_ascii import dump

    >>> grid = RasterModelGrid((3, 4), xy_spacing=2.0)
    >>> print(dump(grid, at="node"))
    NROWS 3
    NCOLS 4
    CELLSIZE 2.0
    XLLCENTER 0.0
    YLLCENTER 0.0
    NODATA_VALUE -9999

    >>> print(dump(grid, at="cell"))
    NROWS 1
    NCOLS 2
    CELLSIZE 2.0
    XLLCENTER 2.0
    YLLCENTER 2.0
    NODATA_VALUE -9999
    """
    grid = validate_grid(grid)

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
    RasterModelGrid
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
    >>> from landlab.io.esri_ascii import loads

    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCORNER -2.0
    ... YLLCORNER 4.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... '''.lstrip()
    >>> grid = loads(contents, at="cell")
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
    >>> grid = loads(contents, at="cell", name="foo")
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
    >>> grid = loads(contents, at="node", name="foo")
    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCENTER 1.0
    ... YLLCENTER 1.0
    ... CELLSIZE 1.0
    ... NODATA_VALUE -9999
    ... 10 20
    ... '''.lstrip()
    >>> grid = loads(contents, at="cell", name="foo", out=grid)
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
        grid = validate_grid(
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
    >>> from landlab.io.esri_ascii import lazy_loads
    >>> from landlab import RasterModelGrid

    >>> contents = '''
    ... NROWS 1
    ... NCOLS 2
    ... XLLCORNER -2.0
    ... YLLCORNER 4.0
    ... CELLSIZE 2.0
    ... NODATA_VALUE -9999
    ... 10 20
    ... '''.lstrip()
    >>> args, data = lazy_loads(contents, at="cell", name="foo")
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
    >>> args, data = lazy_loads(contents, at="node", name="foo", out=grid)
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
            grid = validate_grid(out, *args)
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
    >>> from landlab.io.esri_ascii import parse

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
    >>> parse(contents)
    OrderedDict({'nrows': 2,
                 'ncols': 3,
                 'xllcorner': -2.0,
                 'yllcorner': 4.0,
                 'cellsize': 2.0,
                 'nodata_value': -9999.0})
    >>> info, data = parse(contents, with_data=True)
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


def validate_grid(
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
    if at == "node" or (at == "patch" and ref == "corner"):
        return 0.0
    elif (
        at == "corner"
        or (at == "patch" and ref == "center")
        or (at == "cell" and ref == "corner")
    ):
        return -0.5
    elif at == "cell" and ref == "center":
        return -1.0

    raise RuntimeError(f"unreachable code ({at!r}, {ref!r}")


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


_VALID_HEADER_KEYS = (
    "ncols",
    "nrows",
    "xllcorner",
    "xllcenter",
    "yllcorner",
    "yllcenter",
    "cellsize",
    "nodata_value",
)
_HEADER_KEY_REGEX_PATTERN = re.compile(r"\s*(?P<key>[a-zA-z]\w+)")
_HEADER_REGEX_PATTERN = re.compile(r"\s*(?P<key>[a-zA-Z]\w+)\s+(?P<value>[\w.+-]+)")
_HEADER_VALUE_TESTS = {
    "nrows": (int, lambda x: x > 0),
    "ncols": (int, lambda x: x > 0),
    "cellsize": (float, lambda x: x > 0),
    "xllcorner": (float, lambda x: True),
    "xllcenter": (float, lambda x: True),
    "yllcorner": (float, lambda x: True),
    "yllcenter": (float, lambda x: True),
    "nodata_value": (float, lambda x: True),
}


class Error(Exception):
    """Base class for errors in this module."""

    pass


class BadHeaderLineError(Error):
    """Raise this error for a bad header is line."""

    def __init__(self, line: str) -> None:
        self._line = line

    def __str__(self) -> str:
        return self._line


class MissingRequiredKeyError(Error):
    """Raise this error when a header is missing a required key."""

    def __init__(self, key: str) -> None:
        self._key = key

    def __str__(self) -> str:
        return self._key


class KeyTypeError(Error):
    """Raise this error when a header's key value is of the wrong type."""

    def __init__(self, key: str, expected_type: type | str):
        self._key = key
        self._type = str(expected_type)

    def __str__(self) -> str:
        return f"Unable to convert {self._key} to {self._type}"


class KeyValueError(Error):
    """Raise this error when a header's key value has a bad value."""

    def __init__(self, key: str, message: str) -> None:
        self._key = key
        self._msg = message

    def __str__(self) -> str:
        return f"{self._key}: {self._msg}"  # this line not yet tested


class DataSizeError(Error):
    """Raise this error if the size of data does not match the header."""

    def __init__(self, size: int, expected_size: int) -> None:
        self._actual = size
        self._expected = expected_size

    def __str__(self) -> str:
        return f"{self._actual} != {self._expected}"


class MismatchGridDataSizeError(Error):
    """Raise this error if the data size does not match the grid size."""

    def __init__(self, size: int, expected_size: int) -> None:
        self._actual = size
        self._expected = expected_size

    def __str__(self) -> str:
        return f"(data size) {self._actual} != {self._expected} (grid size)"


class MismatchGridXYSpacing(Error):
    """Raise this error if the file cell size does not match the grid dx."""

    def __init__(
        self, dx: tuple[float, float], expected_dx: tuple[float, float]
    ) -> None:
        self._actual = dx
        self._expected = expected_dx

    def __str__(self):
        return f"(data dx) {self._actual} != {self._expected} (grid dx)"


class MismatchGridXYLowerLeft(Error):
    """Raise this error if the file lower left does not match the grid."""

    def __init__(
        self, llc: tuple[float, float], expected_llc: tuple[float, float]
    ) -> None:
        self._actual = llc
        self._expected = expected_llc

    def __str__(self) -> str:
        return f"(data lower-left) {self._actual} != {self._expected} (grid lower-left)"


def _parse_header_key_value(line: str) -> tuple[str, str] | None:
    """Parse a header line into a key-value pair.

    Parameters
    ----------
    line : str
        Header line.

    Returns
    -------
    (str, str)
        Header key-value pair

    Raises
    ------
    BadHeaderLineError
        There is something wrong with the header line.
    """
    match = _HEADER_KEY_REGEX_PATTERN.match(line)
    if match is None:
        return None
        # raise BadHeaderLineError(line)

    match = _HEADER_REGEX_PATTERN.match(line)
    if match is None:
        raise BadHeaderLineError(line)

    (key, value) = (match.group("key").lower(), match.group("value"))

    if key in _VALID_HEADER_KEYS:
        return (key, value)
    else:
        raise BadHeaderLineError(line)


def _header_lines(asc_file: TextIO) -> Generator[tuple[str, str], None, None]:
    """Iterate over header lines for a ESRI ASCII file.

    Parameters
    ----------
    asc_file : file_like
        File-like object for an ESRI ASCII file.

    Yields
    ------
    str
        Header line.
    """
    pos = asc_file.tell()
    line = asc_file.readline()
    while len(line) > 0:
        if len(line.strip()) > 0:
            item = _parse_header_key_value(line)
            if item:
                yield item
            else:
                asc_file.seek(pos, 0)
                break
        pos = asc_file.tell()
        line = asc_file.readline()


def _header_is_valid(header: dict[str, str]) -> dict[str, int | float]:
    """Check if the ESRI ASCII header is valid.

    Parameters
    ----------
    header : dict
        Header as key-values pairs.

    Raises
    ------
    MissingRequiredKeyError
        The header is missing a required key.
    KeyTypeError
        The header has the key but its values is of the wrong type.
    """
    header_keys = set(header)
    required_keys = {"ncols", "nrows", "cellsize"}

    if not required_keys.issubset(header_keys):
        raise MissingRequiredKeyError(", ".join(required_keys - header_keys))

    for keys in [("xllcenter", "xllcorner"), ("yllcenter", "yllcorner")]:
        if len(set(keys) & header_keys) != 1:
            raise MissingRequiredKeyError("|".join(keys))

    validated_header = {}
    for key, requires in _HEADER_VALUE_TESTS.items():
        to_type, is_valid = requires

        if key not in header:
            continue

        try:
            validated_header[key] = to_type(header[key])
        except ValueError as exc:
            raise KeyTypeError(key, to_type) from exc

        if not is_valid(validated_header[key]):
            raise KeyValueError(key, "Bad value")

    return validated_header


def read_asc_header(asc_file: TextIO) -> dict[str, int | float]:
    """Read header information from an ESRI ASCII raster file.

    The header contains the following variables,

    * ``ncols``: Number of cell columns
    * ``nrows``: Number of cell rows
    * ``xllcenter`` or ``xllcorner``: X (column) coordinate of lower-left
        coordinate of grid (by center or lower-left corner of the cell)
    * ``yllcenter``, ``yllcorner``: Y (row) coordinate of lower-left
        coordinate of grid (by center or lower-left corner of the cell)
    * ``cellsize``: Grid spacing between rows and columns
    * ``nodata_value``: No-data value (optional)

    Parameters
    ----------
    asc_file : file_like
        File-like object from which to read header.

    Returns
    -------
    dict
        Header as key-value pairs.

    Raises
    ------
    :class:`~landlab.io.esri_ascii.MissingRequiredKeyError`
        The header is missing a required key.
    :class:`~landlab.io.esri_ascii.KeyTypeError`
        The header has the key but its values is of the wrong type.

    Examples
    --------
    >>> from io import StringIO
    >>> from landlab.io.esri_ascii import read_asc_header

    >>> contents = '''
    ... nrows 100
    ... ncols 200
    ... cellsize 1.5
    ... xllcenter 0.5
    ... yllcenter -0.5
    ... '''

    >>> hdr = read_asc_header(StringIO(contents))
    >>> hdr["nrows"], hdr["ncols"]
    (100, 200)
    >>> hdr["cellsize"]
    1.5
    >>> hdr["xllcenter"], hdr["yllcenter"]
    (0.5, -0.5)

    :class:`~landlab.io.esri_ascii.MissingRequiredKeyError` is raised if the
    header does not contain all of the necessary keys.

    >>> contents = '''
    ... ncols 200
    ... cellsize 1.5
    ... xllcenter 0.5
    ... yllcenter -0.5
    ... '''
    >>> read_asc_header(StringIO(contents))
    Traceback (most recent call last):
    MissingRequiredKeyError: nrows

    :class:`~landlab.io.esri_ascii.KeyTypeError` is raised if a value is of
    the wrong type. For instance, *nrows* and *ncols* must be ``int``.

    >>> contents = '''
    ... nrows 100.5
    ... ncols 200
    ... cellsize 1.5
    ... xllcenter 0.5
    ... yllcenter -0.5
    ... '''
    >>> read_asc_header(StringIO(contents))
    Traceback (most recent call last):
    KeyTypeError: Unable to convert nrows to <type 'int'>
    """
    header = {}
    for key, value in _header_lines(asc_file):
        header[key] = value

    return _header_is_valid(header)


def _read_asc_data(asc_file: TextIO) -> NDArray:
    """Read gridded data from an ESRI ASCII data file.

    Parameters
    ----------
    asc_file : file-like
        File-like object of the data file pointing to the start of the data.

    .. note::
        First row of the data is at the top of the raster grid, the second
        row is the second from the top, and so on.
    """
    return np.loadtxt(asc_file)


def read_esri_ascii(
    asc_file: str | TextIO,
    grid: RasterModelGrid | None = None,
    reshape: bool = False,
    name: str | None = None,
    halo: int = 0,
) -> tuple[RasterModelGrid, NDArray]:
    """Read :py:class:`~.RasterModelGrid` from an ESRI ASCII file.

    Read data from *asc_file*, an `ESRI ASCII file`_, into a
    :class:`~.RasterModelGrid`.  *asc_file* is either the name of
    the data file or is a file-like object.

    The grid and data read from the file are returned as a tuple
    (*grid*, *data*) where *grid* is an instance of
    :py:class:`~.RasterModelGrid` and *data* is a numpy
    array of doubles with that has been reshaped to have the number of rows
    and columns given in the header.

    .. _ESRI ASCII file: https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm

    Parameters
    ----------
    asc_file : str of file-like
        Data file to read.
    reshape : boolean, optional
        Reshape the returned array, otherwise return a flattened array.
    name : str, optional
        Add data to the grid as a named field.
    grid : *grid* , optional
        Adds data to an existing *grid* instead of creating a new one.
    halo : integer, optional
        Adds outer border of depth halo to the *grid*.

    Returns
    -------
    (grid, data) : tuple
        A newly-created :class:`~.RasterModelGrid` and the associated node data.

    Raises
    ------
    :class:`~landlab.io.esri_ascii.DataSizeError`
        Data are not the same size as indicated by the header file.
    :class:`~landlab.io.esri_ascii.MismatchGridDataSizeError`
        If a grid is passed, and the size of the grid does not agree with the
        size of the data.
    :class:`~landlab.io.esri_ascii.MismatchGridXYSpacing`
        If a grid is passed, and the *cellsize* listed in the heading does not
        match the node spacing of the grid.
    :class:`~landlab.io.esri_ascii.MismatchGridXYLowerLeft`
        If a grid is passed and the *xllcorner* and *yllcorner* do not match that
        of the grid.

    Examples
    --------

    >>> from landlab.io import read_esri_ascii
    >>> from io import StringIO

    >>> contents = '''
    ... ncols         3
    ... nrows         4
    ... xllcorner     1.
    ... yllcorner     2.
    ... cellsize      10.
    ... NODATA_value  -1
    ... 0. 1. 2.
    ... 3. 4. 5.
    ... 6. 7. 8.
    ... 9. 10. 11.
    ... '''

    >>> (grid, data) = read_esri_ascii(StringIO(contents))

    The returned grid is a :class:`~.RasterModelGrid` with 4 rows and 3 columns.

    >>> grid
    RasterModelGrid((4, 3), xy_spacing=(10.0, 10.0), xy_of_lower_left=(1.0, 2.0))

    Note that the first row of values is the bottom-most of the data file.

    >>> data.reshape(grid.shape)
    array([[  9.,  10.,  11.],
           [  6.,   7.,   8.],
           [  3.,   4.,   5.],
           [  0.,   1.,   2.]])

    >>> (grid, data) = read_esri_ascii(StringIO(contents), halo=1)

    Because of the halo, the returned grid now has two more rows and columns than before.

    >>> grid
    RasterModelGrid((6, 5), xy_spacing=(10.0, 10.0), xy_of_lower_left=(-9.0, -8.0))
    >>> data.reshape(grid.shape)
    array([[-1.,  -1.,  -1.,  -1.,  -1.],
           [-1.,   9.,  10.,  11.,  -1.],
           [-1.,   6.,   7.,   8.,  -1.],
           [-1.,   3.,   4.,   5.,  -1.],
           [-1.,   0.,   1.,   2.,  -1.],
           [-1.,  -1.,  -1.,  -1.,  -1.]])
    """  # noqa: B950
    from ..grid import RasterModelGrid

    if halo < 0:
        raise ValueError("negative halo")

    # if the asc_file is provided as a string, open it and pass the pointer to
    # _read_asc_header, and _read_asc_data
    if isinstance(asc_file, (str, pathlib.Path)):
        with open(asc_file) as f:
            header = read_asc_header(f)
            data = _read_asc_data(f)

    # otherwise, pass asc_file directly.
    else:
        header = read_asc_header(asc_file)
        data = _read_asc_data(asc_file)

    shape = (int(header["nrows"]) + 2 * halo, int(header["ncols"]) + 2 * halo)
    nodata_value = header.get("nodata_value", -9999.0)
    if data.size != (shape[0] - 2 * halo) * (shape[1] - 2 * halo):
        raise DataSizeError(shape[0] * shape[1], data.size)

    xy_spacing = (header["cellsize"], header["cellsize"])
    xy_of_lower_left = (
        header["xllcorner"] - halo * header["cellsize"],
        header["yllcorner"] - halo * header["cellsize"],
    )

    data = np.flipud(data)

    if halo > 0:
        data = add_halo(
            data.reshape((int(header["nrows"]), int(header["ncols"]))),
            halo=halo,
            halo_value=nodata_value,
        ).reshape((-1,))

    if not reshape:
        data = data.flatten()

    if grid is not None:
        if (grid.number_of_node_rows != shape[0]) or (
            grid.number_of_node_columns != shape[1]
        ):
            raise MismatchGridDataSizeError(
                shape[0] * shape[1],
                grid.number_of_node_rows * grid.number_of_node_columns,
            )
        if (grid.dx, grid.dy) != xy_spacing:
            raise MismatchGridXYSpacing((grid.dx, grid.dy), xy_spacing)

        if grid.xy_of_lower_left != xy_of_lower_left:
            raise MismatchGridXYLowerLeft(grid.xy_of_lower_left, xy_of_lower_left)

    if grid is None:
        grid = RasterModelGrid(
            shape, xy_spacing=xy_spacing, xy_of_lower_left=xy_of_lower_left
        )
    if name:
        grid.add_field(name, data, at="node")

    return (grid, data)


def write_esri_ascii(
    path: str,
    fields: RasterModelGrid,
    names: Iterable[str] | None = None,
    clobber: bool = False,
) -> list[str]:
    """Write landlab fields to ESRI ASCII.

    Write the data and grid information for *fields* to *path* in the ESRI
    ASCII format.

    Parameters
    ----------
    path : str
        Path to output file.
    fields : field-like
        Landlab field object that holds a grid and associated values.
    names : iterable of str, optional
        Names of the fields to include in the output file. If not provided,
        write all fields.
    clobber : boolean
        If *path* exists, clobber the existing file, otherwise raise an
        exception.

    Examples
    --------
    >>> import numpy as np
    >>> import os
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.esri_ascii import write_esri_ascii

    >>> grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    >>> grid.at_node["air__temperature"] = np.arange(20.0)
    >>> files = write_esri_ascii("test.asc", grid)  # doctest: +SKIP
    >>> [os.path.basename(name) for name in sorted(files)]  # doctest: +SKIP
    ['test.asc']

    >>> _ = grid.add_field("land_surface__elevation", np.arange(20.0), at="node")
    >>> grid.at_node["land_surface__elevation"] = np.arange(20.0)
    >>> files = write_esri_ascii("test.asc", grid)  # doctest: +SKIP
    >>> [os.path.basename(name) for name in sorted(files)]  # doctest: +SKIP
    ['test_air__temperature.asc', 'test_land_surface__elevation.asc']
    """
    if os.path.exists(path) and not clobber:
        raise ValueError("file exists")

    if isinstance(names, (str, pathlib.Path)):
        names = [names]

    names = tuple(names or fields.at_node)

    if len(names) == 1:
        paths = [path]
    elif len(names) > 1:
        (base, ext) = os.path.splitext(path)
        paths = [base + "_" + name + ext for name in names]
    else:
        raise ValueError("no node fields to write")

    bad_names = set(names) - set(fields.at_node.keys())
    if len(bad_names) > 0:
        raise ValueError("unknown field name(s): %s" % ",".join(bad_names))

    header = {
        "ncols": fields.number_of_node_columns,
        "nrows": fields.number_of_node_rows,
        "xllcorner": fields.node_x[0],
        "yllcorner": fields.node_y[0],
        "cellsize": fields.dx,
    }

    for path, name in zip(paths, names):
        header_lines = [f"{key} {str(val)}" for key, val in list(header.items())]
        data = fields.at_node[name].reshape(header["nrows"], header["ncols"])
        np.savetxt(
            path, np.flipud(data), header=os.linesep.join(header_lines), comments=""
        )

    return paths
