#! /usr/bin/env python
"""
Unit tests for landlab.io.esri_ascii module.
"""
from io import StringIO

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.field.errors import FieldError
from landlab.io import BadHeaderLineError
from landlab.io import DataSizeError
from landlab.io import KeyTypeError
from landlab.io import KeyValueError
from landlab.io import MismatchGridDataSizeError
from landlab.io import MismatchGridXYLowerLeft
from landlab.io import MismatchGridXYSpacing
from landlab.io import MissingRequiredKeyError
from landlab.io import read_asc_header
from landlab.io import read_esri_ascii
from landlab.io.esri_ascii import BadHeaderError
from landlab.io.esri_ascii import EsriAsciiError
from landlab.io.esri_ascii import dump
from landlab.io.esri_ascii import load
from landlab.io.esri_ascii import loads


@pytest.mark.parametrize(
    "params",
    (
        {"nrows": 0},
        {"nrows": 3.2},
        {"nrows": -2},
        {"ncols": 0},
        {"ncols": 3.1},
        {"ncols": -2},
        {"cellsize": 0},
        {"cellsize": -3},
        {"cellsize": "foobar"},
        {"xll": "xllcorner", "yll": "yllcenter"},
        {"xll": "xllcenter", "yll": "yllcorner"},
        {"xll": "xllcenter", "yll": "xllcenter"},
        {"xll": "yllcorner", "yll": "yllcorner"},
    ),
)
def test_loads_bad_header(params):
    header = {
        "nrows": 3,
        "ncols": 4,
        "xll": "xllcenter",
        "yll": "yllcenter",
        "cellsize": 10.0,
    }
    header.update(params)
    contents = """\
ncols {ncols}
nrows {nrows}
{xll} 1.0
{yll} 2.0
cellsize {cellsize}
"""
    with pytest.raises(BadHeaderError):
        loads(contents.format(**header))


def test_dump_to_string_no_data():
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    actual = dump(grid)

    assert (
        actual
        == """\
NROWS 4
NCOLS 3
XLLCENTER 1.0
YLLCENTER 2.0
CELLSIZE 10.0
NODATA_VALUE -9999
"""
    )


def test_dump_to_file_no_data(tmpdir):
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))

    expected = """\
NROWS 4
NCOLS 3
XLLCENTER 1.0
YLLCENTER 2.0
CELLSIZE 10.0
NODATA_VALUE -9999
"""
    with tmpdir.as_cwd():
        with open("foo.asc", "w") as fp:
            assert dump(grid, stream=fp) is None
        with open("foo.asc") as fp:
            actual = fp.read()
        assert actual == expected


def test_dump_to_string_data_at_node():
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    grid.at_node["foo"] = np.arange(12, dtype=int)
    actual = dump(grid, at="node", name="foo")

    assert (
        actual
        == """\
NROWS 4
NCOLS 3
XLLCENTER 1.0
YLLCENTER 2.0
CELLSIZE 10.0
NODATA_VALUE -9999
9 10 11
6 7 8
3 4 5
0 1 2
"""
    )


def test_dump_to_string_data_at_cell():
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    grid.at_cell["foo"] = np.arange(2, dtype=int)
    actual = dump(grid, at="cell", name="foo")

    assert (
        actual
        == """\
NROWS 2
NCOLS 1
XLLCENTER 11.0
YLLCENTER 12.0
CELLSIZE 10.0
NODATA_VALUE -9999
1
0
"""
    )


@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
def test_dump_round_trip(at):
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    getattr(grid, f"at_{at}")["foo"] = np.arange(
        grid.number_of_elements(at), dtype=float
    )

    actual = loads(dump(grid, at=at, name="foo"), at=at, name="foo")

    assert actual.shape == grid.shape
    assert actual.spacing == grid.spacing
    assert actual.xy_of_lower_left == grid.xy_of_lower_left
    assert_array_equal(
        getattr(actual, f"at_{at}")["foo"], getattr(grid, f"at_{at}")["foo"]
    )


@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
@pytest.mark.parametrize("name", ("foo", None))
def test_load_hugo(datadir, at, name):
    with open(datadir / "hugo_site.asc") as fp:
        grid = load(fp, at=at, name=name)

    if name:
        assert name in getattr(grid, f"at_{at}")
        assert getattr(grid, f"at_{at}")[name].size == 55 * 76
    else:
        assert len(getattr(grid, f"at_{at}")) == 0
    assert grid.number_of_elements(at) == 55 * 76
    assert grid.dx == 10.0
    assert grid.dy == 10.0


@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
def test_loads(at):
    contents = """\
ncols 3
nrows 4
xllcorner 1.0
yllcorner 2.0
cellsize 10.
0. 1. 2.
3. 4. 5.
6. 7. 8.
9. 10. 11.
"""
    grid = loads(contents, at=at)

    assert grid.number_of_elements(at) == 3 * 4
    assert grid.dx == 10.0
    assert grid.dy == 10.0
    assert len(getattr(grid, f"at_{at}")) == 0


@pytest.mark.parametrize(
    "at,lower_left",
    (
        ("node", "xllcorner 1.0\nyllcorner 2.0"),
        ("node", "xllcenter 1.0\nyllcenter 2.0"),
        ("patch", "xllcorner 1.0\nyllcorner 2.0"),
        ("patch", "xllcenter 6.0\nyllcenter 7.0"),
        ("corner", "xllcorner 6.0\nyllcorner 7.0"),
        ("corner", "xllcenter 6.0\nyllcenter 7.0"),
        ("cell", "xllcorner 6.0\nyllcorner 7.0"),
        ("cell", "xllcenter 11.0\nyllcenter 12.0"),
    ),
)
def test_loads_lower_left(at, lower_left):
    contents = """\
ncols 3
nrows 4
cellsize 10.0
"""
    grid = loads(contents + lower_left, at=at)
    assert grid.xy_of_lower_left == (1.0, 2.0)


@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
def test_dump_missing_field(at):
    grid = RasterModelGrid((4, 6))
    grid.add_ones("foo", at=at)

    with pytest.raises(FieldError):
        dump(grid, at=at, name="bar")


def test_dump_not_raster():
    grid = HexModelGrid((3, 4))

    with pytest.raises(EsriAsciiError):
        dump(grid)


def test_dump_unequal_spacing():
    grid = RasterModelGrid((4, 6), xy_spacing=(3, 4))

    with pytest.raises(EsriAsciiError):
        dump(grid)


def test_hugo_read_file_name(datadir):
    (grid, field) = read_esri_ascii(datadir / "hugo_site.asc")

    assert isinstance(grid, RasterModelGrid)

    assert field.size == 55 * 76
    assert field.shape == (55 * 76,)
    assert (grid.x_of_node.min(), grid.y_of_node.min()) == (0.0, 0.0)


def test_hugo_read_file_like(datadir):
    with open(datadir / "hugo_site.asc") as asc_file:
        (grid, field) = read_esri_ascii(asc_file)

    assert isinstance(grid, RasterModelGrid)

    assert field.size == 55 * 76
    assert field.shape == (55 * 76,)


def test_hugo_reshape(datadir):
    with open(datadir / "hugo_site.asc") as asc_file:
        (grid, field) = read_esri_ascii(asc_file, reshape=True)

    assert isinstance(grid, RasterModelGrid)

    assert field.shape == (55, 76)


def test_4x3_read_file_name(datadir):
    (grid, field) = read_esri_ascii(datadir / "4_x_3.asc")

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert (grid.x_of_node.min(), grid.y_of_node.min()) == (1.0, 2.0)
    assert_array_equal(
        field, np.array([9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0])
    )


def test_4x3_read_file_like(datadir):
    with open(datadir / "4_x_3.asc") as asc_file:
        (grid, field) = read_esri_ascii(asc_file)

    assert isinstance(grid, RasterModelGrid)

    assert_array_equal(
        field, np.array([9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0])
    )


def test_4x3_shape_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4.
5. 6. 7. 8.
9. 10. 11. 12.
        """
    )
    (grid, field) = read_esri_ascii(asc_file)
    assert field.size == 12

    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    (grid, field) = read_esri_ascii(asc_file)
    assert field.size == 12


def test_4x3_size_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10.
        """
    )
    with pytest.raises(DataSizeError):
        read_esri_ascii(asc_file)


def test_4x3_size_mismatch_with_halo():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10.
        """
    )
    with pytest.raises(DataSizeError):
        read_esri_ascii(asc_file, halo=1)


def test_grid_data_size_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    rmg = RasterModelGrid((10, 10), xy_spacing=10.0)
    with pytest.raises(MismatchGridDataSizeError):
        read_esri_ascii(asc_file, grid=rmg)


def test_grid_dx_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    rmg = RasterModelGrid((4, 3), xy_spacing=15.0)
    with pytest.raises(MismatchGridXYSpacing):
        read_esri_ascii(asc_file, grid=rmg)


def test_grid_lower_left_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    rmg = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(10, 15))
    with pytest.raises(MismatchGridXYLowerLeft):
        read_esri_ascii(asc_file, grid=rmg)


def test_header_missing_required_key():
    asc_file = StringIO(
        """
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    with pytest.raises(MissingRequiredKeyError):
        read_asc_header(asc_file)


def test_header_unknown_key():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
invalid_key   1
        """
    )
    with pytest.raises(BadHeaderLineError):
        read_asc_header(asc_file)


def test_header_missing_value():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize
NODATA_value  -9999
invalid_key   1
        """
    )
    with pytest.raises(BadHeaderLineError):
        read_asc_header(asc_file)


def test_header_bad_values():
    asc_file = StringIO(
        """
nrows         -4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    with pytest.raises(KeyValueError):
        read_asc_header(asc_file)


def test_header_missing_mutex_key():
    asc_file = StringIO(
        """
ncols         3
nrows         4
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    with pytest.raises(MissingRequiredKeyError):
        read_asc_header(asc_file)


def test_header_mutex_key():
    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    header = read_asc_header(asc_file)
    assert header["xllcenter"] == 1.0
    with pytest.raises(KeyError):
        header["xllcorner"]

    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    header = read_asc_header(asc_file)
    assert header["xllcorner"] == 1.0
    with pytest.raises(KeyError):
        header["xllcenter"]


def test_header_missing_optional():
    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
        """
    )
    header = read_asc_header(asc_file)
    with pytest.raises(KeyError):
        header["nodata_value"]


def test_header_case_insensitive():
    asc_file = StringIO(
        """
nCoLs         3
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
        """
    )
    header = read_asc_header(asc_file)
    for key in ["ncols", "nrows", "xllcenter", "yllcorner", "cellsize", "nodata_value"]:
        assert key in header


def test_header_wrong_type():
    asc_file = StringIO(
        """
nCoLs         3.5
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
        """
    )
    with pytest.raises(KeyTypeError):
        read_asc_header(asc_file)


def test_name_keyword(datadir):
    (grid, field) = read_esri_ascii(datadir / "4_x_3.asc", name="air__temperature")

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert_array_equal(
        field, np.array([9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0])
    )
    assert_array_almost_equal(grid.at_node["air__temperature"], field)
    assert grid.at_node["air__temperature"] is field


def test_halo_keyword(datadir):
    (grid, field) = read_esri_ascii(datadir / "4_x_3.asc", halo=1)

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert_array_equal(
        field,
        np.array(
            [
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                9.0,
                10.0,
                11.0,
                -9999.0,
                -9999.0,
                6.0,
                7.0,
                8.0,
                -9999.0,
                -9999.0,
                3.0,
                4.0,
                5.0,
                -9999.0,
                -9999.0,
                0.0,
                1.0,
                2.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
            ]
        ),
    )


def test_halo_keyword_no_nodata_value(datadir):
    (grid, field) = read_esri_ascii(datadir / "4_x_3_no_nodata_value.asc", halo=1)

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert_array_equal(
        field,
        np.array(
            [
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                9.0,
                10.0,
                11.0,
                -9999.0,
                -9999.0,
                6.0,
                7.0,
                8.0,
                -9999.0,
                -9999.0,
                3.0,
                4.0,
                5.0,
                -9999.0,
                -9999.0,
                0.0,
                1.0,
                2.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
            ]
        ),
    )
