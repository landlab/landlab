import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.field.errors import FieldError
from landlab.io import esri_ascii


@pytest.mark.parametrize("with_data", (True, False))
def test_parse_no_head(with_data):
    """Check if the file doesn't contain a header."""
    contents = "10 20 30"
    with pytest.raises(esri_ascii.EsriAsciiError):
        esri_ascii.parse(contents, with_data=with_data)


def test_parse_no_body():
    """Check if the file doesn't contain data."""
    contents = """\
NROWS 4
NCOLS 3
CELLSIZE 10.0
XLLCENTER 1.0
YLLCENTER 2.0
NODATA_VALUE -9999
"""

    expected = [
        ("nrows", 4),
        ("ncols", 3),
        ("cellsize", 10.0),
        ("xllcenter", 1.0),
        ("yllcenter", 2.0),
        ("nodata_value", -9999),
    ]

    header = esri_ascii.parse(contents, with_data=False)
    assert list(header.items()) == expected

    with pytest.raises(esri_ascii.EsriAsciiError):
        esri_ascii.parse(contents, with_data=True)


@pytest.mark.parametrize(
    "body",
    (
        "1 2 3\n4 5 6\n7 8 9\n10 11 12",
        "1 2 3 4 5 6 7 8 9 10 11 12",
        "\n\n1 2 3 4 5 6 7 8 9 10 11 12",
        "1 2 3 4 5 6 7 8 9 10 11 12\n\n",
        "1 2\n\n3 4 5 6 7 8 9 10 11 12",
        "1 2\t\n3\t4\n\t5 6 7 8 9 10 11 12",
    ),
)
def test_parse(body):
    """Check whitespace in the data."""
    contents = f"""\
NROWS 4
NCOLS 3
CELLSIZE 10.0
XLLCENTER 1.0
YLLCENTER 2.0
NODATA_VALUE -9999
{body}
"""
    expected = [
        ("nrows", 4),
        ("ncols", 3),
        ("cellsize", 10.0),
        ("xllcenter", 1.0),
        ("yllcenter", 2.0),
        ("nodata_value", -9999),
    ]

    header = esri_ascii.parse(contents, with_data=False)
    assert list(header.items()) == expected

    header, data = esri_ascii.parse(contents, with_data=True)
    assert list(header.items()) == expected
    assert_array_equal(data, [10, 11, 12, 7, 8, 9, 4, 5, 6, 1, 2, 3])


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
    with pytest.raises(esri_ascii.BadHeaderError):
        esri_ascii.loads(contents.format(**header))


@pytest.mark.parametrize("nodata_value", (0, 1, -999))
def test_dump_with_nodata_value(nodata_value):
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    actual = esri_ascii.dump(grid, nodata_value=nodata_value)

    assert (
        actual
        == f"""\
NROWS 4
NCOLS 3
CELLSIZE 10.0
XLLCENTER 1.0
YLLCENTER 2.0
NODATA_VALUE {nodata_value}
"""
    )


def test_dump_to_string_no_data():
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    actual = esri_ascii.dump(grid)

    assert (
        actual
        == """\
NROWS 4
NCOLS 3
CELLSIZE 10.0
XLLCENTER 1.0
YLLCENTER 2.0
NODATA_VALUE -9999
"""
    )


def test_dump_to_file_no_data(tmpdir):
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))

    expected = """\
NROWS 4
NCOLS 3
CELLSIZE 10.0
XLLCENTER 1.0
YLLCENTER 2.0
NODATA_VALUE -9999
"""
    with tmpdir.as_cwd():
        with open("foo.asc", "w") as fp:
            assert esri_ascii.dump(grid, stream=fp) is None
        with open("foo.asc") as fp:
            actual = fp.read()
        assert actual == expected


def test_dump_to_string_data_at_node():
    grid = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(1.0, 2.0))
    grid.at_node["foo"] = np.arange(12, dtype=int)
    actual = esri_ascii.dump(grid, at="node", name="foo")

    assert (
        actual
        == """\
NROWS 4
NCOLS 3
CELLSIZE 10.0
XLLCENTER 1.0
YLLCENTER 2.0
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
    actual = esri_ascii.dump(grid, at="cell", name="foo")

    assert (
        actual
        == """\
NROWS 2
NCOLS 1
CELLSIZE 10.0
XLLCENTER 11.0
YLLCENTER 12.0
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

    actual = esri_ascii.loads(
        esri_ascii.dump(grid, at=at, name="foo"), at=at, name="foo"
    )

    assert actual.shape == grid.shape
    assert actual.spacing == grid.spacing
    assert actual.xy_of_lower_left == grid.xy_of_lower_left
    assert_array_equal(
        getattr(actual, f"at_{at}")["foo"], getattr(grid, f"at_{at}")["foo"]
    )


@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
@pytest.mark.parametrize("name", ("foo", None))
def test_load(tmpdir, at, name):
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
    with tmpdir.as_cwd():
        with open("foo.asc", "w") as fp:
            print(contents, file=fp)

        with open("foo.asc") as fp:
            grid = esri_ascii.load(fp, at=at, name=name)

    if name:
        assert name in getattr(grid, f"at_{at}")
        assert getattr(grid, f"at_{at}")[name].size == 3 * 4
    else:
        assert len(getattr(grid, f"at_{at}")) == 0
    assert grid.number_of_elements(at) == 3 * 4
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
    grid = esri_ascii.loads(contents, at=at)

    assert grid.number_of_elements(at) == 3 * 4
    assert grid.dx == 10.0
    assert grid.dy == 10.0
    assert len(getattr(grid, f"at_{at}")) == 0


@pytest.mark.parametrize(
    "at,lower_left",
    (
        ("node", "xllcorner -4.0\nyllcorner -3.0"),
        ("node", "xllcenter 1.0\nyllcenter 2.0"),
        ("patch", "xllcorner 1.0\nyllcorner 2.0"),
        ("patch", "xllcenter 6.0\nyllcenter 7.0"),
        ("corner", "xllcorner 1.0\nyllcorner 2.0"),
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
    grid = esri_ascii.loads(contents + lower_left, at=at)
    assert grid.xy_of_lower_left == (1.0, 2.0)


@pytest.mark.parametrize("ref", ("center", "corner"))
@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
def test_file_round_trip(ref, at):
    cellsize = 4.0
    grid = esri_ascii.loads(
        f"""\
NROWS 3
NCOLS 4
CELLSIZE {cellsize}
XLL{ref.upper()} 5.0
YLL{ref.upper()} -4.0
NODATA_VALUE -9999
""",
        at=at,
    )
    assert grid.number_of_elements(at) == 12
    actual = esri_ascii.dump(grid, at=at)

    expected = f"""\
NROWS 3
NCOLS 4
CELLSIZE {cellsize}
XLLCENTER {5.0 if ref == "center" else 5.0 + cellsize / 2.0}
YLLCENTER {-4.0 if ref == "center" else -4.0 + cellsize / 2.0}
NODATA_VALUE -9999
"""
    assert actual == expected


@pytest.mark.parametrize("at", ("node", "patch", "corner", "cell"))
def test_dump_missing_field(at):
    grid = RasterModelGrid((4, 6))
    grid.add_ones("foo", at=at)

    with pytest.raises(FieldError):
        esri_ascii.dump(grid, at=at, name="bar")


def test_dump_not_raster():
    grid = HexModelGrid((3, 4))

    with pytest.raises(esri_ascii.EsriAsciiError):
        esri_ascii.dump(grid)


def test_dump_unequal_spacing():
    grid = RasterModelGrid((4, 6), xy_spacing=(3, 4))

    with pytest.raises(esri_ascii.EsriAsciiError):
        esri_ascii.dump(grid)


def test_reference_shift_with_bad_args():
    with pytest.raises(ValueError) as exc_info:
        esri_ascii._get_lower_left_shift(at="foo", ref="center")
    assert str(exc_info.value).startswith("Unrecognized grid location")

    with pytest.raises(ValueError) as exc_info:
        esri_ascii._get_lower_left_shift(at="node", ref="foo")
    assert str(exc_info.value).startswith("Unrecognized reference")
