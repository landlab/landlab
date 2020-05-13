import pathlib

import pytest

from landlab import ExampleData


def test_fetch(tmpdir):
    data = ExampleData("io/shapefile")
    expected = set(data)
    with tmpdir.as_cwd():
        data.fetch()

        fetched = set(p.name for p in pathlib.Path(".").iterdir())
        assert expected == fetched


def test_fetch_with_case(tmpdir):
    data = ExampleData("io/shapefile", case="methow")
    expected = ["methow"]
    with tmpdir.as_cwd():
        data.fetch()

        fetched = set(p.name for p in pathlib.Path(".").iterdir())
        assert expected == fetched


def test_fetch_no_overwrite(tmpdir):
    data = ExampleData("io/shapefile")
    with tmpdir.as_cwd():
        pathlib.Path(list(data).pop()).touch()
        with pytest.raises(FileExistsError):
            data.fetch()
