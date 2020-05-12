import pathlib

import pytest

from landlab import copy_example_data


def test_copy(tmpdir):
    with tmpdir.as_cwd():
        copied = copy_example_data("io/shapefile")
        assert set(copied) == set(p.name for p in pathlib.Path(".").iterdir())


def test_copy_with_case(tmpdir):
    with tmpdir.as_cwd():
        copied = copy_example_data("io/shapefile", case="methow")
        assert set(copied) == set(p.name for p in pathlib.Path(".").iterdir())


def test_copy_dry_run(tmpdir):
    with tmpdir.as_cwd():
        to_copy = copy_example_data("io/shapefile", case="methow", dry_run=True)
        assert len(list(p.name for p in pathlib.Path(".").iterdir())) == 0

        copied = copy_example_data("io/shapefile", case="methow", dry_run=False)
        assert set(to_copy) == set(copied)
        assert set(to_copy) == set(p.name for p in pathlib.Path(".").iterdir())


def test_copy_no_overwrite(tmpdir):
    with tmpdir.as_cwd():
        to_copy = copy_example_data("io/shapefile", case="methow", dry_run=True)
        pathlib.Path(to_copy[0]).touch()
        with pytest.raises(FileExistsError):
            copy_example_data("io/shapefile", case="methow", dry_run=False)
