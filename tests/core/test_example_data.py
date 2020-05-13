import pathlib

import pytest

from landlab import ExampleData


@pytest.mark.parametrize("case", (None, "", "methow"))
def test_fetch(tmpdir, case):
    if case is None:
        data = ExampleData("io/shapefile")
    else:
        data = ExampleData("io/shapefile", case=case)
    expected = set(data)
    with tmpdir.as_cwd():
        data.fetch()

        fetched = set(p.name for p in pathlib.Path(".").iterdir())
        assert expected == fetched


def test_fetch_all_files(tmpdir):
    with tmpdir.as_cwd():
        ExampleData("io/shapefile", case="methow").fetch()
        fetched_with_case = set(p.name for p in pathlib.Path(".").iterdir())

        ExampleData("io/shapefile").fetch()
        fetched_without_case = set(p.name for p in pathlib.Path("methow").iterdir())

        assert fetched_with_case == fetched_without_case


def test_fetch_no_overwrite(tmpdir):
    data = ExampleData("io/shapefile")
    with tmpdir.as_cwd():
        pathlib.Path(list(data).pop()).touch()
        with pytest.raises(FileExistsError):
            data.fetch()
