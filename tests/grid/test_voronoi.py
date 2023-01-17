import pathlib

import numpy as np
import pytest

from landlab import VoronoiDelaunayGrid


def test_save(tmpdir):
    with tmpdir.as_cwd():
        grid = VoronoiDelaunayGrid(np.random.rand(20), np.random.rand(20))
        filename = grid.save("mytestsave.grid")

        assert pathlib.Path(filename).is_file()


def test_save_with_clobber(tmpdir):
    with tmpdir.as_cwd():
        pathlib.Path("mytestsave.grid").touch()

        grid = VoronoiDelaunayGrid(np.random.rand(20), np.random.rand(20))
        with pytest.raises(ValueError):
            grid.save("mytestsave.grid")
        filename = grid.save("mytestsave.grid", clobber=True)

        assert pathlib.Path(filename).is_file()
        assert pathlib.Path(filename).stat().st_size > 0


def test_save_with_extension(tmpdir):
    with tmpdir.as_cwd():
        grid = VoronoiDelaunayGrid(np.random.rand(20), np.random.rand(20))
        filename = grid.save("mytestsave", clobber=True)

        assert filename == "mytestsave.grid"
        assert pathlib.Path(filename).is_file()

        filename = grid.save("mytestsave.foo", clobber=True)
        assert filename == "mytestsave.foo.grid"
        assert pathlib.Path(filename).is_file()
