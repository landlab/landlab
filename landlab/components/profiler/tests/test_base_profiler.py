import pytest

from landlab.components.profiler.base_profiler import _BaseProfiler


@pytest.fixture()
def EmptyBaseProfiler():
    class EmptyBaseProfiler(_BaseProfiler):
        def __init__(self):
            pass

        def _calculate_distances(self):
            pass

    return EmptyBaseProfiler


def test_bad_base():
    class Base(_BaseProfiler):
        pass

    with pytest.raises(TypeError):
        Base()


def test_base_implemented(EmptyBaseProfiler):
    empty_base = EmptyBaseProfiler()
    assert isinstance(empty_base, _BaseProfiler)
    empty_base._calculate_distances()
