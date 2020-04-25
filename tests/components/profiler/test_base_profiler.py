import pytest

from landlab.components.profiler.base_profiler import _BaseProfiler


@pytest.fixture()
def EmptyBaseProfiler():
    class EmptyBaseProfiler(_BaseProfiler):
        def __init__(self):
            pass

        def _create_profile_structure(self):
            pass

        @property
        def distance_along_profile(self):
            return 1

        @property
        def nodes(self):
            return 2

    return EmptyBaseProfiler


def test_bad_base():
    class Base(_BaseProfiler):
        pass

    with pytest.raises(TypeError):
        Base()


def test_base_implemented(EmptyBaseProfiler):
    empty_base = EmptyBaseProfiler()
    assert isinstance(empty_base, _BaseProfiler)
    empty_base._create_profile_structure()
    assert empty_base.distance_along_profile == 1
    assert empty_base.nodes == 2
