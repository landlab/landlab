#! /usr/bin/env python
"""
Unit tests for landlab.collections
"""
import pytest

from landlab import Arena, Implements, ImplementsOrRaise, NoProvidersError, Palette
from landlab.framework.interfaces import BmiBase, BmiNoGrid


@Implements(BmiBase)
class Sample1(object):

    """A sample component."""

    __implements__ = (BmiBase, BmiNoGrid)

    _input_var_names = ["air__temperature", "surface__elevation"]
    _output_var_names = ["deposition__rate"]

    model_name = "Sample 1"
    author_name = "Eric Hutton"
    version = "0.1"
    time_units = "s"
    time_step_type = "fixed"
    step_method = "explicit"
    grid_type = "none"

    _vars = {"deposition__rate": [1.0]}

    def initialize(self, name):
        pass

    def update(self):
        pass

    def finalize(self):
        pass

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_var_rank(self, name):
        return 0

    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        return 0.0

    def get_end_time(self):
        return 100.0

    def get_time_step(self):
        return 1.0

    def get_var_type(self, name):
        return "float64"

    def get_var_units(self, name):
        return "m"

    def set_value(self, name, value):
        pass

    def get_value(self, name):
        return self._vars[name]


@ImplementsOrRaise(BmiBase)
class Sample2(object):

    """A sample component."""

    __implements__ = (BmiBase, BmiNoGrid)

    _input_var_names = ["deposition__rate"]
    _output_var_names = ["air__temperature", "surface__elevation"]

    model_name = "Sample 2"
    author_name = "Eric Hutton"
    version = "0.1"
    time_units = "s"
    time_step_type = "fixed"
    step_method = "explicit"
    grid_type = "none"

    _vars = {"air__temperature": [1.0], "surface__elevation": [1.0]}

    def initialize(self, name):
        pass

    def update(self):
        pass

    def finalize(self):
        pass

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_var_rank(self, name):
        return 0

    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        return 0.0

    def get_end_time(self):
        return 100.0

    def get_time_step(self):
        return 1.0

    def get_var_type(self, name):
        return "float64"

    def get_var_units(self, name):
        return "m"

    def get_value(self, name):
        return self._vars[name]

    def set_value(self, name, value):
        pass


def test_empty_palette():
    """Create a palette without components."""
    palette = Palette()

    assert len(palette) == 0
    assert list(palette.list()) == []
    assert list(palette.keys()) == []
    assert palette.uses() == []
    assert palette.provides() == []

    providers = palette.find_provider("air__temperature")
    assert providers == []

    users = palette.find_user("air__temperature")
    assert users == []

    connections = palette.find_connections()
    assert connections == {}


def test_1_component_create():
    palette = Palette(sample=Sample1)

    assert len(palette) == 1


def test_1_component_dict_interface():
    palette = Palette(sample=Sample1)

    assert dict(sample=Sample1) == palette
    assert len(palette) == 1
    assert list(palette.keys()) == ["sample"]
    assert list(palette.values()) == [Sample1]

    items = list(palette.items())
    assert ("sample", Sample1) == items[0]


def test_1_component_list():
    palette = Palette(sample=Sample1)

    assert ["sample"] == list(palette.list())


def test_1_component_uses():
    palette = Palette(sample=Sample1)

    uses = palette.uses()
    uses.sort()
    assert ["air__temperature", "surface__elevation"] == uses


def test_1_component_provides():
    palette = Palette(sample=Sample1)

    provides = palette.provides()
    provides.sort()
    assert ["deposition__rate"] == provides


def test_1_component_find_providers():
    palette = Palette(sample=Sample1)

    providers = palette.find_provider("air__temperature")
    assert providers == []

    providers = palette.find_provider("deposition__rate")
    assert providers == ["sample"]


def test_1_component_find_users():
    palette = Palette(sample=Sample1)

    users = palette.find_user("air__temperature")
    assert users == ["sample"]


def test_1_component_find_connections():
    palette = Palette(sample=Sample1)

    with pytest.raises(NoProvidersError):
        palette.find_connections()


def test_2_components_create():
    palette = Palette(one=Sample1, two=Sample2)

    assert len(palette) == 2


def test_2_components_dict_interface():
    palette = Palette(one=Sample1, two=Sample2)

    assert dict(one=Sample1, two=Sample2) == palette
    assert len(palette) == 2

    keys = list(palette.keys())
    keys.sort()
    assert ["one", "two"] == keys

    values = palette.values()
    assert 2 == len(values)
    assert Sample1 in values and Sample2 in values

    items = list(palette.items())
    items.sort()
    assert 2 == len(items)
    assert ("one", Sample1) == items[0]
    assert ("two", Sample2) == items[1]


def test_2_components_list():
    palette = Palette(one=Sample1, two=Sample2)

    components = list(palette.list())
    components.sort()
    assert ["one", "two"] == components


def test_2_components_uses():
    palette = Palette(one=Sample1, two=Sample2)

    uses = palette.uses()
    uses.sort()
    assert ["air__temperature", "deposition__rate", "surface__elevation"] == uses


def test_2_components_provides():
    palette = Palette(one=Sample1, two=Sample2)

    provides = palette.provides()
    provides.sort()
    assert ["air__temperature", "deposition__rate", "surface__elevation"] == provides


def test_2_components_find_providers():
    palette = Palette(one=Sample1, two=Sample2)

    providers = palette.find_provider("air__temperature")
    assert ["two"] == providers

    providers = palette.find_provider("deposition__rate")
    assert ["one"] == providers


def test_2_components_find_users():
    palette = Palette(one=Sample1, two=Sample2)

    users = palette.find_user("air__temperature")
    assert ["one"] == users


def test_2_components_find_connections():
    palette = Palette(one=Sample1, two=Sample2)

    connections = {
        "one": {"deposition__rate": ["two"]},
        "two": {"air__temperature": ["one"], "surface__elevation": ["one"]},
    }
    assert connections == palette.find_connections()


def test_arena_instantiate():
    arena = Arena()
    assert dict() == arena

    arena.instantiate(Sample1, "one")
    assert 1 == len(arena)
    assert "one" in arena
    assert isinstance(arena["one"], Sample1)

    arena.instantiate(Sample2, "two")
    assert 2 == len(arena)
    assert "one" in arena and "two" in arena
    assert isinstance(arena["one"], Sample1)
    assert isinstance(arena["two"], Sample2)


def test_arena_connect():
    arena = Arena()
    arena.instantiate(Sample1, "one")
    arena.instantiate(Sample2, "two")

    arena.connect("one", "two", "air__temperature")
    arena.connect("one", "two", "surface__elevation")
    arena.connect("two", "one", "deposition__rate")

    arena["one"].get_value("deposition__rate")
    arena["two"].get_value("air__temperature")
    arena["two"].get_value("surface__elevation")


def test_arena_walk():
    arena = Arena()
    arena.instantiate(Sample1, "one")
    arena.instantiate(Sample2, "two")

    arena.connect("one", "two", "air__temperature")
    arena.connect("one", "two", "surface__elevation")
    arena.connect("two", "one", "deposition__rate")

    tree = arena.walk("one")
    assert ["one", "two"] == tree
