#! /usr/bin/env python
"""
Unit tests for landlab.collections
"""
from nose.tools import assert_equal, assert_true, assert_raises
try:
    from nose.tools import (assert_dict_equal, assert_tuple_equal,
                            assert_list_equal)
except ImportError:
    from landlab.testing.tools import (assert_dict_equal, assert_tuple_equal,
                                       assert_list_equal)

from landlab import Palette, Arena, NoProvidersError

from landlab import Implements, ImplementsOrRaise
from landlab.framework.interfaces import BmiBase, BmiNoGrid


@Implements(BmiBase)
class Sample1(object):

    """A sample component."""

    __implements__ = (BmiBase, BmiNoGrid, )

    _input_var_names = [
        'air__temperature',
        'surface__elevation',
    ]
    _output_var_names = [
        'deposition__rate',
    ]

    model_name = 'Sample 1'
    author_name = 'Eric Hutton'
    version = '0.1'
    time_units = 's'
    time_step_type = 'fixed'
    step_method = 'explicit'
    grid_type = 'none'

    _vars = {
        'deposition__rate': [1.]
    }

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
        return 0.

    def get_current_time(self):
        return 0.

    def get_end_time(self):
        return 100.

    def get_time_step(self):
        return 1.

    def get_var_type(self, name):
        return 'float64'

    def get_var_units(self, name):
        return 'm'

    def set_value(self, name, value):
        pass

    def get_value(self, name):
        return self._vars[name]


@ImplementsOrRaise(BmiBase)
class Sample2(object):

    """A sample component."""

    __implements__ = (BmiBase, BmiNoGrid, )

    _input_var_names = [
        'deposition__rate',
    ]
    _output_var_names = [
        'air__temperature',
        'surface__elevation',
    ]

    model_name = 'Sample 2'
    author_name = 'Eric Hutton'
    version = '0.1'
    time_units = 's'
    time_step_type = 'fixed'
    step_method = 'explicit'
    grid_type = 'none'

    _vars = {
        'air__temperature': [1.],
        'surface__elevation': [1.],
    }

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
        return 0.

    def get_current_time(self):
        return 0.

    def get_end_time(self):
        return 100.

    def get_time_step(self):
        return 1.

    def get_var_type(self, name):
        return 'float64'

    def get_var_units(self, name):
        return 'm'

    def get_value(self, name):
        return self._vars[name]

    def set_value(self, name, value):
        pass


def test_empty_palette():
    """Create a palette without components."""
    palette = Palette()

    assert_equal(len(palette), 0)
    assert_equal(list(palette.list()), [])
    assert_equal(list(palette.keys()), [])
    assert_equal(palette.uses(), [])
    assert_equal(palette.provides(), [])

    providers = palette.find_provider('air__temperature')
    assert_equal(providers, [])

    users = palette.find_user('air__temperature')
    assert_equal(users, [])

    connections = palette.find_connections()
    assert_dict_equal(connections, {})


def test_1_component_create():
    palette = Palette(sample=Sample1)

    assert_equal(len(palette), 1)


def test_1_component_dict_interface():
    palette = Palette(sample=Sample1)

    assert_dict_equal(dict(sample=Sample1), palette)
    assert_equal(len(palette), 1)
    assert_list_equal(list(palette.keys()), ['sample'])
    assert_list_equal(list(palette.values()), [Sample1])

    items = list(palette.items())
    assert_tuple_equal(('sample', Sample1), items[0])


def test_1_component_list():
    palette = Palette(sample=Sample1)

    assert_equal(['sample'], list(palette.list()))


def test_1_component_uses():
    palette = Palette(sample=Sample1)

    uses = palette.uses()
    uses.sort()
    assert_equal(['air__temperature', 'surface__elevation'], uses)


def test_1_component_provides():
    palette = Palette(sample=Sample1)

    provides = palette.provides()
    provides.sort()
    assert_equal(['deposition__rate'], provides)


def test_1_component_find_providers():
    palette = Palette(sample=Sample1)

    providers = palette.find_provider('air__temperature')
    assert_equal(providers, [])

    providers = palette.find_provider('deposition__rate')
    assert_equal(providers, ['sample'])


def test_1_component_find_users():
    palette = Palette(sample=Sample1)

    users = palette.find_user('air__temperature')
    assert_equal(users, ['sample'])


def test_1_component_find_connections():
    palette = Palette(sample=Sample1)

    assert_raises(NoProvidersError, palette.find_connections)


def test_2_components_create():
    palette = Palette(one=Sample1, two=Sample2)

    assert_equal(len(palette), 2)


def test_2_components_dict_interface():
    palette = Palette(one=Sample1, two=Sample2)

    assert_dict_equal(dict(one=Sample1, two=Sample2), palette)
    assert_equal(len(palette), 2)

    keys = list(palette.keys())
    keys.sort()
    assert_list_equal(['one', 'two'], keys)

    values = palette.values()
    assert_equal(2, len(values))
    assert_true(Sample1 in values and Sample2 in values)

    items = list(palette.items())
    items.sort()
    assert_equal(2, len(items))
    assert_tuple_equal(('one', Sample1), items[0])
    assert_tuple_equal(('two', Sample2), items[1])


def test_2_components_list():
    palette = Palette(one=Sample1, two=Sample2)

    components = list(palette.list())
    components.sort()
    assert_list_equal(['one', 'two'], components)


def test_2_components_uses():
    palette = Palette(one=Sample1, two=Sample2)

    uses = palette.uses()
    uses.sort()
    assert_list_equal(['air__temperature',
                       'deposition__rate',
                       'surface__elevation'], uses)


def test_2_components_provides():
    palette = Palette(one=Sample1, two=Sample2)

    provides = palette.provides()
    provides.sort()
    assert_list_equal(['air__temperature',
                       'deposition__rate',
                       'surface__elevation'], provides)


def test_2_components_find_providers():
    palette = Palette(one=Sample1, two=Sample2)

    providers = palette.find_provider('air__temperature')
    assert_list_equal(['two'], providers)

    providers = palette.find_provider('deposition__rate')
    assert_equal(['one'], providers)


def test_2_components_find_users():
    palette = Palette(one=Sample1, two=Sample2)

    users = palette.find_user('air__temperature')
    assert_list_equal(['one'], users)


def test_2_components_find_connections():
    palette = Palette(one=Sample1, two=Sample2)

    connections = {
        'one': {'deposition__rate': ['two']},
        'two': {'air__temperature': ['one'],
                'surface__elevation': ['one']},
    }
    assert_dict_equal(connections, palette.find_connections())


def test_arena_instantiate():
    arena = Arena()
    assert_dict_equal(dict(), arena)

    arena.instantiate(Sample1, 'one')
    assert_equal(1, len(arena))
    assert_true('one' in arena)
    assert_true(isinstance(arena['one'], Sample1))

    arena.instantiate(Sample2, 'two')
    assert_equal(2, len(arena))
    assert_true('one' in arena and 'two' in arena)
    assert_true(isinstance(arena['one'], Sample1))
    assert_true(isinstance(arena['two'], Sample2))


def test_arena_connect():
    arena = Arena()
    arena.instantiate(Sample1, 'one')
    arena.instantiate(Sample2, 'two')

    arena.connect('one', 'two', 'air__temperature')
    arena.connect('one', 'two', 'surface__elevation')
    arena.connect('two', 'one', 'deposition__rate')

    dz = arena['one'].get_value('deposition__rate')
    t = arena['two'].get_value('air__temperature')
    z = arena['two'].get_value('surface__elevation')


def test_arena_walk():
    arena = Arena()
    arena.instantiate(Sample1, 'one')
    arena.instantiate(Sample2, 'two')

    arena.connect('one', 'two', 'air__temperature')
    arena.connect('one', 'two', 'surface__elevation')
    arena.connect('two', 'one', 'deposition__rate')

    tree = arena.walk('one')
    assert_list_equal(['one', 'two'], tree)
