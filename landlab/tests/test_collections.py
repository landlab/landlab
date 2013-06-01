#! /usr/bin/env python
"""
Unit tests for landlab.collections
"""

import unittest

from landlab import Palette, Arena, NoProvidersError
from landlab.components.sample import Sample1, Sample2

class TestLandlabPalette(unittest.TestCase):
    def test_empty_palette(self):
        """
        Create a palette without components.
        """
        palette = Palette()

        self.assertEqual(len(palette), 0)
        self.assertEqual(palette.list(), [])
        self.assertEqual(palette.keys(), [])
        self.assertEqual(palette.uses(), [])
        self.assertEqual(palette.provides(), [])

        providers = palette.find_provider('air__temperature')
        self.assertEqual(providers, [])

        users = palette.find_user('air__temperature')
        self.assertEqual(users, [])

        connections = palette.find_connections()
        self.assertDictEqual(connections, {})


class TestLandlabOneComponentPalette(unittest.TestCase):
    def test_create(self):
        palette = Palette(sample=Sample1)

        self.assertEqual(len(palette), 1)

    def test_dict_interface(self):
        palette = Palette(sample=Sample1)

        self.assertDictEqual(dict(sample=Sample1), palette)
        self.assertEqual(len(palette), 1)
        self.assertEqual(palette.keys(), ['sample'])
        self.assertEqual(palette.values(), [Sample1])

        items = palette.items()
        self.assertTupleEqual(('sample', Sample1), items[0])

    def test_list(self):
        palette = Palette(sample=Sample1)

        self.assertEqual(['sample'], palette.list())

    def test_uses(self):
        palette = Palette(sample=Sample1)

        uses = palette.uses()
        uses.sort()
        self.assertEqual(['air__temperature', 'surface__elevation'], uses)

    def test_provides(self):
        palette = Palette(sample=Sample1)

        provides = palette.provides()
        provides.sort()
        self.assertEqual(['deposition__rate'], provides)

    def test_find_providers(self):
        palette = Palette(sample=Sample1)

        providers = palette.find_provider('air__temperature')
        self.assertEqual(providers, [])

        providers = palette.find_provider('deposition__rate')
        self.assertEqual(providers, ['sample'])

    def test_find_users(self):
        palette = Palette(sample=Sample1)

        users = palette.find_user('air__temperature')
        self.assertEqual(users, ['sample'])

    def test_find_connections(self):
        palette = Palette(sample=Sample1)

        with self.assertRaises(NoProvidersError):
            connections = palette.find_connections()


class TestLandlabTwoComponentPalette(unittest.TestCase):
    def test_create(self):
        palette = Palette(one=Sample1, two=Sample2)

        self.assertEqual(len(palette), 2)

    def test_dict_interface(self):
        palette = Palette(one=Sample1, two=Sample2)

        self.assertDictEqual(dict(one=Sample1, two=Sample2), palette)
        self.assertEqual(len(palette), 2)

        keys = palette.keys()
        keys.sort()
        self.assertListEqual(['one', 'two'], keys)

        values = palette.values()
        self.assertEqual(2, len(values))
        self.assertTrue(Sample1 in values and Sample2 in values)

        items = palette.items()
        items.sort()
        self.assertEqual(2, len(items))
        self.assertTupleEqual(('one', Sample1), items[0])
        self.assertTupleEqual(('two', Sample2), items[1])

    def test_list(self):
        palette = Palette(one=Sample1, two=Sample2)

        list = palette.list()
        list.sort()
        self.assertListEqual(['one', 'two'], list)

    def test_uses(self):
        palette = Palette(one=Sample1, two=Sample2)

        uses = palette.uses()
        uses.sort()
        self.assertListEqual(['air__temperature',
                              'deposition__rate',
                              'surface__elevation'], uses)

    def test_provides(self):
        palette = Palette(one=Sample1, two=Sample2)

        provides = palette.provides()
        provides.sort()
        self.assertListEqual(['air__temperature',
                              'deposition__rate',
                              'surface__elevation'], provides)

    def test_find_providers(self):
        palette = Palette(one=Sample1, two=Sample2)

        providers = palette.find_provider('air__temperature')
        self.assertListEqual(['two'], providers)

        providers = palette.find_provider('deposition__rate')
        self.assertEqual(['one'], providers)

    def test_find_users(self):
        palette = Palette(one=Sample1, two=Sample2)

        users = palette.find_user('air__temperature')
        self.assertListEqual(['one'], users)

    def test_find_connections(self):
        palette = Palette(one=Sample1, two=Sample2)

        connections = {
            'one': {'deposition__rate': ['two']},
            'two': {'air__temperature': ['one'],
                    'surface__elevation': ['one']},
        }
        self.assertDictEqual(connections, palette.find_connections())


class TestLandlabArena(unittest.TestCase):
    def test_instantiate(self):
        arena = Arena()
        self.assertDictEqual(dict(), arena)

        arena.instantiate(Sample1, 'one')
        self.assertEqual(1, len(arena))
        self.assertTrue('one' in arena)
        self.assertTrue(isinstance(arena['one'], Sample1))

        arena.instantiate(Sample2, 'two')
        self.assertEqual(2, len(arena))
        self.assertTrue('one' in arena and 'two' in arena)
        self.assertTrue(isinstance(arena['one'], Sample1))
        self.assertTrue(isinstance(arena['two'], Sample2))

    def test_connect(self):
        arena = Arena()
        arena.instantiate(Sample1, 'one')
        arena.instantiate(Sample2, 'two')

        arena.connect('one', 'two', 'air__temperature')
        arena.connect('one', 'two', 'surface__elevation')
        arena.connect('two', 'one', 'deposition__rate')

        dz = arena['one'].get_value('deposition__rate')
        t = arena['two'].get_value('air__temperature')
        z = arena['two'].get_value('surface__elevation')

    def test_walk(self):
        arena = Arena()
        arena.instantiate(Sample1, 'one')
        arena.instantiate(Sample2, 'two')

        arena.connect('one', 'two', 'air__temperature')
        arena.connect('one', 'two', 'surface__elevation')
        arena.connect('two', 'one', 'deposition__rate')

        tree = arena.walk('one')
        self.assertListEqual(['one', 'two'], tree)

if __name__ == '__main__':
    unittest.main()
