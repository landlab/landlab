#! /usr/bin/env python

from nose.tools import assert_true, assert_false, assert_raises
try:
    from nose.tools import assert_is_not, assert_is, assert_set_equal, assert_dict_equal
except ImportError:
    from landlab.testing.tools import (assert_is_not, assert_is, assert_set_equal,
                                       assert_dict_equal)
import numpy as np
from numpy.testing import assert_array_equal

from landlab.field.graph_field import GraphFields as ModelDataFields
from landlab.field import GroupError, FieldError


def test_init():
    fields = ModelDataFields()
    assert_set_equal(set(), fields.groups)


def test_new_field_location():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    assert_set_equal(set(['node']), fields.groups)


def test_add_existing_group():
    fields = ModelDataFields()
    fields.new_field_location('node', size=12)
    assert_raises(ValueError, fields.new_field_location, 'node', size=24)


def test_add_multiple_groups():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    fields.new_field_location('cell', 2)
    fields.new_field_location('face', 7)
    fields.new_field_location('link', 7)
    assert_set_equal(set(['node', 'cell', 'face', 'link']), fields.groups)


def test_ones():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    fields.new_field_location('cell', 2)

    value_array = fields.ones('node')
    assert_array_equal(np.ones(12), value_array)

    value_array = fields.ones('cell')
    assert_array_equal(np.ones(2), value_array)


def test_add_ones():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    fields.new_field_location('cell', 2)

    fields.add_ones('z', at='node')
    assert_array_equal(np.ones(12), fields['node']['z'])
    assert_array_equal(np.ones(12), fields.field_values('node', 'z'))

    fields.add_ones('z', at='cell')
    assert_array_equal(np.ones(2), fields['cell']['z'])
    assert_array_equal(np.ones(2), fields.field_values('cell', 'z'))


def test_add_ones_return_value():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    fields.new_field_location('cell', 2)

    rtn_value = fields.add_ones('z', at='node')
    assert_array_equal(rtn_value, np.ones(12))
    assert_is(rtn_value, fields['node']['z'])
    assert_is(rtn_value, fields.field_values('node', 'z'))

    rtn_value = fields.add_ones('z', at='cell')
    assert_array_equal(rtn_value, np.ones(2))
    assert_is(rtn_value, fields['cell']['z'])
    assert_is(rtn_value, fields.field_values('cell', 'z'))


def test_add_existing_field_default():
    """Test default is to not replace existing field."""
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    fields.add_empty('z', at='node')

    assert_raises(FieldError, fields.add_empty, 'z', at='node')
    assert_raises(FieldError, fields.add_ones, 'z', at='node')
    assert_raises(FieldError, fields.add_zeros, 'z', at='node')


def test_add_existing_field_with_noclobber():
    """Test noclobber raises an error with an existing field."""
    fields = ModelDataFields()
    fields.new_field_location('node', 12)
    fields.add_empty('z', at='node')

    assert_raises(FieldError, fields.add_empty, 'z', at='node', noclobber=True)
    assert_raises(FieldError, fields.add_ones, 'z', at='node', noclobber=True)
    assert_raises(FieldError, fields.add_zeros, 'z', at='node', noclobber=True)


def test_add_field_with_noclobber():
    """Test noclobber does not raise an error with an new field."""
    fields = ModelDataFields()
    fields.new_field_location('node', 12)

    fields.add_empty('a', at='node', noclobber=True)
    assert_true('a' in fields['node'])

    fields.add_ones('b', at='node', noclobber=True)
    assert_true('b' in fields['node'])

    fields.add_zeros('c', at='node', noclobber=True)
    assert_true('c' in fields['node'])


def test_add_field_with_clobber():
    """Test adding a field with clobber on."""
    fields = ModelDataFields()
    fields.new_field_location('node', 12)

    assert_is_not(fields.add_empty('a', at='node'),
                  fields.add_empty('a', at='node', noclobber=False))
    assert_is_not(fields.add_ones('b', at='node'),
                  fields.add_ones('b', at='node', noclobber=False))
    assert_is_not(fields.add_zeros('c', at='node'),
                  fields.add_zeros('c', at='node', noclobber=False))


def test_getitem():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)

    assert_dict_equal(dict(), fields['node'])
    assert_raises(GroupError, lambda k: fields[k], 'cell')
    assert_raises(KeyError, lambda k: fields[k], 'cell')


def test_at_attribute():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)

    assert_dict_equal(dict(), fields.at_node)
    assert_raises(AttributeError, lambda: fields.at_cell)

    fields.add_ones('z', at='node')
    assert_array_equal(np.ones(12), fields.at_node['z'])


def test_has_group():
    fields = ModelDataFields()
    fields.new_field_location('node', 12)

    assert_true(fields.has_group('node'))
    assert_false(fields.has_group('cell'))


def test_delete_field():
    fields = ModelDataFields()
    fields.new_field_location('link', 17)

    assert_dict_equal(dict(), fields.at_link)
    assert_raises(AttributeError, lambda: fields.at_node)

    fields.add_zeros('link', 'vals')
    assert_array_equal(np.zeros(17), fields.at_link['vals'])

    assert_raises(KeyError, lambda: fields.delete_field('node', 'vals'))
    fields.delete_field('link', 'vals')
    assert_raises(KeyError, lambda: fields.field_units('link', 'vals'))
    assert_raises(KeyError, lambda: fields.at_link['vals'])


def test_scalar_field():
    """Test adding a generic scalar field."""
    fields = ModelDataFields()
    fields.new_field_location('all_over_the_place', 1)

    assert_dict_equal(dict(), fields.at_all_over_the_place)
    assert_raises(AttributeError, lambda: fields.at_cell)

    fields.at_all_over_the_place['const'] = 1.
    assert_array_equal(np.array(1.), fields.at_all_over_the_place['const'])
    
    val = np.array(2.)
    fields.at_all_over_the_place['const'] = val
    assert_is(val, fields.at_all_over_the_place['const'])


def test_grid_field_as_array():
    """Test adding an array as a grid field."""
    fields = ModelDataFields()
    fields.new_field_location('grid', 1)

    fields.at_grid['const'] = [1., 2.]
    assert_array_equal(np.array([1., 2.]), fields.at_grid['const'])

    val = np.array([1., 2.])
    fields.at_grid['const'] = val
    assert_is(val, fields.at_grid['const'])

    val.shape = (1, 1, 2, 1)
    fields.at_grid['const'] = val
    assert_array_equal(np.array([1., 2.]), fields.at_grid['const'])
    assert_is(val, fields.at_grid['const'])


def test_grid_field_add_zeros_ones_empty():
    """Test creating scalar fields with add_zeros, add_empty, and add_ones."""
    fields = ModelDataFields()
    fields.new_field_location('grid', 1)
    
    assert_raises(ValueError, fields.add_zeros, 'value', at='grid')
    assert_raises(ValueError, fields.add_empty, 'value', at='grid')
    assert_raises(ValueError, fields.add_ones, 'value', at='grid')

def test_grid_field_zeros_ones_empty():
    """Test creating scalar fields with zeros, empty, and ones."""
    fields = ModelDataFields()
    fields.new_field_location('grid', 1)    
    assert_raises(ValueError, fields.zeros, 'grid')
    assert_raises(ValueError, fields.empty, 'grid')
    assert_raises(ValueError, fields.ones, 'grid')

def test_nd_field():
    """Test creating fields that are nd in shape."""
    fields = ModelDataFields()
    fields.new_field_location('node', 12)

    fields.add_field('new_value', np.ones((12,4,5)), at='node')
    fields.add_field('newer_value', np.ones((12,4)), at='node')

    assert_raises(ValueError, fields.add_field, 'newest_value', np.ones((13,4,5)), at='node')
    assert_raises(ValueError, fields.add_field, 'newestest_value', np.ones((13)), at='node')

    
