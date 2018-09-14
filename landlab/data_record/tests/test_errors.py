
# -*- coding: utf-8 -*-
"""

Unit tests for landlab.data_record.data_record.DataRecord
Test errors/exceptions

Last updated 8/27/2018


"""

import pytest
import numpy as np
from landlab import (RasterModelGrid,
                     HexModelGrid,
                     RadialModelGrid,
                     VoronoiDelaunayGrid)
from landlab.data_record import DataRecord

grid = RasterModelGrid((3,3))


def dr_bad_dim():
    time=[0.]
    data_vars={'mean_elev' : (['time'], [100.]), 'test' : (['bad_dim'], [12])}
    with pytest.raises(ValueError):
        DataRecord(grid,
                   time=time,
                   data_vars=data_vars)
# should return ValueError('Data variable dimensions must be time and/or'
#                              'item_id')
def dr_bad_time():
    bad_time = 'bad_time'
    with pytest.raises(TypeError):
        DataRecord(grid=grid, time=bad_time)
#should return TypeError: Time must be a list or an array of length 1

def dr_bad_datavars():
    time=[0.]
    data_vars=['not a dict']
    with pytest.raises(TypeError):
        DataRecord(grid, time=time, data_vars=data_vars)
# should return TypeError(('Data variables (data_vars) passed to'
#                                 ' DataRecord must be a dictionary (see '
#                                 'documentation for valid structure)'))

def dr_bad_items():
    time=[0.]
    data_vars={'mean_elevation' : (['time'], np.array([100]))}
    attrs={'time_units' : 'y'}
    my_items_bad = ['not a dict']
    with pytest.raises(KeyError):
        DataRecord(grid,
                   time=time,
                   items=my_items_bad,
                   data_vars=data_vars,
                   attrs=attrs)
#Should return KeyError(('You must provide an ''items'' dictionary '
#                                 '(see documentation for required format)'))

def dr_bad_items_keys():
    time=[0.]
    data_vars={'mean_elevation' : (['time'], np.array([100]))}
    attrs={'time_units' : 'y'}
    my_items_bad = {'grid_element':np.array([['node'], ['link']]),
                    'bad_key': np.array([[1],[3]])}
    with pytest.raises(KeyError):
        DataRecord(grid,
                   time=time,
                   items=my_items_bad,
                   data_vars=data_vars,
                   attrs=attrs)
#Should return KeyError(('You must provide an ''items'' dictionary '
#                                 '(see documentation for required format)'))

def dr_bad_attrs():
    time=[0.]
    data_vars={'mean_elevation' : (['time'], np.array([100]))}
    attrs=['not a dict']
    my_items3 = {'grid_element':np.array([['node'], ['link']]),
                 'element_id': np.array([[1],[3]])}
    with pytest.raises(AttributeError):
        DataRecord(grid,
                   time=time,
                   items=my_items3,
                   data_vars=data_vars,
                   attrs=attrs)
#Should return except AttributeError:
#                raise TypeError(('Attributes (attrs) passed to DataRecord'
#                                'must be a dictionary'))

def dr_bad_loc_and_id():
    rmg = RasterModelGrid((3,3)),
    hmg = HexModelGrid(3, 2, 1.0),
    radmg = RadialModelGrid(num_shells=1, dr=1., origin_x=0., origin_y=0.)
    vdmg = VoronoiDelaunayGrid(np.random.rand(25), np.random.rand(25))
    my_items_bad_loc = {'grid_element':np.array(['node', 'bad_loc']),
                        'element_id': np.array([1,3])}
    my_items_bad_loc2 = {'grid_element': 'bad_loc',
                         'element_id': np.array([1,3])}
    my_items_bad_loc3 = {'grid_element': np.array(['node', 'node', 'node']),
                         'element_id': np.array([1,3])}
    my_items_bad_id = {'grid_element':np.array(['node', 'link']),
                       'element_id': np.array([1,300])}
    my_items_bad_id2 = {'grid_element':np.array(['node', 'link']),
                        'element_id': np.array([1,-300])}
    my_items_bad_id3 = {'grid_element':np.array(['node', 'link']),
                        'element_id': np.array([1,2.])}
    my_items_bad_id4 = {'grid_element':np.array(['node', 'link']),
                        'element_id': np.array([1, 2, 3])}
    with pytest.raises(ValueError):
        DataRecord(rmg,items=my_items_bad_loc)
        DataRecord(hmg,items=my_items_bad_loc)
        DataRecord(radmg,items=my_items_bad_loc)
        DataRecord(vdmg,items=my_items_bad_loc)
#should return ValueError(('One or more of the grid elements provided is/are'
#                            ' not permitted location for this grid type.'))
        DataRecord(rmg,items=my_items_bad_loc2)
#should return ValueError: Location provided: bad_loc is not a permitted
# location for this grid type.
        DataRecord(rmg,items=my_items_bad_loc3)
#should return ValueError(('grid_element passed to DataRecord must be '
# ' the same length as the number of items or 1.'))
        DataRecord(rmg,items=my_items_bad_id)
        DataRecord(hmg,items=my_items_bad_id)
        DataRecord(radmg,items=my_items_bad_id)
        DataRecord(vdmg,items=my_items_bad_id)
#should return ValueError(('An item residing at ' + at + ' has an '
# 'element_id larger than the number of'+ at + 'on the grid.'))
        DataRecord(rmg,items=my_items_bad_id2)
# should return ValueError(('An item residing at ' + at + ' has '
#'an element id below zero. This is not permitted.'))
        DataRecord(rmg,items=my_items_bad_id3)
#should return ValueError(('You have passed a non-integer element_id to '
#                             'DataRecord, this is not permitted.'))
        DataRecord(rmg,items=my_items_bad_id4)
#should return ValueError(('The number of grid_element passed '
# ' to Datarecord must be 1 or equal  to the number of element_id '))

########### TIME ONLY #########################################################
def dr_time_bad_add_record(dr_time):
    with pytest.raises(TypeError):
        dr_time.add_record(time = 'bad_time', new_record=None)
# should return TypeError: You have passed a time that is not permitted,
# must be list or array
    with pytest.raises(KeyError):
        dr_time.add_record(time = [300.], item_id=[0], new_record=None)
# should return KeyError: 'This Datarecord does not hold items'

def dr_time_bad_get_data_(dr_time):
    with pytest.raises(KeyError):
        dr_time.get_data(time=0., data_variable='bad_variable')
#should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
        dr_time.get_data(time=0.,
                         item_id=0.,
                         data_variable='mean_elevation')
#should return KeyError: 'This Datarecord does not hold items'

def dr_time_bad_set_data_(dr_time):
    with pytest.raises(KeyError):
        dr_time.set_data(time=0.,
                         data_variable='bad_variable',
                         new_value=105.)
#should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
        dr_time.set_data(time=10.,
                         data_variable='mean_elevation',
                         new_value=105.)
#shoulde return KeyError: 'The time you passed is not currently in the
# Datarecord, you must the value you pass or first create the new time
# coordinate using the add_record method.'
        dr_time.set_data(time=0.,
                         item_id=0.,
                         data_variable='mean_elevation',
                         new_value=105.)
# should return KeyError: 'This datarecord does not hold items'
def dr_time_bad_properties(dr_time):
    with pytest.raises(AttributeError):
        dr_time.number_of_items
        dr_time.item_coordinates

###############################################################################

########### ITEM ONLY #########################################################
def dr_item_bad_add_record(dr_item):
    with pytest.raises(KeyError):
        dr_item.add_record(time=[0.],
                           item_id=[0],
                           new_record={'new_var' : ['new_data']})
#should return KeyError: 'This Datarecord does not record time'
    with pytest.raises(ValueError):
        dr_item.add_record(item_id=[10],
                           new_record={'new_var' : ['new_data']})
#should return ValueError: There is no item with item_id 14, modify
# the value(s) you pass as item_id or create a new item using the method
#add_item.

def dr_item_bad_add_item(dr_item):
    with pytest.raises(KeyError):
        dr_item.add_item(time = [0.],
                         new_item={'grid_element' : np.array(['node', 'link']),
                                   'element_id' : np.array([1,5])},
                         new_item_spec={'size' : (
                                 ['item_id'], np.array([.2,.3]))})
# should return KeyError: This Datarecord does not record time
def dr_item_bad_get_data(dr_item):
    with pytest.raises(KeyError):
        dr_item.get_data(time=0.,
                         item_id=0,
                         data_variable='element_id')
#should return KeyError: 'This Datarecord does not record time.'
        dr_item.get_data(item_id=0,
                         data_variable='bad_variable')
#should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
def dr_item_bad_set_data(dr_item):
    with pytest.raises(KeyError):
        dr_item.set_data(time=0.,
                         item_id=1,
                         data_variable='element_id',
                         new_value=2)
#should return KeyError: 'This Datarecord does not record time.'
        dr_item.get_data(time=0.,
                         item_id=0,
                         data_variable='bad_variable')
#should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
    with pytest.raises(IndexError):
        dr_item.set_data(item_id=3,
                     data_variable='element_id',
                     new_value=2)
#should return IndexError: The item_id you passed does not exist in this
# Datarecord
def dr_item_bad_properties(dr_item):
    with pytest.raises(AttributeError):
        dr_item.number_of_timesteps
        dr_item.time_coordinates
        dr_item.earliest_time
        dr_item.latest_time
        dr_item.prior_time
###############################################################################

########### ITEM AND TIME #####################################################
def dr_2dim_bad_add_record(dr_2dim):
    with pytest.raises(TypeError):
        dr_2dim.add_record(time=10.,
                       item_id=[1],
                       new_item_loc={'grid_element' : np.array([['cell']]),
                                     'element_id' : np.array([[0]])})
#should return TypeError: You have passed a time that is not permitted,
#must be list or array.
    with pytest.raises(ValueError):
        dr_2dim.add_record(time=[10.],
                       item_id=[10],
                       new_item_loc={'grid_element' : np.array([['cell']]),
                                     'element_id' : np.array([[0]])})
#should return ValueError: There is no item with item_id 10, modify the
#value(s) you pass as item_id or create a new item using the method add_item.
    with pytest.raises(KeyError):
        dr_2dim.add_record(time=[10.],
                       item_id=[0],
                       new_item_loc={'grid_element' : np.array([['cell']])})
#should return KeyError: 'You must provide a new_item_loc dictionnary with
# both grid_element and element_id'
def dr_2dim_bad_add_item(dr_2dim):
    with pytest.raises(ValueError):
        dr_2dim.add_item(new_item={'grid_element' : np.array(
                                             [['node'], ['cell']]),
                               'element_id' : np.array([[2],[0]])},
                     new_item_spec={'size': (['item_id'], [10,5])})
#should return ValueError: The items previously defined in this Datarecord have
# dimensions "time" and "item_id", please provide a "time" for the
#new item(s)
    with pytest.raises(AttributeError):
        dr_2dim.add_item(time=[10.],
                         new_item=('not a dict'),
                         new_item_spec={'size': (['item_id'], [10,5])})
#should return AttributeError: You must provide an new_item dictionnary (see
# documentation for required format)
    with pytest.raises(KeyError):
        dr_2dim.add_item(time=[10.],
                     new_item={'grid_element' : np.array(
                                             [['node'], ['cell']])},
                     new_item_spec={'size': (['item_id'], [10,5])})
#should return AttributeError: You must provide an new_item dictionnary (see
# documentation for required format)
    with pytest.raises(TypeError):
        dr_2dim.add_item(time=10.,
                         new_item={'grid_element' : np.array(
                                                 [['node'], ['cell']]),
                                   'element_id' : np.array([[2],[0]])},
                         new_item_spec={'size': (['item_id'], [10,5])})
#should return TypeError: You have passed a time that is not permitted,
#must be list or array.
def dr_2dim_bad_get_data(dr_2dim):
    with pytest.raises(KeyError):
        dr_2dim.get_data(time=0.,
                            item_id=1,
                            data_variable='bad_variable')
        dr_2dim.get_data(data_variable='bad_variable')
#should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
    with pytest.raises(IndexError):
        dr_2dim.get_data(time=0.,
                            item_id=10,
                            data_variable='grid_element')
#should return IndexError: The item_id you passed does not exist in this
# Datarecord
def dr_2dim_bad_set_data(dr_2dim):
    with pytest.raises(KeyError):
        dr_2dim.set_data(time=0.,
                     item_id=1,
                     data_variable='bad_variable',
                     new_value='node')
    with pytest.raises(ValueError):
        dr_2dim.set_data(time=0.,
                     item_id=1,
                     data_variable='grid_element',
                     new_value='cell')
        dr_2dim.set_data(time=0.,
                     item_id=-1,
                     data_variable='grid_element',
                     new_value='node')
        dr_2dim.set_data(time=0.,
                     item_id=1.,
                     data_variable='grid_element',
                     new_value='node')
    with pytest.raises(IndexError):
        dr_2dim.set_data(time=130.,
                     item_id=1,
                     data_variable='grid_element',
                     new_value='node')
