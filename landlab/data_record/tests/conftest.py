import pytest
import numpy as np
from landlab import RasterModelGrid
from landlab.data_record import DataRecord

grid = RasterModelGrid((3,3))

@pytest.fixture
def dr_time():
    # Dimension = time
    time=[0.]
    data_vars={'mean_elevation' : (['time'], np.array([100]))}
    attrs={'time_units' : 'y'}
    return DataRecord(grid=grid,
                      time=time,
                      data_vars=data_vars,
                      attrs=attrs)


@pytest.fixture
def dr_item():
    my_items2 = {'grid_element': np.array(('node', 'link'), dtype=str),
                 'element_id': np.array([1, 3])}
    return DataRecord(grid=grid,
                      items=my_items2)

@pytest.fixture
def dr_2dim():
    time=[0.]
    my_items3 = {'grid_element':np.array([['node'], ['link']]),
                 'element_id': np.array([[1],[3]])}
    my_data_vars = {'mean_elevation' : (['time'], [110.]),
                'item_size' : (['item_id', 'time'],
                               np.array([[0.3], [0.4]]))}
    return DataRecord(grid=grid,
                      time=time,
                      items=my_items3,
                      data_vars=my_data_vars)


# test time=[1,3,1]

#
######### To test errors:
#@pytest.fixture
#def dr_bad_dim():
#    time=[0.]
#    data_vars={'mean_elev' : (['time'], [100.]), 'test' : (['bad_dim'], [12])}
#    return DataRecord(grid,
#                      time=time,
#                      data_vars=data_vars)
## should return ValueError('Data variable dimensions must be time and/or'
##                              'item_id')
#
#@pytest.fixture
#def dr_bad_datavars():
#    time=[0.]
#    data_vars=['not a dict']
#    return DataRecord(grid, time=time, data_vars=data_vars)
## should return TypeError(('Data variables (data_vars) passed to'
##                                 ' DataRecord must be a dictionary (see '
##                                 'documentation for valid structure)'))
#
#@pytest.fixture
#def dr_bad_items():
#    time=[0.]
#    data_vars={'mean_elevation' : (['time'], np.array([100]))}
#    attrs={'time_units' : 'y'}
#    my_items_bad = ['not a dict']
#    return DataRecord(grid,
#                      time=time,
#                      items=my_items_bad,
#                      data_vars=data_vars,
#                      attrs=attrs)
##Should return KeyError(('You must provide an ''items'' dictionary '
##                                 '(see documentation for required format)'))
#
#def dr_bad_items_keys():
#    time=[0.]
#    data_vars={'mean_elevation' : (['time'], np.array([100]))}
#    attrs={'time_units' : 'y'}
#    my_items_bad = {'grid_element':np.array([['node'], ['link']]),
#                    'bad_key': np.array([[1],[3]])}
#    return DataRecord(grid,
#                      time=time,
#                      items=my_items_bad,
#                      data_vars=data_vars,
#                      attrs=attrs)
##Should return KeyError(('You must provide an ''items'' dictionary '
##                                 '(see documentation for required format)'))
#
#
#@pytest.fixture
#def dr_bad_attrs():
#    time=[0.]
#    data_vars={'mean_elevation' : (['time'], np.array([100]))}
#    attrs=['not a dict']
#    my_items3 = {'grid_element':np.array([['node'], ['link']]),
#                 'element_id': np.array([[1],[3]])}
#    return DataRecord(grid,
#                      time=time,
#                      items=my_items3,
#                      data_vars=data_vars,
#                      attrs=attrs)
##Should return except AttributeError:
##                raise TypeError(('Attributes (attrs) passed to DataRecord'
##                                'must be a dictionary'))
#
#
#
#
