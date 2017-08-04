#! /usr/bin/env python
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
import pickle 
from nose.tools import assert_equal
import os
from landlab.io.native_landlab import save_grid, load_grid

def compare_dictionaries(dict_1, dict_2, dict_1_name, dict_2_name, path=""):
    """Compare two dictionaries recursively to find non mathcing elements

    Args:
        dict_1: dictionary 1
        dict_2: dictionary 2

    Returns:

    """
    import numpy as np
    err = ''
    key_err = ''
    value_err = ''
    old_path = path
    for k in dict_1.keys():
        path = old_path + "[%s]" % k
        if not k in dict_2:
            key_err += "Key %s%s not in %s\n" % (dict_2_name, path, dict_2_name)
        else:
            if isinstance(dict_1[k], dict) and isinstance(dict_2[k], dict):
                err += compare_dictionaries(dict_1[k],dict_2[k],'d1','d2', path)
            else:
                o1 = dict_1[k]
                o2 = dict_2[k]
                try:
                    if o1  != o2:
                        value_err += "Value of %s%s (%s) not same as %s%s (%s)\n"\
                            % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])
                            
                except ValueError:
                    if np.array_equal(np.asarray(o1), np.asarray(o1)) == False:
                        value_err += "Value of %s%s (%s) not same as %s%s (%s)\n"\
                            % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])
                    
    for k in dict_2.keys():
        path = old_path + "[%s]" % k
        if not k in dict_1:
            key_err += "Key %s%s not in %s\n" % (dict_2_name, path, dict_1_name)

    return key_err + value_err + err


def test_pickle():
    # Make a simple-ish grid
    mg1 = RasterModelGrid(10,10,2.)
    z = mg1.add_zeros('node', 'topographic__elevation')
    z += mg1.node_x.copy()
    fa = FlowAccumulator(mg1, flow_director='D8')
    fa.run_one_step()
    
    # save it with pickle
    with open('testsavedgrid.grid', 'wb') as f:
        pickle.dump(mg1, f)
    
    # load it with pickle
    with open('testsavedgrid.grid', 'rb') as f:
        mg2 = pickle.load(f)
     
    # compare the two
    len(mg1.__dict__) == len(mg2.__dict__)
    mg1keys = sorted(list(mg1.__dict__.keys()))
    mg2keys = sorted(list(mg2.__dict__.keys()))
    
    for i in range(len(mg1keys)):        
        assert_equal(mg1keys[i], mg2keys[i])
        
    a = compare_dictionaries(mg1.__dict__,mg2.__dict__,'m1','m2')
    assert_equal(a, '')
    
    os.remove('testsavedgrid.grid')

def test_save():
    # Make a simple-ish grid
    mg1 = RasterModelGrid(10,10,2.)
    z = mg1.add_zeros('node', 'topographic__elevation')
    z += mg1.node_x.copy()
    fa = FlowAccumulator(mg1, flow_director='D8')
    fa.run_one_step()
    
    save_grid(mg1, 'testsavedgrid.grid')
    
    mg2 = load_grid('testsavedgrid.grid')
    
    # compare the two
    len(mg1.__dict__) == len(mg2.__dict__)
    mg1keys = sorted(list(mg1.__dict__.keys()))
    mg2keys = sorted(list(mg2.__dict__.keys()))
    
    for i in range(len(mg1keys)):        
        assert_equal(mg1keys[i], mg2keys[i])
        
    a = compare_dictionaries(mg1.__dict__,mg2.__dict__,'m1','m2')
    assert_equal(a, '')
    
    os.remove('testsavedgrid.grid')