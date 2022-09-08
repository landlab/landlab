#distutils: language = c++
#distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
""" Contains a tool function necessary to init the FlowRouter component.
"""

import numpy as np
cimport numpy as cnp
cimport cython

def _get_start_end_indexes_in_sorted_array(cnp.int_t [:] sorted_array, 
    cnp.int_t n_values, cnp.int_t max_value):
    """ Get the start and end indexes of each value in a nd.array sorted by the value.
    Starts are in return[0, :] and ends in return[1, :].
    
    Parameters
    ----------
    sorted_array: memoryview(long).
        Array ordered by its increasing values
    n_values: long
        Number of unique values in the array
    max_value: int
        Maximal value + 1 possible in the array
    
    Return
    ------
    void
    """
    cdef:
        cnp.int_t [:, :] idx2
        cnp.int_t i, value 
    idx2 = np.empty((2, n_values), dtype=int)
    idx2[1, :] = -max_value
    idx2[0, :] = max_value
    for i in range(len(sorted_array)):
            value = sorted_array[i]
            if i < idx2[0, value]: idx2[0, value] = i #min
            if i > idx2[1, value]: idx2[1, value] = i #max
    return idx2
