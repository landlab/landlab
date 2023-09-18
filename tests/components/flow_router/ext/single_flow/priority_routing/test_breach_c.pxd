# distutils: language = c++
import numpy as np

cimport cython
cimport numpy as cnp
from libcpp cimport bool
from libcpp.pair cimport pair


cdef extern from "_priority_queue.hpp" nogil:
    cdef cppclass _priority_queue:
        _priority_queue(...) except +
        void push(pair[cnp.int_t, cnp.float_t])
        pair[cnp.int_t, cnp.float_t] top() except +
        void pop()
        bool empty()
        cnp.int_t size()
