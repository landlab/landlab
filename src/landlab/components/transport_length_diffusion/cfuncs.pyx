cimport numpy as np
cimport cython
import numpy as np

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t

@cython.boundscheck(False)
cpdef non_local_Depo(DTYPE_FLOAT_t dx,
                     np.ndarray[DTYPE_INT_t, ndim=1] rev_stack,
                     np.ndarray[DTYPE_INT_t, ndim=1] r,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] qs_out,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] L,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] ero,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] depo):
    """
    
    """
    # define internal variables    
    cdef int node
    
    for node in rev_stack:
        depo[node] = qs_out[node]/L[node]     
        qs_out[r[node]] += qs_out[node]+ (ero[node] -depo[node])*dx
                
   
   
@cython.boundscheck(False)
cpdef depo_loop(DTYPE_INT_t nb,                      
                     np.ndarray[DTYPE_INT_t, ndim=1] cores,
                     np.ndarray[DTYPE_INT_t, ndim=1] r,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] flux_in,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] flux_out):
    """
    
    """
    # define internal variables    
    cdef int node
    
    for node in cores:
        flux_in[r[node]] +=  flux_out[node]      