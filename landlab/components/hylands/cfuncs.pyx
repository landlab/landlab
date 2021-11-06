cimport numpy as np
cimport cython
import numpy as np

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t

@cython.boundscheck(False)
cpdef non_local_Depo(DTYPE_INT_t nb, DTYPE_INT_t q, 
                               DTYPE_FLOAT_t dx,DTYPE_FLOAT_t dx2,DTYPE_FLOAT_t phi,
                               np.ndarray[DTYPE_INT_t, ndim=1] stack,
                               np.ndarray[DTYPE_INT_t, ndim=2] receivers,
                               np.ndarray[DTYPE_FLOAT_t, ndim=2] fract,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] Qs_in,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] L_Hill,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] Qs_out,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] dH_Hill,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] H_i_temp,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] max_D,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] max_dH,
                               np.ndarray[DTYPE_FLOAT_t, ndim=1] Grid_Status):
    """
    
    """
    # define internal variables    
    cdef int donor, rcvr, i, r
    cdef double accum, proportion
    
           
    # Iterate backward through the stack, which means we work from upstream to
    # downstream.
    for i in range(nb-1, -1, -1):
        donor = stack[i]
        # No deposition at boundaries
        if Grid_Status[donor] != 1:
            # Following Carretier 2016
            dH =min(max_dH[donor],
                    max(0,
                        min(((Qs_in[donor]/dx)/L_Hill[donor])/(1-phi),
                        max_D[donor])
                        ))
            dH_Hill[donor] = dH_Hill[donor] + dH
            H_i_temp[donor] = H_i_temp[donor] + dH
            
            Qs_in[donor] = Qs_in[donor] - dH*dx2*(1-phi)
            Qs_out[donor] = Qs_out[donor] + Qs_in[donor]
            
            for r in range(0,q):
                rcvr = receivers[donor,r]
                max_D[rcvr] = max(max_D[rcvr] , H_i_temp[donor] - H_i_temp[rcvr])
                
                proportion = fract[donor,r]
                if proportion > 0.:
                    if donor != rcvr:
                        Qs_in[rcvr] = Qs_in[rcvr] + Qs_out[donor] * proportion
                        Qs_in[donor] =Qs_in[donor] - Qs_out[donor] * proportion
                        
                # Take care of grain sizes
                

