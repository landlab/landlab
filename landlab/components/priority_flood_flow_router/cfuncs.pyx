import numpy as np
cimport numpy as np
cimport cython


DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


cpdef D8_flowDir(np.ndarray[DTYPE_INT_t, ndim=1] receivers,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] distance_receiver,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] steepest_slope,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] el,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] el_ori,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] dist,
             np.ndarray[DTYPE_INT_t, ndim=1] ngb,
             np.ndarray[DTYPE_INT_t, ndim=1] activeCores,
             np.ndarray[DTYPE_INT_t, ndim=1] activeCells,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] el_d,
             DTYPE_INT_t r, DTYPE_INT_t c,
             DTYPE_FLOAT_t dx,
             np.ndarray[DTYPE_INT_t, ndim=2] adj_link ,
             np.ndarray[DTYPE_INT_t, ndim=1] rec_link):
    
    """
    Calcualte D8 flow dirs
    """
    cdef int idx, i
    
    # for i in range(0,r*c):      
    for i in activeCores:    
        ngb[0] = i + 1
        ngb[1] = i + c 
        ngb[2] = i - 1
        ngb[3] = i -c
        ngb[4] = i + c + 1
        ngb[5] = i + c - 1 
        ngb[6] = i -c - 1
        ngb[7] = i -c + 1
            
        # Differences after filling can be very small, *1e3 to exaggerate those
        # Set to -1 at boundaries (active cell ==0)
        el_d[0] = (el[i] - el[i + 1])*1e3*activeCells[i + 1]-1+activeCells[i + 1] 
        el_d[1] = (el[i] - el[i + c])*1e3 *activeCells[i + c]-1+activeCells[i + c]
        el_d[2] = (el[i] - el[i -1])*1e3  *activeCells[i -1]-1+activeCells[i -1]
        el_d[3] = (el[i] - el[i - c])*1e3 *activeCells[i - c]-1+activeCells[i - c]
        el_d[4] = (el[i] - el[i + c + 1])*1e3/(np.sqrt(2)) *activeCells[i + c + 1]-1+activeCells[i + c + 1]
        el_d[5] = (el[i] - el[i + c - 1])*1e3/(np.sqrt(2)) *activeCells[i + c - 1]-1+activeCells[i + c - 1]
        el_d[6] = (el[i] - el[i -c - 1])*1e3/(np.sqrt(2)) *activeCells[i -c - 1]-1+activeCells[i -c - 1]
        el_d[7] = (el[i] - el[i -c + 1])*1e3/(np.sqrt(2)) *activeCells[i -c + 1]-1+activeCells[i -c + 1]
      
     
        # check to see if pixel is not a lake, if a lake, drain to itself
        if np.max(el_d)>=0:       
            idx = np.argmax(el_d)
            receivers[i] = ngb[idx]
            distance_receiver[i]=dist[idx]
            rec_link[i]= adj_link[i][idx]
            # Slope over original dem can have negative values, but we set it to zero here in analogy to the other Landlab LakeFiller
            steepest_slope[i]  = np.maximum(0,(el_ori[i]-el_ori[receivers[i]])/distance_receiver[i])
        else:
            receivers[i] = i
            rec_link[i]= -1
            steepest_slope[i] = 0
            
            
          
    



cpdef _FA_D8(DTYPE_INT_t np,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] a,
             np.ndarray[DTYPE_FLOAT_t, ndim=1] q,
             np.ndarray[DTYPE_INT_t, ndim=1] stack,
             np.ndarray[DTYPE_INT_t, ndim=1] receivers):
    """
    Accumulates drainage area and discharge, permitting transmission losses.
    """
    cdef int donor, recvr, i

    # Iterate backward through the list, which means we work from upstream to
    # downstream.
    for i in range (np-1,-1,-1):
            donor = stack[i]
            if receivers[donor] != -1:
                rcvr = receivers[donor]
                a[rcvr] += a[donor]
                q[rcvr] += q[donor]            
    
