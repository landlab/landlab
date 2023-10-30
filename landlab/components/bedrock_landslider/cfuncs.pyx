cimport cython
cimport numpy as np

import numpy as np

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t


@cython.boundscheck(False)
cpdef _landslide_runout(
    DTYPE_FLOAT_t dx,
    DTYPE_FLOAT_t phi,
    DTYPE_FLOAT_t min_deposition_slope,
    np.ndarray[DTYPE_INT_t, ndim=1] stack_rev_sel,
    np.ndarray[DTYPE_INT_t, ndim=2] receivers,
    np.ndarray[DTYPE_FLOAT_t, ndim=2] fract,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Qs_in,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] L_Hill,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Qs_out,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] dH_Hill,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] H_i_temp,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] max_D,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] length_adjacent_cells,
):
    """
    Calculate landslide runout using a non-local deposition algorithm, see:

    * Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.
      (2020) HyLands 1.0: a hybrid landscape evolution model to simulate
      the impact of landslides and landslide-derived sediment on landscape
      evolution. Geosci Model Dev: 13(9):3863–86.
    """
    # define internal variables
    cdef int donor, rcvr, r
    cdef double proportion, dH

    # Iterate backward through the stack, which means we work from upstream to
    # downstream.
    for donor in stack_rev_sel:
        dH = max(
                0,
                min(((Qs_in[donor] / dx) / L_Hill[donor]) / (1 - phi), max_D[donor])
            )
        dH_Hill[donor] += dH
        H_i_temp[donor] += dH

        Qs_in[donor] -= dH * dx * dx * (1 - phi)
        Qs_out[donor] += Qs_in[donor]

        for r in range(receivers.shape[1]):
            rcvr = receivers[donor, r]

            max_D_angle = (
                H_i_temp[donor]
                - min_deposition_slope * length_adjacent_cells[r]
                - H_i_temp[rcvr]
            )
            max_D[rcvr] = min(
                max(max_D[rcvr], H_i_temp[donor] - H_i_temp[rcvr]), max_D_angle
            )

            proportion = fract[donor, r]
            if proportion > 0. and donor != rcvr:
                Qs_in[rcvr] += Qs_out[donor] * proportion
                Qs_in[donor] -= Qs_out[donor] * proportion
