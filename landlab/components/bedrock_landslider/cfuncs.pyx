cimport cython

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
cpdef _landslide_runout(
    double dx,
    double phi,
    double min_deposition_slope,
    id_t [:] stack_rev_sel,
    id_t [:, :] receivers,
    cython.floating [:, :] fract,
    cython.floating [:] Qs_in,
    cython.floating [:] L_Hill,
    cython.floating [:] Qs_out,
    cython.floating [:] dH_Hill,
    cython.floating [:] H_i_temp,
    cython.floating [:] max_D,
    cython.floating [:] length_adjacent_cells,
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
