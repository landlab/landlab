import numpy as np

cimport numpy as np


cdef extern from "math.h":
    double exp(double x) nogil

from libc.math cimport exp
from libc.math cimport log

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t

DTYPE_UINT8 = np.uint8
ctypedef np.uint8_t DTYPE_UINT8_t


def _sequential_ero_depo(
    np.ndarray[DTYPE_INT_t, ndim=1] stack_flip_ud_sel,
    np.ndarray[DTYPE_INT_t, ndim=1] flow_receivers,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] cell_area,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] q,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] qs,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] qs_in,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Es,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Er,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Q_to_the_m,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] slope,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] H,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] br,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] sed_erosion_term,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] bed_erosion_term,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] K_sed,
    DTYPE_FLOAT_t v,
    DTYPE_FLOAT_t phi,
    DTYPE_FLOAT_t F_f,
    DTYPE_FLOAT_t H_star,
    DTYPE_FLOAT_t dt,
    DTYPE_FLOAT_t thickness_lim,
):

    """Calculate and qs and qs_in."""
    # define internal variables
    cdef unsigned int node_id
    cdef double H_Before
    cdef double vol_SSY_riv
    vol_SSY_riv =0.0

    for node_id in stack_flip_ud_sel:
        qs_out = (
            qs_in[node_id]
            + Es[node_id] * cell_area[node_id]
            + (1.0 - F_f) * Er[node_id] * cell_area[node_id]
        ) / (1.0 + (v * cell_area[node_id] / q[node_id]))

        depo_rate = v*qs_out/q[node_id]
        H_loc = H[node_id]
        H_Before = H[node_id]
        slope_loc = slope[node_id]
        sed_erosion_loc = sed_erosion_term[node_id]
        bed_erosion_loc = bed_erosion_term[node_id]

        # Correct for thick soils where soil thickness can grow to inf
        if (H_loc > thickness_lim or slope_loc <= 0 or   sed_erosion_loc==0):
            H_loc += (depo_rate / (1 - phi) - sed_erosion_loc/ (1 - phi)) * dt
        else:
            # Blowup
            if (depo_rate == (K_sed[node_id] * Q_to_the_m[node_id] * slope_loc)) :
                H_loc = H_loc * log(
                    ((sed_erosion_loc / (1 - phi)) / H_star) * dt + exp(H_loc / H_star)
                )
            # No blowup
            else:
                H_loc = H_star * np.log(
                    (1 / ((depo_rate / (1 - phi)) / (sed_erosion_loc / (1 - phi)) - 1))
                    * (
                        exp(
                            (depo_rate / (1 - phi) - (sed_erosion_loc / (1 - phi)))
                            * (dt / H_star)
                        )
                        * (
                            (
                                (depo_rate / (1 - phi) / (sed_erosion_loc / (1 - phi)))
                                - 1
                            )
                            * exp(H_loc / H_star)
                            + 1
                        )
                        - 1
                    )
                )
            # In case soil depth evolves to infinity, fall back to no entrainment
            if H_loc == np.inf:
                H_loc = (
                    H[node_id]
                    + (depo_rate / (1 - phi) - sed_erosion_loc/ (1 - phi)) * dt
                )

        H_loc = max(0, H_loc)
        ero_bed = bed_erosion_loc* (exp(-H_loc / H_star))

        # should always be bigger than 0
        qs_out_adj = (
            qs_in[node_id]
            - ((H_loc - H_Before) * (1 - phi) * cell_area[node_id] / dt)
            + (1.0 - F_f) * ero_bed * cell_area[node_id]
        )

        qs[node_id] = qs_out_adj
        qs_in[flow_receivers[node_id]] += qs[node_id]

        H[node_id] = H_loc
        br[node_id] += -dt * ero_bed
        vol_SSY_riv += F_f*ero_bed* cell_area[node_id]

    return vol_SSY_riv
