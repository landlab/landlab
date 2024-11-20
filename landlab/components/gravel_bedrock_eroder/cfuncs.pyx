import numpy as np
from cython.parallel import prange

cimport cython
cimport numpy as np

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t

DTYPE = np.double
ctypedef np.double_t DTYPE_t


@cython.boundscheck(False)
def _calc_sediment_influx(
    DTYPE_INT_t num_sed_classes,
    DTYPE_INT_t num_core_nodes,
    np.ndarray[DTYPE_t, ndim=1] sediment_influx,
    np.ndarray[DTYPE_t, ndim=2] sed_influxes,
    np.ndarray[DTYPE_t, ndim=1] sediment_outflux,
    np.ndarray[DTYPE_t, ndim=2] sed_outfluxes,
    np.ndarray[DTYPE_INT_t, ndim=1] core_nodes,
    np.ndarray[DTYPE_INT_t, ndim=1] receiver_node,
):
    cdef int i, c, r

    sediment_influx[:] = 0.0
    sed_influxes[:, :] = 0.0
    for i in range(num_core_nodes):  # send sediment downstream
        c = core_nodes[i]
        r = receiver_node[c]
        if num_sed_classes > 1:
            sediment_influx[r] += sediment_outflux[c]
        for i in range(num_sed_classes):
            sed_influxes[i, r] += sed_outfluxes[i, c]

@cython.boundscheck(False)
def _estimate_max_time_step_size_ext(
    DTYPE_t upper_limit_dt,
    DTYPE_INT_t num_nodes,
    np.ndarray[DTYPE_t, ndim=1] sed_thickness,
    np.ndarray[DTYPE_t, ndim=1] elev,
    np.ndarray[DTYPE_t, ndim=1] dHdt,
    np.ndarray[DTYPE_t, ndim=1] dzdt,
    np.ndarray[DTYPE_INT_t, ndim=1] receiver_node,
):
    cdef int i, r
    cdef float min_time_to_exhaust_sed, min_time_to_flatten_slope, rate_diff, height_above_rcvr

    min_time_to_exhaust_sed = upper_limit_dt
    min_time_to_flatten_slope = upper_limit_dt

    for i in range(num_nodes):
        if dHdt[i] < 0.0 and sed_thickness[i] > 0.0:
            min_time_to_exhaust_sed = min(min_time_to_exhaust_sed, -sed_thickness[i] / dHdt[i])
        r = receiver_node[i]
        rate_diff = dzdt[r] - dzdt[i]
        height_above_rcvr = elev[i] - elev[r]
        if rate_diff > 0.0 and height_above_rcvr > 0.0:
            min_time_to_flatten_slope = min(min_time_to_flatten_slope, height_above_rcvr / rate_diff)

    return 0.5 * min(min_time_to_exhaust_sed, min_time_to_flatten_slope)

@cython.boundscheck(False)
def _calc_sediment_rate_of_change(
    DTYPE_INT_t num_sed_classes,
    DTYPE_INT_t num_core_nodes,
    DTYPE_t porosity_factor,
    DTYPE_t area_of_cell,
    np.ndarray[DTYPE_INT_t, ndim=1] core_nodes,
    np.ndarray[DTYPE_t, ndim=1] pluck_coarse_frac,
    np.ndarray[DTYPE_t, ndim=1] dHdt,
    np.ndarray[DTYPE_t, ndim=1] pluck_rate,
    np.ndarray[DTYPE_t, ndim=2] dHdt_by_class,
    np.ndarray[DTYPE_t, ndim=2] sed_influxes,
    np.ndarray[DTYPE_t, ndim=2] sed_outfluxes,
    np.ndarray[DTYPE_t, ndim=2] sed_abr_rates,
    np.ndarray[DTYPE_t, ndim=2] br_abr_rates,
):
    cdef int c, i, j

    for j in range(num_core_nodes):
        c = core_nodes[j]
        dHdt[c] = 0.0
        for i in range(num_sed_classes):
            dHdt_by_class[i, c] = porosity_factor * (
                (sed_influxes[i, c] - sed_outfluxes[i, c])
                / area_of_cell
                + (pluck_rate[c] * pluck_coarse_frac[i])
                - sed_abr_rates[i, c]
            )
            dHdt[c] += dHdt_by_class[i, c]
