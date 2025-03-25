import numpy as np
from cython.parallel import prange

cimport cython
cimport numpy as np

DTYPE_INT = int
ctypedef np.intp_t DTYPE_INT_t

DTYPE = np.double
ctypedef np.double_t DTYPE_t

DTYPE_FLOAT = np.double
DTYPE_complex = np.complexfloating
ctypedef np.double_t DTYPE_FLOAT_t


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
            sed_influxes[r, i ] += sed_outfluxes[c, i]

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
    np.ndarray[DTYPE_t, ndim=2] fractions_from_plucking,
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
            dHdt_by_class[c, i] = porosity_factor * (
                (sed_influxes[c, i ] - sed_outfluxes[c, i])
                / area_of_cell
                + (pluck_rate[c] * fractions_from_plucking[c, i])
                - sed_abr_rates[c, i]
            )
            dHdt[c] += dHdt_by_class[c, i]



@cython.boundscheck(False)
def _calc_sed_abrs_rate(
    DTYPE_INT_t num_sed_classes,
    DTYPE_INT_t num_core_nodes,
    DTYPE_t flow_link_length_over_cell_area,
    np.ndarray[DTYPE_INT_t, ndim=1] core_nodes,
    np.ndarray[DTYPE_t, ndim=2] sed_abr_rates,
    np.ndarray[DTYPE_t, ndim=2] sed_abr_coeff,
    np.ndarray[DTYPE_t, ndim=2] sed_influxes,
    np.ndarray[DTYPE_t, ndim=2] sed_outfluxes,
):
    cdef int c, i
    cdef float flow_link_length
    flow_link_length=flow_link_length_over_cell_area

    for j in range(num_core_nodes):
        c = core_nodes[j]
        for i in range(num_sed_classes):
            sed_abr_rates[c,i] = (sed_abr_coeff[c,i] *
                                  0.5 *
                                  (sed_outfluxes[c,i] + sed_influxes[c,i]) *
                                  flow_link_length)



@cython.boundscheck(False)
def _calc_bedrock_abrs_rate(
    DTYPE_INT_t num_sed_classes,
    DTYPE_INT_t num_core_nodes,
    DTYPE_t flow_link_length_over_cell_area,
    np.ndarray[DTYPE_t, ndim=1] rock_exposure_fraction,
    np.ndarray[DTYPE_INT_t, ndim=1] core_nodes,
    np.ndarray[DTYPE_t, ndim=1] bedrock_abr_rates,
    np.ndarray[DTYPE_t, ndim=1] bedrock_abr_coeff,
    np.ndarray[DTYPE_t, ndim=2] sed_influxes,
    np.ndarray[DTYPE_t, ndim=2] sed_outfluxes,
):
    cdef int c, i
    cdef float flow_link_length
    flow_link_length = flow_link_length_over_cell_area
    for j in range(num_core_nodes):
        c = core_nodes[j]
        for i in range(num_sed_classes):
            bedrock_abr_rates[c] = (bedrock_abr_rates[c] +
                                    (bedrock_abr_coeff[c] *
                                     0.5 *
                                     rock_exposure_fraction[c] *
                                     (sed_outfluxes[c,i] + sed_influxes[c,i]) *
                                     flow_link_length))




@cython.boundscheck(False)
@cython.wraparound(False)
def _get_classes_fractions(
        DTYPE_INT_t num_classes,
        DTYPE_INT_t num_core_nodes,
        np.ndarray[DTYPE_INT_t, ndim=1] core_nodes,
        cython.floating[:, :] value_at_node_per_class,
        cython.floating[:, :] out,
):

    cdef int col, row, node
    cdef float sum_at_node

    for row in prange(num_core_nodes, nogil=True, schedule="static",num_threads=32):
        node = core_nodes[row]
        sum_at_node = 0.0
        for col in range(num_classes):
            sum_at_node  = sum_at_node + value_at_node_per_class[node, col]
        for col in range(num_classes):
            if sum_at_node>0:
                out[node,col] = value_at_node_per_class[node, col] / sum_at_node
            else:
                out[node, col] = 0

    return out.base




@cython.boundscheck(False)
@cython.wraparound(False)
def _calc_pluck_rate(
        DTYPE_INT_t num_classes,
        DTYPE_INT_t num_core_nodes,
        intermittency_factor,
        flow_link_length_over_cell_area,
        np.ndarray[DTYPE_INT_t, ndim = 1] core_nodes,
        cython.floating[:] plucking_coef,
        cython.floating[:] channel_width,
        cython.floating[:] rock_exposure_fraction,
        cython.floating[:, :] classes_fractions,
        cython.floating[:, :] excess_stress,
        cython.floating[:] out,
):


    cdef int col, node, row
    cdef float sum_at_node
    cdef float intermittency_f
    cdef float flow_length

    intermittency_f=intermittency_factor
    flow_length = flow_link_length_over_cell_area
    for row in prange(num_core_nodes, nogil=True, schedule="static",num_threads=32):
        node = core_nodes[row]
        out[node] = 0.0
        for col in range(num_classes):
            out[node]  = out[node] + (intermittency_f*
                                      plucking_coef[node]*
                                      excess_stress[node,col]**1.5*
                                      channel_width[node]*
                                      flow_length*
                                      rock_exposure_fraction[node]*classes_fractions[node,col])

