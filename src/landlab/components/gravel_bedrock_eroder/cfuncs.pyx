cimport cython

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
def _calc_sediment_influx(
    const int num_sed_classes,
    const int num_core_nodes,
    cython.floating [:] sediment_influx,
    cython.floating [:, :] sed_influxes,
    const cython.floating [:] sediment_outflux,
    const cython.floating [:, :] sed_outfluxes,
    const id_t [:] core_nodes,
    const id_t [:] receiver_node,
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
    const double upper_limit_dt,
    const int num_nodes,
    const cython.floating [:] sed_thickness,
    const cython.floating [:] elev,
    const cython.floating [:] dHdt,
    const cython.floating [:] dzdt,
    const id_t [:] receiver_node,
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
    const int num_sed_classes,
    const int num_core_nodes,
    const double porosity_factor,
    const double area_of_cell,
    const id_t [:] core_nodes,
    const cython.floating [:, :] pluck_coarse_frac,
    cython.floating [:] dHdt,
    const cython.floating [:] pluck_rate,
    cython.floating [:, :] dHdt_by_class,
    const cython.floating [:, :] sed_influxes,
    const cython.floating [:, :] sed_outfluxes,
    const cython.floating [:, :] sed_abr_rates,
):
    cdef int c, i, j

    for j in range(num_core_nodes):
        c = core_nodes[j]
        dHdt[c] = 0.0
        for i in range(num_sed_classes):
            dHdt_by_class[i, c] = porosity_factor * (
                (sed_influxes[i, c] - sed_outfluxes[i, c])
                / area_of_cell
                + (pluck_rate[c] * pluck_coarse_frac[i, c])
                - sed_abr_rates[i, c]
            )
            dHdt[c] += dHdt_by_class[i, c]
