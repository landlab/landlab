cimport cython
from libc.math cimport exp
from libc.math cimport isinf
from libc.math cimport log

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


def _sequential_ero_depo(
    const id_t [:] stack_flip_ud_sel,
    const id_t [:] flow_receivers,
    const cython.floating [:] cell_area,
    const cython.floating [:] q,
    cython.floating [:] qs,
    cython.floating [:] qs_in,
    const cython.floating [:] Es,
    const cython.floating [:] Er,
    const cython.floating [:] Q_to_the_m,
    const cython.floating [:] slope,
    cython.floating [:] H,
    cython.floating [:] br,
    const cython.floating [:] sed_erosion_term,
    const cython.floating [:] bed_erosion_term,
    const cython.floating [:] K_sed,
    cython.floating [:] ero_sed_effective,
    cython.floating [:] depo_effective,
    const double v,
    const double phi,
    const double F_f,
    const double H_star,
    const double dt,
    const double thickness_lim,
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
                H_loc = H_star * log(
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
            if isinf(H_loc):
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

        # Update deposition rate based on adjusted fluxes
        Hd = H_loc - H_Before
        depo_effective[node_id] = (v*qs_out_adj/q[node_id])/(1 - phi)
        # Deposition should be larger or equal to increase in soil depth
        depo_effective[node_id] = max(depo_effective[node_id], Hd/dt)
        ero_sed_effective[node_id] = depo_effective[node_id] - Hd/dt
    return vol_SSY_riv
