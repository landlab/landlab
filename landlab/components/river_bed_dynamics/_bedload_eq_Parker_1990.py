"""
Implements the Parker 1990 bed load transport equations to calculate bed load rates
at links.

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager


Examples
--------
This is the same base example described extensively in river bed dynamics, so
we removed comments that are already available in the main component

>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import RiverBedDynamics
>>> import copy

>>> grid = RasterModelGrid((5, 5))
>>> grid.at_node["topographic__elevation"] = [
...     [1.07, 1.06, 1.00, 1.06, 1.07],
...     [1.08, 1.07, 1.03, 1.07, 1.08],
...     [1.09, 1.08, 1.07, 1.08, 1.09],
...     [1.09, 1.09, 1.08, 1.09, 1.09],
...     [1.09, 1.09, 1.09, 1.09, 1.09],
... ]
>>> z = grid.at_node["topographic__elevation"].copy()
>>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])
>>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
>>> grid.at_node["surface_water__velocity"] = np.full(grid.number_of_nodes, 0.25)
>>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.102)
>>> grid.at_link["surface_water__velocity"] = np.full(grid.number_of_links, 0.25)

>>> gsd_loc = [
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
... ]

>>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="Parker1990",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.001598, -0.004793,  0.001598,  0.      ],
       [ 0.      ,  0.      ,  0.001597,  0.      ,  0.      ],
       [ 0.      ,  0.      , -0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Let's take a look to an example in which the gsd contains sand.

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> gsd = [[32, 100, 100], [16, 60, 50], [8, 30, 20], [2, 20, 10], [1, 0, 0]]
>>> sed_transp__bedload_rate_fix_link = np.full(grid.number_of_links, 0.0)
>>> sed_transp__bedload_rate_fix_link[1] = 0.001

Notice that there is sand in gsd, therefore after removing sand thre will 3
equivalent grain sizes

>>> sed_transp__bedload_gsd_fix_link = np.full([grid.number_of_links, 3], 0.0)
>>> sed_transp__bedload_gsd_fix_link[1] = [0.20, 0.30, 0.50]
>>> sed_transp__bedload_gsd_fix_link[3] = [0.20, 0.30, 0.50]
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="Parker1990",
...     bed_surf__gsd_loc_node=gsd_loc,
...     sed_transp__bedload_rate_fix_link=sed_transp__bedload_rate_fix_link,
...     sed_transp__bedload_gsd_fix_link=sed_transp__bedload_gsd_fix_link,
... )
>>> rbd.run_one_step()

The sand-free gsd is:

>>> rbd._gsd
array([[ 32, 100, 100],
       [ 16,  50,  44],
       [  8,  12,  11],
       [  2,   0,   0]])

and the ned bed load at nodes is:

>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)

Here we only show nodes different than zero bed load

>>> print(format(qb[1, 1], ".6f"))
0.001473

>>> print(format(qb[1, 2], ".6f"))
-0.004421

>>> print(format(qb[1, 3], ".6f"))
0.001473

>>> print(format(qb[2, 2], ".6f"))
0.001472

>>> print(format(qb[3, 1], ".6f"))
0.000001

>>> print(format(qb[3, 3], ".6f"))
0.000001

"""

import numpy as np
from scipy.interpolate import interp1d


def strain_functions(phi_sgo):
    """Determine the strain functions values for the Parker 1990 bed load equation

    Examples
    --------
    Here we test that the interpolations are correct. We will test for different phi_sgo values.

    >>> import numpy as np
    >>> from . import _bedload_eq_Parker_1990 as Parker1990

    range 1: phi_sgo <= 0.7639

    >>> phi_sgo = np.array([0.1, 0.2, 0.3, 0.7639])
    >>> (omega0, sigma0) = Parker1990.strain_functions(phi_sgo)
    >>> omega0
    array([ 1.011,  1.011,  1.011,  1.011])

    >>> sigma0
    array([ 0.8157,  0.8157,  0.8157,  0.8157])

    range 2: (phi_sgo > 0.7639) & (phi_sgo < 231.2)

    >>> phi_sgo = np.array([0.7639, 1, 100, 231.1])
    >>> (omega0, sigma0) = Parker1990.strain_functions(phi_sgo)
    >>> omega0
    array([ 1.011     ,  0.9997    ,  0.45601387,  0.45409851])

    >>> sigma0
    array([ 0.8157    ,  0.8439    ,  1.49812569,  1.49900013])

    range 3: phi_sgo >= 231.2

    >>> phi_sgo = np.array([231.1, 1000, 2000, 2320])
    >>> (omega0, sigma0) = Parker1990.strain_functions(phi_sgo)
    >>> omega0
    array([ 0.45409851,  0.45358472,  0.45291448,  0.4527    ])

    >>> sigma0
    array([ 1.49900013,  1.49936806,  1.4998468 ,  1.5       ])

    """

    po = np.array(
        [
            [0.6684, 0.7639, 0.8601, 0.9096, 0.9615, 1.000],
            [1.055, 1.108, 1.197, 1.302, 1.407, 1.529],
            [1.641, 1.702, 1.832, 1.937, 2.044, 2.261],
            [2.499, 2.732, 2.993, 3.477, 4.075, 4.469],
            [5.016, 6.158, 7.821, 10.06, 14.38, 19.97],
            [25.79, 38.57, 68.74, 91.95, 231.2, 2320],
        ]
    ).flatten()

    oo = np.array(
        [
            [1.011, 1.011, 1.01, 1.008, 1.004, 0.9997],
            [0.9903, 0.9789, 0.9567, 0.9273, 0.8964, 0.8604],
            [0.8287, 0.8123, 0.7796, 0.7554, 0.7326, 0.6928],
            [0.6585, 0.6345, 0.615, 0.5877, 0.564, 0.5523],
            [0.5395, 0.5209, 0.5045, 0.4917, 0.479, 0.4712],
            [0.4668, 0.462, 0.4578, 0.4564, 0.4541, 0.4527],
        ]
    ).flatten()

    so = np.array(
        [
            [0.8157, 0.8157, 0.8182, 0.8233, 0.8333, 0.8439],
            [0.8621, 0.8825, 0.9214, 0.9723, 1.025, 1.083],
            [1.13, 1.153, 1.196, 1.225, 1.25, 1.287],
            [1.313, 1.333, 1.352, 1.38, 1.403, 1.414],
            [1.426, 1.444, 1.458, 1.469, 1.48, 1.486],
            [1.49, 1.493, 1.497, 1.498, 1.499, 1.5],
        ]
    ).flatten()

    omega0 = np.zeros_like(phi_sgo)
    sigma0 = np.zeros_like(phi_sgo)

    # There are three intervals where we can interpolate Omega and Sigma
    (link_id_phi_sgo,) = np.where(phi_sgo <= 0.7639)
    if link_id_phi_sgo.shape[0] > 0:
        omega0[link_id_phi_sgo] = oo[0]
        sigma0[link_id_phi_sgo] = so[0]

    (link_id_phi_sgo,) = np.where(phi_sgo >= 231.2)
    if link_id_phi_sgo.shape[0] > 0:
        omega0[link_id_phi_sgo] = np.interp(phi_sgo, po, oo)[link_id_phi_sgo]
        sigma0[link_id_phi_sgo] = np.interp(phi_sgo, po, so)[link_id_phi_sgo]

    (link_id_phi_sgo,) = np.where((phi_sgo > 0.7639) & (phi_sgo < 231.2))
    if link_id_phi_sgo.shape[0] > 0:
        foo = interp1d(po, oo, kind="cubic")
        fso = interp1d(po, so, kind="cubic")
        omega0[link_id_phi_sgo] = foo(phi_sgo[link_id_phi_sgo])
        sigma0[link_id_phi_sgo] = fso(phi_sgo[link_id_phi_sgo])

    return omega0, sigma0


def bedload_equation(self):
    """Surface-based bedload transport equation of Parker 1990

    G. Parker (1990) Surface-based bedload transport relation for gravel rivers,
    Journal of Hydraulic Research, 28:4, 417-436, DOI: 10.1080/00221689009499058
    """

    # Variables definition
    shear_stress = self._surface_water__shear_stress_link
    shear_stress_signed = self._shear_stress
    rho = self._rho
    R = self._R
    g = self._g
    gs = self._gs
    gs_freq = self._bed_surf__gsd_link
    gs_geom_mean = self._bed_surf__geom_mean_size_link
    gs_geo_std = self._bed_surf__geo_std_size_link
    bedload_rate_fix_link = self._sed_transp__bedload_rate_fix_link
    bedload_rate_fix_link_id = self._sed_transp__bedload_rate_fix_link_id
    bedload_gsd_fix_link = self._sed_transp__bedload_gsd_fix_link
    bedload_gsd_fix_link_id = self._sed_transp__bedload_gsd_fix_link_id

    tau_star_rsgo = 0.0386  # Reference dimensionless shear
    beta = 0.0951  # stress for the median size

    shear_stress_star_sg = shear_stress / (rho * R * g * (gs_geom_mean / 1000))
    phi_sgo = shear_stress_star_sg / tau_star_rsgo

    (omega0, sigma0) = strain_functions(phi_sgo)
    omega = 1 + (np.log2(gs_geo_std) / sigma0) * (omega0 - 1)

    gs_Di_Dsg = np.tile(gs, (self._grid.number_of_links, 1)) / (
        np.reshape(gs_geom_mean, [self._grid.number_of_links, 1])
    )
    phi_sgo = np.reshape(phi_sgo, [phi_sgo.shape[0], 1])
    omega = np.reshape(omega, [omega.shape[0], 1])
    phi_i = omega * phi_sgo * (gs_Di_Dsg) ** -beta

    qb_G = np.zeros_like(phi_i)

    # There are three intervals where qb_G is evaluated
    (id0, id1) = np.where(phi_i > 1.59)
    if id0.shape[0] > 0:
        qb_G[id0, id1] = 5474 * (1 - 0.853 / phi_i[id0, id1]) ** 4.5

    (id0, id1) = np.where((phi_i >= 1) & (phi_i <= 1.59))
    if id0.shape[0] > 0:
        qb_G[id0, id1] = np.exp(
            14.2 * (phi_i[id0, id1] - 1) - 9.28 * (phi_i[id0, id1] - 1) ** 2
        )

    (id0, id1) = np.where(phi_i < 1)
    if id0.shape[0] > 0:
        qb_G[id0, id1] = phi_i[id0, id1] ** 14.2

    qb_Gf = gs_freq * qb_G
    qb_Gf_sum = np.sum(qb_Gf, axis=1)

    # Total bedload transport rate
    qb_W_star_s = 0.00218 * qb_Gf_sum
    sed_transp__bedload_rate_link = (
        ((np.sqrt(shear_stress / rho)) ** 3 * qb_W_star_s) / (R * g)
    ) * np.sign(shear_stress_signed)

    # Restores fixed bed load rates at fixed links
    if bedload_rate_fix_link_id.size > 0:
        sed_transp__bedload_rate_link[bedload_rate_fix_link_id] = bedload_rate_fix_link[
            bedload_rate_fix_link_id
        ]

    # During the calculation of fractional sediment transport, there might
    # be instances where certain fractions will be zero. This situation
    # could trigger a division by zero warning. However, we are aware that
    # in these specific cases, there's no transport occurring. Thus, such a
    # warning is not indicative of an actual problem in the calculation.
    # To avoid unnecessary warnings, we've deactivated them.

    np.seterr(divide="ignore", invalid="ignore")

    # Frational bedload transport rate
    p = np.zeros_like(qb_Gf)
    qb_Gf_sum = np.transpose(np.tile(qb_Gf_sum, (gs_freq.shape[1], 1)))
    p = qb_Gf / qb_Gf_sum
    id0 = np.isnan(p)
    p[id0] = 0
    sed_transp__bedload_gsd_link = p

    # Now that bedload GSD has been calculated in all links we replace those that
    # are fixed. If there is no imposed bedload GSD this will do nothing.
    if bedload_gsd_fix_link_id.size > 0:
        sed_transp__bedload_gsd_link[bedload_gsd_fix_link_id, :] = bedload_gsd_fix_link[
            bedload_gsd_fix_link_id, :
        ]

    return sed_transp__bedload_rate_link, sed_transp__bedload_gsd_link
