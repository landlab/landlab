"""
Implements the Wiclock and Crowe 2003 bed load transport equations to calculate bed load
rates at links.

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager


Examples
--------
This is the same base example described extensively in river bed dynamics, so
we removed comments that are already available in the main component

>>> import numpy as np
>>> from landlab import RasterModelGrid, imshow_grid
>>> from landlab.components import RiverBedDynamics

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
...     bedload_equation="WilcockAndCrowe",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.000531, -0.001592,  0.000531,  0.      ],
       [ 0.      ,  0.      ,  0.00053 ,  0.      ,  0.      ],
       [ 0.      ,  0.      , -0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Let's take a look to an example in which the gsd contains sand.

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> gsd = [[32, 100, 100], [16, 60, 50], [8, 30, 20], [2, 20, 10], [1, 0, 0]]
>>> sed_transp__bedload_rate_fix_link = np.full(grid.number_of_links, 0.0)
>>> sed_transp__bedload_rate_fix_link[1] = 0.001
>>> sed_transp__bedload_gsd_fix_link = np.full([grid.number_of_links, 4], 0.0)
>>> sed_transp__bedload_gsd_fix_link[1] = [0.20, 0.30, 0.30, 0.20]
>>> sed_transp__bedload_gsd_fix_link[3] = [0.20, 0.30, 0.30, 0.20]
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="WilcockAndCrowe",
...     bed_surf__gsd_loc_node=gsd_loc,
...     sed_transp__bedload_rate_fix_link=sed_transp__bedload_rate_fix_link,
...     sed_transp__bedload_gsd_fix_link=sed_transp__bedload_gsd_fix_link,
... )
>>> rbd.run_one_step()

The gsd retains the sand because it is used in the hiding function in the
Wilcock and Crowe Equation

>>> rbd._gsd
array([[ 32, 100, 100],
       [ 16,  60,  50],
       [  8,  30,  20],
       [  2,  20,  10],
       [  1,   0,   0]])

and the ned bed load at nodes is:

>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)

Here we only show nodes different than zero bed load

>>> print(format(qb[1, 1], ".6f"))
0.001313

>>> print(format(qb[1, 2], ".6f"))
-0.003979

>>> print(format(qb[1, 3], ".6f"))
0.001313

>>> print(format(qb[2, 2], ".6f"))
0.001286

>>> print(format(qb[3, 1], ".6f"))
0.000027

>>> print(format(qb[3, 3], ".6f"))
0.000027

"""

import numpy as np


def bedload_equation(self):
    """Surface-based bedload transport equation of Wilcock and Crowe 2003

    Wilcock, P. R., & Crowe, J. C. (2003). Surface-based transport model
    for mixed-size sediment. Journal of hydraulic engineering, 129(2),
    120-128.
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
    bedload_rate_fix_link = self._sed_transp__bedload_rate_fix_link
    bedload_rate_fix_link_id = self._sed_transp__bedload_rate_fix_link_id
    bedload_gsd_fix_link = self._sed_transp__bedload_gsd_fix_link
    bedload_gsd_fix_link_id = self._sed_transp__bedload_gsd_fix_link_id
    sand_fraction = self._bed_surf__sand_fract_link

    shear_stress_star_sg = shear_stress / (rho * R * g * (gs_geom_mean / 1000))
    shear_stress_star_rsg0 = 0.021 + 0.015 * np.exp(-20 * sand_fraction)
    phi_sg0 = shear_stress_star_sg / shear_stress_star_rsg0

    phi_i = np.zeros((phi_sg0.shape[0], gs.shape[0]))

    qb_G = np.zeros((phi_sg0.shape[0], gs.shape[0]))
    p = np.zeros((phi_sg0.shape[0], gs.shape[0]))

    for i in np.arange(0, gs.shape[0]):
        b = 0.67 / (1 + np.exp(1.5 - gs[i] / gs_geom_mean))
        phi_i[:, i] = phi_sg0 * (gs[i] / (gs_geom_mean)) ** (-b)

    # There are two intervals where qb_G is evaluated
    (id0, id1) = np.where(phi_i >= 1.35)
    if id0.shape[0] > 0:
        qb_G[id0, id1] = 14 * (1 - 0.894 / (phi_i[id0, id1] ** 0.5)) ** 4.5

    (id0, id1) = np.where(phi_i < 1.35)
    if id0.shape[0] > 0:
        qb_G[id0, id1] = 0.002 * phi_i[id0, id1] ** 7.5

    qb_W_star_i = gs_freq * qb_G
    qb_W_star = np.sum(qb_W_star_i, axis=1)

    # Total bedload transport rate is calculated at each link
    # Notice that qb_W_star includes fi (Eq 2 in the paper)
    sed_transp__bedload_rate_link = (
        ((np.sqrt(shear_stress / rho)) ** 3 * qb_W_star) / (R * g)
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
    p = qb_W_star_i / np.reshape(qb_W_star, [qb_W_star.shape[0], 1])
    id0 = np.isnan(p)
    p[id0] = 0
    sed_transp__bedload_gsd_link = p

    # Now that bedload GSD has been calculated in all links we replace those that
    # were imposed. If there is no imposed bedload GSD this will do nothing.
    if bedload_gsd_fix_link_id.size > 0:
        sed_transp__bedload_gsd_link[bedload_gsd_fix_link_id, :] = bedload_gsd_fix_link[
            bedload_gsd_fix_link_id, :
        ]

    return sed_transp__bedload_rate_link, sed_transp__bedload_gsd_link
