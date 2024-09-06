"""
Implements several bed load transport equations to calculate bed load rates at
links. All equations are of the style qb =  a * (tau* - tau*c) **b

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

Examples
--------
This is the same base example described extensively in river bed dynamics, so
we removed comments that are already available in the main component

>>> import numpy as np
>>> from landlab import RasterModelGrid, imshow_grid
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
>>> z = copy.deepcopy(grid.at_node["topographic__elevation"])

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

Case 1, we calculate using the default equation, MPM

>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.002297, -0.006891,  0.002297,  0.      ],
       [ 0.      ,  0.      ,  0.002297,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 1b, MPM with variable_critical_shear_stress=True

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bed_surf__gsd_loc_node=gsd_loc,
...     variable_critical_shear_stress=True,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.00065 , -0.001949,  0.00065 ,  0.      ],
       [ 0.      ,  0.      ,  0.00065 ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 1c, MPM with a bed load rate fixed in one link

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> qb_fix_link = np.zeros([grid.number_of_links, 1])
>>> qb_fix_link[[29]] = 0.01
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bed_surf__gsd_loc_node=gsd_loc,
...     variable_critical_shear_stress=True,
...     sed_transp__bedload_rate_fix_link=qb_fix_link,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.00065 , -0.001949,  0.00065 ,  0.      ],
       [ 0.      ,  0.      ,  0.00065 ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.01    , -0.01    ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 2, we use Fernandez Luque and Van Beek

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="FLvB",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.001682, -0.005047,  0.001682,  0.      ],
       [ 0.      ,  0.      ,  0.001682,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 2b, FLvB with variable_critical_shear_stress=True

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="FLvB",
...     variable_critical_shear_stress=True,
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.000463, -0.001389,  0.000463,  0.      ],
       [ 0.      ,  0.      ,  0.000463,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 3, we use Wong and Parker

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="WongAndParker",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.001133, -0.003398,  0.001133,  0.      ],
       [ 0.      ,  0.      ,  0.001133,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 3b, WongAndParker with variable_critical_shear_stress=True

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="WongAndParker",
...     variable_critical_shear_stress=True,
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.000295, -0.000884,  0.000295,  0.      ],
       [ 0.      ,  0.      ,  0.000295,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 4, we use Huang

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="Huang",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.001188, -0.003564,  0.001188,  0.      ],
       [ 0.      ,  0.      ,  0.001188,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Case 4b, Huang with variable_critical_shear_stress=True

>>> grid.at_node["topographic__elevation"] = z.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="Huang",
...     variable_critical_shear_stress=True,
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.000292, -0.000876,  0.000292,  0.      ],
       [ 0.      ,  0.      ,  0.000292,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])
"""

import numpy as np


def bedload_equation(self):
    # Variables definition
    shear_stress = self._surface_water__shear_stress_link
    shear_stress_signed = self._shear_stress
    rho = self._rho
    R = self._R
    g = self._g
    gs_D50 = self._bed_surf__median_size_link
    dz_ds = self._dz_ds
    selected_bedload_equation = self._bedload_equation
    var_cr_shear_stress = self._variable_critical_shear_stress
    bedload_rate_fix_link = self._sed_transp__bedload_rate_fix_link
    bedload_rate_fix_link_id = self._sed_transp__bedload_rate_fix_link_id

    tau_star = shear_stress / (
        rho * R * g * (gs_D50 / 1000)
    )  # dimensionless shear stress

    # calculates spatially variable dimensionless critical shear stress
    # Asssign 0 where slope is below 0.03 which wil be replaced at each equation
    bed_slope = dz_ds * np.sign(shear_stress_signed)  # slope considers flow direction
    tau_star_cr_0 = np.where(bed_slope > 0.03, 2.18 * bed_slope + 0.021, 0)
    qb_star = np.zeros_like(tau_star)

    if selected_bedload_equation == "MPM":
        qb_star = MeyerPeter_Muller(
            qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0
        )
    elif selected_bedload_equation == "FLvB":
        qb_star = FernandezLuque_VanBeek(
            qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0
        )
    elif selected_bedload_equation == "WongAndParker":
        qb_star = Wong_Parker(qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0)
    elif selected_bedload_equation == "Huang":
        qb_star = Huang(qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0)

    qb_star = np.where(
        np.isnan(qb_star),
        0,
        qb_star,
    )

    qb = (qb_star * (np.sqrt(R * g * (gs_D50 / 1000)) * (gs_D50 / 1000))) * np.sign(
        shear_stress_signed
    )

    # Restores fixed bed load rates at fixed links
    if bedload_rate_fix_link_id.size > 0:
        qb[bedload_rate_fix_link_id] = bedload_rate_fix_link[bedload_rate_fix_link_id]

    return qb  # qb = sed_transp__bedload_rate_link


def MeyerPeter_Muller(qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0):
    """Surface-based bedload transport equation of Meyer-Peter and Müller

    Meyer-Peter, E. and Müller, R., 1948, Formulas for Bed-Load Transport,
    Proceedings, 2nd Congress, International Association of Hydraulic
    Research, Stockholm: 39-64.
    """
    tau_star_cr = np.full(tau_star.shape, 0.047)
    qb_star_coeff = 8
    qb_star_exp = 3 / 2

    if var_cr_shear_stress is True:
        tau_star_cr = np.where(tau_star_cr_0 < 0.021, tau_star_cr, tau_star_cr_0)

    mask = tau_star > tau_star_cr
    qb_star[mask] = qb_star_coeff * (tau_star[mask] - tau_star_cr[mask]) ** (
        qb_star_exp
    )

    return qb_star


def FernandezLuque_VanBeek(qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0):
    """Surface-based bedload transport equation of Fernandez Luque and
    van Beek

    Fernandez Luque, R. and R. van Beek, 1976, Erosion and transport of
    bedload sediment, Journal of Hydraulic Research, 14(2): 127-144.
    """
    tau_star_cr = np.full(tau_star.shape, 0.045)
    qb_star_coeff = 5.7
    qb_star_exp = 3 / 2

    if var_cr_shear_stress is True:
        tau_star_cr = np.where(tau_star_cr_0 < 0.021, tau_star_cr, tau_star_cr_0)

    mask = tau_star > tau_star_cr
    qb_star[mask] = qb_star_coeff * (tau_star[mask] - tau_star_cr[mask]) ** (
        qb_star_exp
    )

    return qb_star


def Wong_Parker(qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0):
    """Surface-based bedload transport equation of Wong and Parker

    Wong and Parker 2006, Reanalysis and Correction of Bed-Load Relation of
    Meyer-Peter and Müller Using Their Own Database. Journal of Hydraulic
    Engineering. Volume 132, Issue 11

    """
    tau_star_cr = np.full(tau_star.shape, 0.047)
    qb_star_coeff = 4.93
    qb_star_exp = 1.6

    if var_cr_shear_stress is True:
        tau_star_cr = np.where(tau_star_cr_0 < 0.021, tau_star_cr, tau_star_cr_0)

    mask = tau_star > tau_star_cr
    qb_star[mask] = qb_star_coeff * (tau_star[mask] - tau_star_cr[mask]) ** (
        qb_star_exp
    )

    return qb_star


def Huang(qb_star, tau_star, var_cr_shear_stress, tau_star_cr_0):
    """Surface-based bedload transport equation of He Qing Huang

    Huang, H. Q. (2010), Reformulation of the bed load equation of Meyer‐Peter
    and Müller in light of the linearity theory for alluvial channel flow,
    Water Resour. Res., 46, W09533, doi:10.1029/2009WR008974.

    """
    tau_star_cr = np.full(tau_star.shape, 0.047)
    qb_star_coeff = 6.0
    qb_star_exp = 5 / 3

    if var_cr_shear_stress is True:
        tau_star_cr = np.where(tau_star_cr_0 < 0.021, tau_star_cr, tau_star_cr_0)

    mask = tau_star > tau_star_cr
    qb_star[mask] = qb_star_coeff * (tau_star[mask] - tau_star_cr[mask]) ** (
        qb_star_exp
    )

    return qb_star
