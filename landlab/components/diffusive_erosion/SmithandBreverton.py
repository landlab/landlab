# -*- coding: utf-8 -*-
"""
DiffusiveErosion Component

@author: D Litwin
"""


import numpy as np

from landlab import Component
from landlab.grid.base import INACTIVE_LINK


class DiffusiveErosion(Component):
    """
    Solves the Smith and Breverton 1972 landscape evolution equation, which
    combines hillslope and channel erosion in a single term.

    .. math::
        \frac{\partial z}{\partial t} = \nabla [(\kappa + c q_w^n ) \nabla z]

    where :math:`z` is topographic elevation, :math:`\kappa` is the hillslope
    linear diffusivity, :math:`q_w` is the surface water discharge, and c and n
    are constants. A primary advantage of this method is that multiple flow direction
    algorithms can be efficiently used to route the streamflow that generates
    erosion, unlike present implementations of the FastScape method, as in
    Armitage 2019 (ESurf).

    """

    _name = "DiffusiveErosion"

    _info = {
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "sediment__unit_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2/s",
            "mapping": "link",
            "doc": "Volume flux per unit width along links",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of land surface in link direction",
        },

    def __init__(
        self,
        grid,
        kappa = 0.01,
        c=1e-4,
        n=1.5,
    ):
