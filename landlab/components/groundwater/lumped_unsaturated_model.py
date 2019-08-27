# -*- coding: utf-8 -*-
"""
LumpedUnsaturatedZone Component

@author: D Litwin
"""


import numpy as np
from landlab import RasterModelGrid, Component
from landlab.utils import return_array_at_node

def EvapotranspirationLoss(p,S,Sw,St,ETmax,dt):
    qet = np.zeros_like(S)

    for i in range(len(qet)):
        if St[i] > 0:
            if p == 0.0 and S[i] > Sw[i]+(S[i]-Sw[i])/(St[i]-Sw[i])*ETmax[i]*dt:
                qet[i] = (S[i]-Sw[i])/(St[i]-Sw[i])*ETmax[i]
            elif p == 0.0 and (S[i]-Sw[i])/(St[i]-Sw[i])*ETmax[i] > (S[i] - Sw[i])/dt:
                qet[i] = (S[i] - Sw[i])/dt
            else:
                qet[i] = 0
        else:
            qet[i] = 0
    return qet

def LeakageLoss(S,f,qet,Sf,Ksat,dt):
    ql = np.zeros_like(S)

    for i in range(len(ql)):
        if S[i] + f[i]*dt- qet[i]*dt - Ksat[i]*dt >= Sf[i]:
            ql[i] = Ksat[i]
        elif S[i] + f[i]*dt - qet[i]*dt >= Sf[i] and S[i] + f[i]*dt- qet[i]*dt - Ksat[i]*dt < Sf[i]:
            ql[i] = (S[i] + f[i]*dt - qet[i]*dt - Sf[i])/dt
        else:
            ql[i] = 0

    return ql

# def RootWaterLoss():


class LumpedUnsaturatedZone(Component):
    """
    Simulate storage and loss of water from the unsaturated zone above
    a Boussinesq aquifer.

    Parameters
    ----------
    grid: ModelGrid
            Landlab ModelGrid object
    soil_field_capacity: float, field name, or array of float
            the dimensionless soil field capacity, when free
            drainage begins
            Default = 0.3
    soil_wilting_point: float, field name, or array of float
            the dimensionless soil wilting point, when ET no longer
            removes water from the soil
            Default = 0.18
    plant_rooting_depth: float, field name, or array of float
            The depth to which plant roots will draw water from the
            water table. If water table is below this depth, plant root
            water uptake is zero.
            Default = 1 m
    ETmax:

    Notes
    -----
    This model takes  a lumped approach to the unsaturated zone above the
    water table, where storage increases with infiltration from precipitation:
        f = min(P,Ksat)
    and with declining water table elevation
        r = -sw*dH/dt
    where r is the rate of "inflow" as the unsaturated zone thickness increases
    at a rate -dH/dt.

    Storage in the reservoir decreases with evapotranspiration and leakage:
        qet = (s-sw)/(1-sw)*ETmax
        qr =
    where s is the relative water content at a node [-], and sw is the soil
    wilting point.

    """

    _name = "LumpedUnsaturatedZone"

    _input_var_names = set(("topographic__elevation",
                            "water_table__elevation"))

    _output_var_names = set(
        ("unsaturated__thickness", "infiltration__rate","relative_saturation",
         "soil_water__storage","")
    )

    _var_units = {
        "topographic__elevation": "m",
        "water_table__elevation": "m",
        "unsaturated__thickness": "m",
        "infiltration__rate": "m/s",
        "relative_saturation": "-",
        "soil_water__storage": "m"
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "water_table__elevation": "node",
        "unsaturated__thickness": "node",
        "infiltration__rate": "node",
        "relative_saturation": "node",
        "soil_water__storage": "node"
    }

    _var_doc = {
        "topographic__elevation": "elevation of land surface",
        "water_table__elevation": "elevation of water table",
        "unsaturated__thickness": "thickness of zone above Boussinesq aquifer",
        "infiltration__rate": "rate of water leaving unsaturated zone to aquifer",
        "relative_saturation": "relative saturation of unsat zone [0<=sat<=1]",
        "soil_water__storage": "storage of water in the unsat zone"
    }



    def __init__(self, grid, porosity=0.5, soil_wilting_point=0.18,
            soil_field_capacity=0.3, plant_rooting_depth=1.0, ETmax=2.0e-7,
            hydraulic_conductivity=0.005):

        # Store grid
        self._grid = grid

        # Convert parameters to fields if needed, and store a reference
        self.n = return_array_at_node(grid, porosity)
        self.sw = return_array_at_node(grid, soil_wilting_point)
        self.sf = return_array_at_node(grid, soil_field_capacity)
        self.dr = return_array_at_node(grid, plant_rooting_depth)
        self.ETmax = return_array_at_node(grid, ETmax)
        self.Ksat = return_array_at_node(grid, hydraulic_conductivity)

        # Create fields:
        if "topographic__elevation" in self.grid.at_node:
            self.elev = self.grid.at_node["topographic__elevation"]
        else:
            self.elev = self.grid.add_ones("node", "topographic__elevation")

        if "water_table__elevation" in self.grid.at_node:
            self.wtable = self.grid.at_node["water_table__elevation"]
        else:
            self.wtable = self.grid.add_zeros("node", "water_table__elevation")
            self.wtable[grid.closed_boundary_nodes] = 0

        if "unsaturated__thickness" in self.grid.at_node:
            self.thickness = self.grid.at_node["unsaturated__thickness"]
        else:
            self.thickness = self.elev - self.wtable
            self.thickness[grid.closed_boundary_nodes] = 0

        if "infiltration__rate" in self.grid.at_node:
            self.f = self.grid.at_node["infiltration__rate"]
        else:
            self.f = self.grid.add_zeros("node","infiltration__rate")

        if "relative_saturation" in self.grid.at_node:
            self.sat = self.grid.at_node["relative_saturation"]
        else:
            self.sat = self.sf*self.grid.add_ones("node", "relative_saturation")
            self.sat[grid.closed_boundary_nodes] = 0

        if "soil_water__storage" in self.grid.at_node:
            self.S = self.grid.at_node["soil_water__storage"]
        else:
            self.S = self.grid.add_zeros("node", "soil_water__storage")
            self.S[:] = self.sat*self.thickness*self.n

    def correct_soil_water__storage(self,dt):
        dhdt = self._grid.at_node["water_table__velocity"]
        self.r = self.n*self.sf*dhdt
        self.S += self.r*dt

    def run_one_step(self,duration,intensity):
        Sw = self.sw*self.thickness*self.n #m storage at the wilting point
        Sf = self.sf*self.thickness*self.n #m storage at the field capacity
        St = 1*self.thickness*self.n #m total available storage

        dt = duration
        self.p = intensity
        self.f = np.minimum(self.p,self.Ksat)
        self.ho = np.maximum(self.p-self.Ksat,0)

        self.qet = EvapotranspirationLoss(intensity,self.S,Sw,St,self.ETmax,dt)
        self.ql = LeakageLoss(self.S,self.f,self.qet,Sf,self.Ksat,dt)

        self.S += self.f*dt - self.qet*dt - self.ql*dt

        self.sat[self.thickness > 0.] = self.S[self.thickness > 0.]/self.thickness[self.thickness > 0.]
        self.sat[self.thickness == 0] = 1.0

    def get_recharge_rate(self):
        return self.ql

    def calc_hortonian_overland_flux(self):
        return np.sum(self.ho)

    def calc_et_flux(self):
        return np.sum(self.qet)
