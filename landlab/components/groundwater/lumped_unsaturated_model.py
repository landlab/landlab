# -*- coding: utf-8 -*-
"""
LumpedUnsaturatedZone Component

@author: D Litwin
"""

import numpy as np
from landlab import RasterModelGrid, Component
from landlab.utils import return_array_at_node

def _regularize_R(u):
    return u*np.greater_equal(u,0)

def _EvapotranspirationLoss(p,S,Sw,Sf,St,ETmax):
    qet = np.zeros_like(S)
    i = St>0.0
    qet[i] = (p==0.0)*np.minimum(
                _regularize_R((S[i]-Sw[i])/(Sf[i]-Sw[i]))*ETmax[i],ETmax[i])
    return qet


def _CapillaryRise(S,Sw,Sf,St,ETmax,Zwt,Z,hroot):
    qc = np.zeros_like(S)
    i = St>0.0
    qc[i] = (Z[i]-Zwt[i] < hroot[i])*(ETmax[i] - np.minimum(
                _regularize_R((S[i]-Sw[i])/(Sf[i]-Sw[i]))*ETmax[i],ETmax[i]))

    qc[~i] = ETmax[~i]
    return qc

def _LeakageLoss(S,Sf,St,Ksat,dt):
    #leakage under unit hydraulic gradient and hydraulic
    #conductivity that varies linearly from 0 at field capacity to Ksat
    #at full saturation.
    ql = np.zeros_like(S)
    i = St>0
    ql[i] =  np.minimum(_regularize_R((S[i]-Sf[i])/(St[i]-Sf[i]))*Ksat[i], np.maximum((S[i]-Sf[i])/dt,np.zeros_like(S)[i]))
    return ql

#untested
def _LeakageLoss_MVG(S,Sf,Sw,St,Ko,n,L):
    #leakage loss under unit hydraulic gradient and hydraulic
    #conductivity given by the Mualem-van Genuchten model K(Se)
    Se = np.zeros_like(S)
    i = St>0
    Se[i] = _regularize_R((S[i]-Sw[i])/(St[i]-Sw[i]))
    ql =  1*Ko*(Se**L)*(1 - (1 - Se**(n/(n-1)) )**(1-1/n) )**2
    return ql

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
    Qcapmax:
    Qrootmax:

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
            soil_field_capacity=0.3, plant_rooting_depth=1.0, ETmax=2.0e-8,
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
            self.wtable[grid.closed_boundary_nodes] = 0.0

        if "unsaturated__thickness" in self.grid.at_node:
            self.thickness = self.grid.at_node["unsaturated__thickness"]
        else:
            self.thickness = self.elev - self.wtable
            self.thickness[grid.closed_boundary_nodes] = 0.0

        if "infiltration__rate" in self.grid.at_node:
            self.f = self.grid.at_node["infiltration__rate"]
        else:
            self.f = self.grid.add_zeros("node","infiltration__rate")

        if "relative_saturation" in self.grid.at_node:
            self.sat = self.grid.at_node["relative_saturation"]
        else:
            self.sat = self.sf*self.grid.add_ones("node", "relative_saturation")
            self.sat[grid.closed_boundary_nodes] = 0.0

        if "soil_water__storage" in self.grid.at_node:
            self.S = self.grid.at_node["soil_water__storage"]
        else:
            self.S = self.grid.add_zeros("node", "soil_water__storage")
            self.S[:] = self.sat*self.thickness*self.n

    def correct_soil_water__storage(self,dt):
        dhdt = self._grid.at_node["water_table__velocity"]
        self.r = (dhdt<0)*self.n*self.sf*(-dhdt) - (dhdt>0)*self.n*self.sat*dhdt
        self.S += self.r*dt

    def run_one_step(self,duration,intensity):

        self.thickness = self.elev - self.wtable

        Sw = self.sw*self.thickness*self.n #m storage at the wilting point
        Sf = self.sf*self.thickness*self.n #m storage at the field capacity
        St = 1*self.thickness*self.n #m total available storage

        dt = duration
        self.p = intensity
        self.f = np.minimum(self.p,self.Ksat)
        self.ho = np.maximum(self.p-self.Ksat,0.)

        self.qet = _EvapotranspirationLoss(intensity,self.S,Sw,Sf,St,self.ETmax)
        self.ql = _LeakageLoss(self.S,Sf,St,self.Ksat,dt)
        self.qc = _CapillaryRise(self.S,Sw,Sf,St,self.ETmax,self.wtable,self.elev,self.dr)
        self.S += self.f*dt + self.qc*dt - self.qet*dt - self.ql*dt

        j = self.thickness > 0.0
        self.sat[j] = self.S[j]/(self.thickness[j]*self.n[j])
        self.sat[~j] = 0.0

    def get_recharge_rate(self):
        return self.ql

    def get_capillary_loss(self):
        return self.qc

    def calc_hortonian_overland_flux(self):
        return np.sum(self.ho[self._grid.core_nodes]*self._grid.area_of_cell)

    def calc_et_flux(self):
        return np.sum(self.qet[self._grid.core_nodes]*self._grid.area_of_cell)

    def calc_cap_flux(self):
        return np.sum(self.qc[self._grid.core_nodes]*self._grid.area_of_cell)

    def calc_recharge_flux(self):
        return np.sum(self.ql[self._grid.core_nodes]*self._grid.area_of_cell)

    def calc_total_storage(self):
        return np.sum(self.S[self._grid.core_nodes]*self._grid.area_of_cell)
