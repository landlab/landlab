# -*- coding: utf-8 -*-
"""
GroundwaterDupuitPercolator Component

@author: G Tucker, D Litwin
"""

import numpy as np
from landlab import Component, BAD_INDEX_VALUE
from landlab.components import FlowAccumulator
from landlab.utils import return_array_at_node, return_array_at_link
from landlab.grid.mappers import map_mean_of_link_nodes_to_link, \
        map_max_of_node_links_to_node, map_value_at_max_node_to_link, map_link_head_node_to_link
from landlab.grid.raster_mappers import map_sum_of_outlinks_to_node


ACTIVE_LINK = 0


#regularization functions used to deal with numerical demons of seepage
def _regularize_G(u,reg_factor):
    return np.exp(-(1-u)/reg_factor)


def _regularize_R(u):
    return u*np.greater_equal(u,0)

def get_link_hydraulic_conductivity(grid,K):
    #where grid is a landlab grid object, and K is a 2D hydraulic conductivity tensor
    u = grid.unit_vector_at_link
    K_link = np.zeros(len(u))
    for i in range(len(u)):
        K_link[i] = np.dot(np.dot(u[i,:],K),u[i,:])
    return K_link

class GroundwaterDupuitPercolator(Component):
    """
    Simulate groundwater flow in a shallow unconfined aquifer.

    The GroundwaterDupuitPercolator uses the Dupuit approximation that the
    hydraulic gradient is equal to the slope of the water table.

    Parameters
    ----------
    grid: ModelGrid
            Landlab ModelGrid object
    hydraulic_conductivity: float, field name, or array of float
            saturated hydraulic conductivity, m/s
            Default = 0.01 m/s
    recharge_rate: float, field name, or array of float
            Rate of recharge, m/s
            Default = 1.0e-8 m/s
    regularization_f: float
            factor controlling the smoothness of the transition between
            surface and subsurface flow

       Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> gdp = GroundwaterDupuitPercolator(mg)

    Notes
    -----
    Groundwater discharge per unit length, q, is calculated as:

        q = - K H dw/dx,

    where K is hydraulic conductivity, H is aquifer thickness, w is water table
    height, and x is horizontal distance.

    An explicit forward-in-time finite-volume method is used to implement a
    numerical solution. Flow discharge between neighboring nodes is calculated
    using the average depth at the nodes.

    Note that the current version does NOT handle surface seepage.
    """

    _name = "GroundwaterDupuitPercolator"

    _input_var_names = set(("topographic__elevation",
                            "aquifer_base__elevation"))

    _output_var_names = set(
        ("aquifer__thickness", "water_table__elevation","aquifer_base__gradient",
         "hydraulic__gradient", "groundwater__specific_discharge",
         "groundwater__velocity", "surface_water__specific_discharge","")
    )

    _var_units = {
        "topographic__elevation": "m",
        "aquifer_base__elevation": "m",
        "aquifer__thickness": "m",
        "aquifer_base__gradient": "m/m",
        "water_table__elevation": "m",
        "hydraulic__gradient": "m/m",
        "groundwater__specific_discharge": "m2/s",
        "groundwater__velocity": "m/s",
        "surface_water__specific_discharge": "m/s",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "aquifer_base__elevation": "node",
        "aquifer__thickness": "node",
        "water_table__elevation": "node",
        "hydraulic__gradient": "link",
        "aquifer_base__gradient":"link",
        "groundwater__specific_discharge": "link",
        "groundwater__velocity": "link",
        "surface_water__specific_discharge": "node",
    }

    _var_doc = {
        "topographic__elevation": "elevation of land surface",
        "aquifer_base__elevation": "elevation of impervious layer",
        "aquifer__thickness": "thickness of saturated zone",
        "water_table__elevation": "elevation of water table",
        "hydraulic__gradient": "gradient of water table in link direction",
        "aquifer_base__gradient":"gradient of the aquifer base in the link direction",
        "groundwater__specific_discharge": "discharge per width in link dir",
        "groundwater__velocity": "velocity of groundwater in link direction",
        "surface_water__specific_discharge": "rate of seepage to surface",
    }

    def __init__(self, grid, hydraulic_conductivity=0.01,porosity=0.5,
                 recharge_rate=1.0e-8,regularization_f=1E-3,unsaturated_zone=False):
        """Initialize the GroundwaterDupuitPercolator.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        hydraulic_conductivity: float, field name, or array of float
                saturated hydraulic conductivity, m/s
                Default = 0.01 m/s
        porosity: float, field name or array of float
                the porosity of the aquifer [-]
                Default = 0.5
        recharge_rate: float, field name, or array of float
                Rate of recharge, m/s
                Default = 1.0e-8 m/s
        regularization_f: float
                factor controlling the smoothness of the transition between
                surface and subsurface flow
        """
        # Store grid
        self._grid = grid

        # Shorthand
        self.cores = grid.core_nodes
        self.inactive_links = np.where(grid.status_at_link != ACTIVE_LINK)[0]

        # Convert parameters to fields if needed, and store a reference
        self.K = return_array_at_link(grid, hydraulic_conductivity)
        self.recharge = return_array_at_node(grid, recharge_rate)
        self.n = return_array_at_node(grid,porosity)
        self.n_link = map_mean_of_link_nodes_to_link(self._grid,self.n)
        self.r = regularization_f
        self.unsat = unsaturated_zone
        # Create fields:

        if "topographic__elevation" in self.grid.at_node:
            self.elev = self.grid.at_node["topographic__elevation"]
        else:
            self.elev = self.grid.add_ones("node", "topographic__elevation")

        if "aquifer_base__elevation" in self.grid.at_node:
            self.base = self.grid.at_node["aquifer_base__elevation"]
        else:
            self.base = self.grid.add_zeros("node", "aquifer_base__elevation")

        if "water_table__elevation" in self.grid.at_node:
            self.wtable = self.grid.at_node["water_table__elevation"]
        else:
            self.wtable = self.grid.add_zeros("node", "water_table__elevation")

        self.wtable[grid.closed_boundary_nodes] = 0

        if "aquifer__thickness" in self.grid.at_node:
            self.thickness = self.grid.at_node["aquifer__thickness"]
        else:
            self.thickness = self.grid.add_zeros("node", "aquifer__thickness")
            self.thickness[:] = self.wtable - self.base

        if "hydraulic__gradient" in self.grid.at_link:
            self.hydr_grad = self.grid.at_link["hydraulic__gradient"]
        else:
            self.hydr_grad = self.grid.add_zeros("link", "hydraulic__gradient")

        if "aquifer_base__gradient" in self.grid.at_link:
            self.base_grad = self.grid.at_link["aquifer_base__gradient"]
        else:
            self.base_grad = self.grid.add_zeros("link", "aquifer_base__gradient")

        if "groundwater__specific_discharge" in self.grid.at_link:
            self.q = self.grid.at_link["groundwater__specific_discharge"]
        else:
            self.q = self.grid.add_zeros("link",
                                         "groundwater__specific_discharge")

        if "groundwater__velocity" in self.grid.at_link:
            self.vel = self.grid.at_link["groundwater__velocity"]
        else:
            self.vel = self.grid.add_zeros("link", "groundwater__velocity")

        if "surface_water__specific_discharge" in self.grid.at_node:
            self.qs = self.grid.at_node["surface_water__specific_discharge"]
        else:
            self.qs = self.grid.add_zeros("node", "surface_water__specific_discharge")

        if "water_table__velocity" in self.grid.at_node:
            self.dhdt = self.grid.at_node["water_table__velocity"]
        else:
            self.dhdt = self.grid.add_zeros("node", "water_table__velocity")

        #just shoving these in here temporarily
        self.S = abs(grid.calc_grad_at_link(self.elev))
        self.S_node = map_max_of_node_links_to_node(grid,self.S)


    def calc_rechage_flux_in(self):
        #  calculate flux into the domain from recharge. Includes recharge that
        #  may immediately become saturation excess overland flow.
        return np.sum(self._grid.area_of_cell*self.recharge[self.cores])

    def calc_gw_flux_out(self):

        """
        Groundwater flux through open boundaries may be positive (out of
        the domain) or negative (into the domain). This function determines the
        correct sign for specific discharge based upon this convention,
        and sums the flux across the boundary faces.

        """

        open_nodes = self._grid.open_boundary_nodes
        links_at_open = self._grid.links_at_node[open_nodes]
        link_dirs_at_open = self._grid.active_link_dirs_at_node[open_nodes]

        active_links_at_open = links_at_open[link_dirs_at_open!=0]
        active_link_dirs_at_open = link_dirs_at_open[link_dirs_at_open!=0]

        q_at_open_links = self._grid.at_link['groundwater__specific_discharge'][active_links_at_open]

        faces = self._grid.face_at_link[active_links_at_open]
        face_widths = self._grid.width_of_face[faces]

        gw_volume_flux_rate_out = q_at_open_links*active_link_dirs_at_open*face_widths

        return np.sum(gw_volume_flux_rate_out)

    def calc_sw_flux_out(self):
        # surface water flux out of the domain through seepage and saturation excess.
        # Note that model does not allow for reinfiltration.

        return np.sum(self._grid.at_node['surface_water__discharge'][self._grid.open_boundary_nodes] )

    def calc_gw_flux_at_node(self):
        # the sum of the flux of groundwater leaving a node
        gw = self._grid.at_link['groundwater__specific_discharge'][self._grid.links_at_node]*self._grid.link_dirs_at_node
        gw_out = -np.sum(gw*(gw<0),axis=1)
        return gw_out

        # Old definition of gw flux at node.
        # return map_max_of_node_links_to_node(self._grid,self._grid.dx* abs(self._grid.at_link['groundwater__specific_discharge']))

    def calc_sw_flux_at_node(self):
        return self._grid.at_node['surface_water__discharge']

    def calc_shear_stress_at_node(self,n_manning=0.05):
        rho = 1000
        g = 9.81
        return rho*g*self.S_node *( (n_manning*self._grid.at_node['surface_water__discharge']/3600)/(self._grid.dx*np.sqrt(self.S_node)) )**(3/5)


    def calc_total_storage(self):
        # calculate the current water storage in the aquifer
        return np.sum(self.n[self.cores] * self._grid.area_of_cell *
                        self._grid.at_node['aquifer__thickness'][self.cores])


    def run_one_step(self, dt, s=None,sf=None,**kwds):
        """
        Advance component by one time step of size dt.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        """

        #Calculate base gradient
        self.base_grad[:] = self._grid.calc_grad_at_link(self.base)
        self.base_grad[self.inactive_links] = 0.0

        # Calculate hydraulic gradient
        self.hydr_grad[:] = self._grid.calc_grad_at_link(self.thickness)
        self.hydr_grad[self.inactive_links] = 0.0

        # Calculate groundwater velocity
        self.vel[:] = -self.K *(self.hydr_grad*np.cos(np.arctan(abs(self.base_grad)))
                                    + np.sin(np.arctan(self.base_grad)))
        self.vel[self._grid.status_at_link != 0] = 0.0

        # Aquifer thickness at links (upwind)
        hlink = map_value_at_max_node_to_link(self._grid,
                               'water_table__elevation','aquifer__thickness')

        # Calculate specific discharge
        self.q[:] = hlink * self.vel

        # Groundwater flux divergence
        dqdx = self._grid.calc_flux_div_at_node(self.q)

        # Calculate surface discharge at nodes
        self.qs[:] = _regularize_G(self.thickness/(self.elev-self.base),
                                    self.r)*_regularize_R(self.recharge - dqdx)

        # Mass balance
        if self.unsat==True:
            # unsaturated method interfaces with the LumpedUnsaturatedZone model
            # accounting for the water gained and lost as the water table moves
            # into and out of the unsaturated zone.
            dsdt = self.recharge - dqdx - self.qs
            self.dhdt = (dsdt>=0)*(1/(self.n*(1-s)))*dsdt \
                        + (dsdt<0)*(1/(self.n*(1+sf)))*dsdt
        else:
            self.dhdt = (1/self.n)*(self.recharge - dqdx - self.qs)

        # Update
        self.thickness[self._grid.core_nodes] += self.dhdt[self.cores] * dt

        # Recalculate water surface height
        self.wtable[self._grid.core_nodes] = (self.base[self.cores]
                                              + self.thickness[self.cores])

    def run_with_adaptive_time_step_solver(self, dt, courant_coefficient=0.5,s=None,sf=None,**kwds):
        """
        Advance component by one time step of size dt, subdividing the timestep
        into substeps as necessary to meet a Courant condition.
        Note this method only returns the fluxes at the last subtimestep.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        courant_coefficient: float
            The muliplying factor on the condition that the timestep is
            smaller than the minimum link length over groundwater flow velocity
        """

        remaining_time = dt
        self.num_substeps = 0

        while remaining_time > 0.0:
            #Calculate base gradient
            self.base_grad[:] = self._grid.calc_grad_at_link(self.base)
            self.base_grad[self.inactive_links] = 0.0

            # Calculate hydraulic gradient
            self.hydr_grad[:] = self._grid.calc_grad_at_link(self.thickness)
            self.hydr_grad[self.inactive_links] = 0.0

            # Calculate groundwater velocity
            self.vel[:] = -self.K *(self.hydr_grad*np.cos(np.arctan(self.base_grad))
                                        + np.sin(np.arctan(self.base_grad)))
            self.vel[self._grid.status_at_link != 0] = 0.0

            # Aquifer thickness at links (upwind)
            hlink = map_value_at_max_node_to_link(self._grid,
                                   'water_table__elevation','aquifer__thickness')

            # Calculate specific discharge
            self.q[:] = hlink * self.vel

            # Groundwater flux divergence
            dqdx = self._grid.calc_flux_div_at_node(self.q)

            # Calculate surface discharge at nodes
            self.qs[:] = _regularize_G(self.thickness/(self.elev-self.base),
                                        self.r)*_regularize_R(self.recharge - dqdx)

            # Mass balance
            if self.unsat==True:
                # unsaturated method interfaces with the LumpedUnsaturatedZone model
                # accounting for the water gained and lost as the water table moves
                # into and out of the unsaturated zone.
                dsdt = self.recharge - dqdx - self.qs
                self.dhdt = (dsdt>=0)*(1/(self.n*(1-s)))*dsdt \
                            + (dsdt<0)*(1/(self.n*(1+sf)))*dsdt
            else:
                self.dhdt = (1/self.n)*(self.recharge - dqdx - self.qs)

            #calculate criteria for timestep
            max_vel = max(abs(self.vel/self.n_link))
            grid_dist = min(self._grid.length_of_link)
            substep_dt = np.nanmin([courant_coefficient*grid_dist/max_vel,remaining_time])

            # Update
            self.thickness[self._grid.core_nodes] += self.dhdt[self.cores] * substep_dt

            # Recalculate water surface height
            self.wtable[self._grid.core_nodes] = (self.base[self.cores]
                                                  + self.thickness[self.cores])

            remaining_time -= substep_dt
            self.num_substeps += 1
