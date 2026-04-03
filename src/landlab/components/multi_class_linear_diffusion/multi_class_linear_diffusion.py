#! /usr/env/python
"""
"""


import numpy as np
from landlab.components import LinearDiffuser
from landlab import Component
from landlab import LinkStatus
from landlab import NodeStatus
from landlab import RasterModelGrid


class MultiClassLinearDiffusion(LinearDiffuser):

    _name = "MultiClassLinearDiffusion"

    _unit_agnostic = True

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Soil depth",
        },
        "grains__weight": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "kg/m2",
            "mapping": "node",
            "doc": "Mass per unit area of grains in each size class stored at node",
        },
    }


    def __init__(self,
                 grid,
                 linear_diffusivity=0.01,
                 phi=0.35,
                 rho_sed=2650,
                 **kwds):
        """
        Parameters
        ----------
        grid : ModelGrid
            A grid.
        linear_diffusivity : float
            Sediment diffusivity. If an array or field name, these must be the
            diffusivities on either nodes or links - the component will
            distinguish which based on array length. Values on nodes will be
            mapped to links using an upwind scheme in the simple case.
        phi : float  (default 0.35)
            Bulk porosity of bed sediment
        rho_sed : float (default 2650)
            Density of the sediment layer
        """

        # if "kd" not in grid.at_node:
        #     grid.add_zeros("kd", at="node")
        #
        # kwds.setdefault("linear_diffusivity", "kd")
        # grid.at_node["kd"][:] = linear_diffusivity

        kwds.setdefault("linear_diffusivity", linear_diffusivity)
        self._kd = linear_diffusivity
        self._rho_sed = rho_sed
        self._phi = phi

        super().__init__(grid,**kwds)

        if np.ndim(self._grid.at_node['grains__weight']) > 1:
            self._n_classes = np.shape(self._grid.at_node['grains__weight'])[1]
        else:
            self._n_classes = 1
        self._zeros_at_node = self._grid.zeros(at="node")
        self._zeros_at_link = self._grid.zeros(at="link")
        self._zeros_at_link_for_fractions = np.zeros(
            (np.shape(self._zeros_at_link)[0], self._n_classes))


    def _spread_mass_between_classes(self, dt,
                                     ):

        # Pointers
        sediment_flux = self._grid.at_link["hillslope_sediment__unit_volume_flux"] * dt
        bed = self._grid.at_node['bedrock__elevation']
        soil = self._grid.at_node['soil__depth']
        topo = self._grid.at_node['topographic__elevation']
        grain_weight_node = self._grid.at_node['grains__weight']
        grain_weight_node[grain_weight_node < 0] = 0
        n_classes = self._n_classes

        # Init arrays
        g_total_link = self._zeros_at_link.copy()  # Total mass at link
        g_state_link = self._zeros_at_link_for_fractions.copy()  # Grain mass of each class
        sed_flux_at_link = self._zeros_at_link.copy()
        self._sum_dzdt = self._zeros_at_node.copy()

        # Calc topographic gradient at link
        tgrad = self.grid.calc_grad_at_link(self.grid.at_node['topographic__elevation'])

        # Get the upwind node ID at link
        upwind_node_id_at_link = self._grid.map_value_at_max_node_to_link('topographic__elevation',
                                                                          self._grid.nodes.flatten() )
        # Get only valid upwind nodes
        nonzero_upwind_node_ids = np.uint(upwind_node_id_at_link[np.nonzero(upwind_node_id_at_link)])
        
        # Get the LINKS IDs 
        nonzero_link_ids = np.nonzero(upwind_node_id_at_link)[0]

        # Get the outlinks in each node
        outlinks_at_node = self.grid.link_at_node_is_downwind(tgrad)
        if n_classes > 1:
            g_total_dt_node = np.sum(grain_weight_node, 1).copy()  # Total grains mass at node
            g_total_link[nonzero_link_ids] = g_total_dt_node[
                nonzero_upwind_node_ids]  # Total sediment mass for of each up-wind node mapped to link.
            g_state_link[nonzero_link_ids, :] = grain_weight_node[nonzero_upwind_node_ids,
            :]  # Sediment mass for all size-fraction, mapped to link
            g_fraction_link = np.divide(g_state_link,
                                        g_total_link.reshape(-1, 1),
                                        out=np.zeros_like(g_state_link),
                                        where=g_total_link.reshape(-1, 1) != 0)
            self._sed_flux_at_link_class = np.multiply(np.abs(sediment_flux.reshape([-1, 1])),
                                                       g_fraction_link) * -np.sign(tgrad[:, np.newaxis])
        else:
            g_total_link[nonzero_link_ids] = grain_weight_node[
                nonzero_upwind_node_ids]  # Total sediment mass for of each up-wind node mapped to link.
            self._sed_flux_at_link_class = np.abs(sediment_flux.reshape([-1, 1])) * -np.sign(tgrad[:, np.newaxis])

        fluxes = (self._grid.link_dirs_at_node[:, :, np.newaxis] *
                  self._sed_flux_at_link_class[self._grid.links_at_node, :])
        outlinks_fluxes_at_node = np.copy(fluxes)
        outlinks_fluxes_at_node[~outlinks_at_node, :] = 0
        outlinks_id = self.grid.links_at_node[outlinks_at_node]
        inlinks_fluxes_at_node = np.copy(fluxes)
        inlinks_fluxes_at_node[outlinks_at_node, :] = 0

        # Sum all outfluxes per node to check if transport rate is greater than the mass
        # stored in the upwind node
        sum_fluxes_out = np.abs(np.sum(np.abs(outlinks_fluxes_at_node), 1))
        dz_per_grainsize_at_node = np.abs(grain_weight_node / (1 - self._phi)) * self.grid.area_of_cell[0]
        if n_classes > 1:
            indices_to_correct_flux = np.where((sum_fluxes_out / self._grid.length_of_link[0]) > dz_per_grainsize_at_node)
        else:
            indices_to_correct_flux = np.where(
                (sum_fluxes_out[:, 0] / self._grid.length_of_link[0]) > dz_per_grainsize_at_node)

        if np.any(indices_to_correct_flux):

            # If the outflux is greater than the mass in the upstream node we should fix the outflux.
            if n_classes > 1:
                ratios = np.divide(dz_per_grainsize_at_node[indices_to_correct_flux],
                                   sum_fluxes_out[indices_to_correct_flux])
                outlinks_fluxes_at_node[indices_to_correct_flux[0], :, indices_to_correct_flux[1]] *= ratios[:,
                np.newaxis]
            else:
                ratios = np.divide(dz_per_grainsize_at_node[indices_to_correct_flux][:, np.newaxis],
                                   sum_fluxes_out[indices_to_correct_flux])

                outlinks_fluxes_at_node[indices_to_correct_flux[0], :] *= ratios[:, np.newaxis]
            # Update the flux.
            self._sed_flux_at_link_class[outlinks_id] = np.abs(outlinks_fluxes_at_node[outlinks_at_node, :]) * np.sign(
                self._sed_flux_at_link_class[outlinks_id])

        sed_flux_at_link[:] = np.sum(self._sed_flux_at_link_class, axis=1)
        dz_soil_sediment = -self._grid.calc_flux_div_at_node(sed_flux_at_link)

        # Update mass across grains classes
        if n_classes > 1:
            for size_class in range(n_classes):
                dz = -self._grid.calc_flux_div_at_node(self._sed_flux_at_link_class[:, size_class])
                grain_weight_node[:, size_class] += (dz) * (1 - self._phi) * self._rho_sed  # in kg/m2
        else:
            dz = -self._grid.calc_flux_div_at_node(self._sed_flux_at_link_class[:, 0])
            grain_weight_node[:] += (dz) * (1 - self._phi) * self._rho_sed  # in kg/m2

        grain_weight_node[grain_weight_node < 0] = 0
        bedrock_sediment_flux_at_link = np.abs((np.abs(sediment_flux) - np.abs(sed_flux_at_link))) * np.sign(
            sed_flux_at_link)

        # Return the depth change in soil layer
        return dz_soil_sediment


    def calc_rock_exposure_fraction(self):
        """Update the bedrock exposure fraction.
        """
        self._rock_exposure_fraction = np.exp(-self._grid.at_node['soil__depth'] /
                                              self._depth_decay_scale)



    def run_one_step(self, dt):
        """
        Advance by one time step.

        Parameters
        ----------
        dt : float
            Time-step duration (y)
        """
        topo_before = self.grid.at_node["topographic__elevation"].copy()
        soil_before = self.grid.at_node["soil__depth"].copy()

        super().run_one_step(dt)
        dz_soil_sediment = self._spread_mass_between_classes(dt=dt)

        # Update topography after soil diffusion
        self._grid.at_node['soil__depth'][:] = soil_before + dz_soil_sediment
        self._grid.at_node['soil__depth'][self._grid.at_node['soil__depth'] <= 0] = 0
        self._grid.at_node['topographic__elevation'][:] = (self._grid.at_node['soil__depth'][:] +
                                                           self._grid.at_node['bedrock__elevation'][:])