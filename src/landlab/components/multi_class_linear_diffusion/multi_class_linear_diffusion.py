#! /usr/env/python
"""
"""

import numpy as np
from landlab.components import LinearDiffuser
from landlab import Component
from landlab import LinkStatus
from landlab import NodeStatus
from landlab import RasterModelGrid

_ALPHA = 0.15


class MultiClassLinearDiffusion(LinearDiffuser):
    """Calculate linear difffusion and spread sediment mass across grains classes.
    This component call the LinearDiffuser component to update sediment flux by diffusion and
    then update the mass accordingly.

    The primary method of this class is : `run_one_step`.

        Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import MultiClassLinearDiffusion
    >>> porosity = 0.4
    >>> rho_sed = 2650.0
    >>> grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> bedrock = grid.add_zeros("bedrock__elevation", at="node")
    >>> grains = np.zeros((grid.number_of_nodes, 2))
    >>> grains[grid.core_nodes] = [500.0, 1500.0]  # in kg/m2
    >>> grains[7] = [3000.0, 1000.0]
    >>> _ = grid.add_field("grains__mass", grains, at="node")
    >>> soil = grid.add_field(
    ...     "soil__depth",
    ...     grains.sum(axis=1) / ((1 - porosity) * rho_sed),
    ...     at="node",
    ... )
    >>> z = grid.add_field("topographic__elevation", bedrock + soil, at="node")
    >>> soil[grid.core_nodes].round(3)
    >>> total_mass = grid.at_node["grains__mass"][grid.core_nodes].sum()
    >>> diffuser = MultiClassLinearDiffusion(
    ...     grid, linear_diffusivity=0.1, phi=porosity, rho_sed=rho_sed
    ... )
    >>> diffuser.run_one_step(dt=500.0)
    >>> soil[grid.core_nodes].round(3)
    array([1.613, 1.805, 1.613])
    """

    _name = "MultiClassLinearDiffusion"

    _unit_agnostic = True

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Topographic elevation at node",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Bedrock elevation at node",
        },
        "hillslope_sediment__unit_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2/s",
            "mapping": "link",
            "doc": "Volume flux per unit width along links",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "grains__mass": {
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
            Bulk porosity of the sediment layer
        rho_sed : float (default 2650)
            Sediment density [kg/m^3]
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

        super().__init__(grid, **kwds)

        if np.ndim(self._grid.at_node['grains__mass']) > 1:
            self._n_classes = np.shape(self._grid.at_node['grains__mass'])[1]
        else:
            self._n_classes = 1
        face_width_at_link = np.zeros(self._grid.number_of_links)
        has_face = self._grid.face_at_link != -1
        face_width_at_link[has_face] = self._grid.length_of_face[self._grid.face_at_link[has_face]]
        links = self._grid.links_at_node
        self._face_width_at_node_links = np.where(links != -1, face_width_at_link[links], 0.0)
        area = self._grid.cell_area_at_node
        self._inv_cell_area_at_node = np.divide(1.0, area, out=np.zeros_like(area), where=area > 0)
        self._dx_min = self._grid.length_of_link[self._grid.active_links].min()
        self._zeros_at_node = self._grid.zeros(at="node")
        self._zeros_at_link = self._grid.zeros(at="link")
        self._zeros_at_link_for_fractions = np.zeros(
            (np.shape(self._zeros_at_link)[0], self._n_classes))

    def _spread_mass_between_classes(self, dt,
                                     topo):

        # Pointers
        sediment_flux = self._grid.at_link["hillslope_sediment__unit_volume_flux"] * dt
        bed = self._grid.at_node['bedrock__elevation']
        soil = self._grid.at_node['soil__depth']
        grain_mass_node = self._grid.at_node['grains__mass']
        grain_mass_node[grain_mass_node < 0] = 0
        n_classes = self._n_classes

        # Init arrays
        g_total_link = self._zeros_at_link.copy()  # Total mass at link
        g_state_link = self._zeros_at_link_for_fractions.copy()  # Grain mass of each class
        sed_flux_at_link = self._zeros_at_link.copy()
        self._sum_dzdt = self._zeros_at_node.copy()

        # Calc topographic gradient at link
        tgrad = self.grid.calc_grad_at_link(topo)

        # Get the upwind node ID at link
        upwind_node_id_at_link = self._grid.map_value_at_max_node_to_link(topo,
                                                                          self._grid.nodes.flatten())
        # Get only valid upwind nodes
        nonzero_upwind_node_ids = upwind_node_id_at_link.astype(int)

        # Get the LINKS IDs 
        nonzero_link_ids = np.arange(self._grid.number_of_links)

        # Get the outlinks in each node
        outlinks_at_node = self.grid.link_at_node_is_downwind(tgrad)
        if n_classes > 1:
            g_total_dt_node = np.sum(grain_mass_node, 1).copy()  # Total grains mass at node
            g_total_link[nonzero_link_ids] = g_total_dt_node[
                nonzero_upwind_node_ids]  # Total sediment mass for of each up-wind node (mapped to link).
            g_state_link[nonzero_link_ids, :] = (
                grain_mass_node)[nonzero_upwind_node_ids, :]  # Sediment mass for all size-fraction, mapped to link
            g_fraction_link = np.divide(g_state_link,
                                        g_total_link.reshape(-1, 1),
                                        out=np.zeros_like(g_state_link),
                                        where=g_total_link.reshape(-1, 1) != 0)
            self._sed_flux_at_link_class = np.multiply(np.abs(sediment_flux.reshape([-1, 1])),
                                                       g_fraction_link) * -np.sign(tgrad[:, np.newaxis])
        else:
            g_total_link[nonzero_link_ids] = grain_mass_node[
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
        # sum_fluxes_out = np.abs(np.sum(np.abs(outlinks_fluxes_at_node), 1))
        sum_fluxes_out = np.abs(
            np.sum(np.abs(outlinks_fluxes_at_node) * self._face_width_at_node_links[:, :, np.newaxis], 1))

        dz_per_grainsize_at_node = grain_mass_node / ((1 - self._phi) * self._rho_sed)
        if n_classes > 1:
            indices_to_correct_flux = np.where(
                (sum_fluxes_out * self._inv_cell_area_at_node[:, np.newaxis]) > dz_per_grainsize_at_node)
        else:
            indices_to_correct_flux = np.where(
                (sum_fluxes_out[:, 0] * self._inv_cell_area_at_node) > dz_per_grainsize_at_node)

        if np.any(indices_to_correct_flux):
            # If the outflux is greater than the mass in the upstream node we should fix the outflux.
            if n_classes > 1:
                ratios = np.divide(dz_per_grainsize_at_node[indices_to_correct_flux],
                                   sum_fluxes_out[indices_to_correct_flux] * self._inv_cell_area_at_node[
                                       indices_to_correct_flux[0]])
                outlinks_fluxes_at_node[indices_to_correct_flux[0], :, indices_to_correct_flux[1]] *= ratios[:,
                np.newaxis]
            else:
                ratios = np.divide(dz_per_grainsize_at_node[indices_to_correct_flux][:, np.newaxis],
                                   sum_fluxes_out[indices_to_correct_flux] *
                                   self._inv_cell_area_at_node[indices_to_correct_flux[0]][:, np.newaxis])

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
                grain_mass_node[:, size_class] += (dz) * (1 - self._phi) * self._rho_sed  # in kg/m2
        else:
            dz = -self._grid.calc_flux_div_at_node(self._sed_flux_at_link_class[:, 0])
            grain_mass_node[:] += (dz) * (1 - self._phi) * self._rho_sed  # in kg/m2

        grain_mass_node[grain_mass_node < 0] = 0

        # Return the depth change in soil layer
        return dz_soil_sediment

    def run_one_step(self, dt):
        """
        Advance by one time step

        Parameters
        ----------
        dt : float
            Time-step duration (y)
        """
        # dt_max = _ALPHA * self._grid.length_of_link[0] ** 2 / self._kd
        dt_max = _ALPHA * self._dx_min ** 2 / self._kd
        n = int(dt // dt_max) + 1
        for _ in range(n):
            self._one_substep(dt / n)

    def _one_substep(self, dt):
        """
        Advance by one substep

        Parameters
        ----------
        dt : float
            Time-step duration (y)
        """
        topo_before = self.grid.at_node["topographic__elevation"].copy()
        soil_before = self.grid.at_node["soil__depth"].copy()

        super().run_one_step(dt)
        dz_soil_sediment = self._spread_mass_between_classes(dt=dt,
                                                             topo=topo_before)

        # Update topography after soil diffusion
        self._grid.at_node['soil__depth'][:] = soil_before + dz_soil_sediment
        self._grid.at_node['soil__depth'][self._grid.at_node['soil__depth'] <= 0] = 0
        self._grid.at_node['topographic__elevation'][:] = (self._grid.at_node['soil__depth'][:] +
                                                           self._grid.at_node['bedrock__elevation'][:])


