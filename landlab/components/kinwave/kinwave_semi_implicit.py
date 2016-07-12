# -*- coding: utf-8 -*-
"""
kinwave_semi_implicit.py: Kinematic wave overland flow component using
semi-implicit numerical solution scheme.

Created on Tue Jul 12 13:42:45 2016

@author: gtucker
"""

from landlab import Component, CORE_NODE, INACTIVE_LINK
import numpy as np

class KinwaveSemiImplicit(Component):
    """
    Calculate water flow over topography.

    Landlab component that implements a two-dimensional 
    kinematic wave model using a semi-implicit numerical scheme.
    
    Construction:
    
        KinwaveSemiImplicit(grid, infilt_rate=0.0,
                            roughness=0.01, **kwds)
    
    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid object.
    infilt_rate : float, optional (defaults to 0)
        Maximum rate of infiltration, mm/hr
    roughnes : float, defaults to 0.01
        Manning roughness coefficient, s/m^1/3
    """

    _name = 'KinwaveOverlandFlowModel'

    _input_var_names = (
        'topographic__elevation',
        'topographic__gradient',
    )

    _output_var_names = (
        'water__depth',
        'water__velocity',
        'water__specific_discharge',
    )

    _var_units = {
        'topographic__elevation' : 'm',
        'topographic__slope' : 'm/m',
        'water__depth' : 'm',
        'water__velocity' : 'm/s',
        'water__specific_discharge' : 'm2/s',
    }

    _var_mapping = {
        'topographic__elevation' : 'node',
        'topographic__gradient' : 'link',
        'water__depth' : 'node',
        'water__velocity' : 'link',
        'water__specific_discharge' : 'link',
    }

    _var_doc = {
        'topographic__elevation':
            'elevation of the ground surface relative to some datum',
        'topographic__gradient':
            'gradient of the ground surface',
        'water__depth':
            'depth of water',
        'water__velocity':
            'flow velocity component in the direction of the link',
        'water__specific_discharge':
            'flow discharge component in the direction of the link',    
    }

    def __init__(self, grid, infilt_rate=0.0, roughness=0.01, **kwds):
        """Initialize the KinwaveOverlandFlowModel.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        infilt_rate : float, optional (defaults to 0)
            Maximum rate of infiltration, mm/hr
        roughnes : float, defaults to 0.01
            Manning roughness coefficient, s/m^1/3
        """

        # Store grid and parameters and do unit conversion
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        self.infilt = infilt_rate / 3600000.0 # convert to m/s
        self.vel_coef = 1.0 / roughness  # do division now to save time

        # Create fields...
        #   Elevation
        if 'topographic__elevation' in grid.at_node:
            self.elev = grid.at_node['topographic__elevation']
        else:
            self.elev = grid.add_zeros('node',
                                       'topographic__elevation')
        #  Water depth
        if 'water__depth' in grid.at_node:
            self.depth = grid.at_node['water__depth']
        else:
            self.depth = grid.add_zeros('node', 'water__depth')
        #   Slope
        if 'topographic__gradient' in grid.at_link:
            self.slope = grid.at_link['topographic__gradient']
        else:
            self.slope = grid.add_zeros('link', 'topographic__gradient')
        #  Velocity
        if 'water__velocity' in grid.at_link:
            self.vel = grid.at_link['water__velocity']
        else:
            self.vel = grid.add_zeros('link', 'water__velocity')
        #  Discharge
        if 'water__specific_discharge' in grid.at_link:
            self.disch = grid.at_link['water__specific_discharge']
        else:
            self.disch = grid.add_zeros('link',
                                        'water__specific_discharge')
        #  Inflow discharge at nodes
        if 'water__specific_discharge' in grid.at_node:
            self.q_in_at_node = grid.at_node['water__inflow_specific_discharge']
        else:
            self.q_in_at_node = grid.add_zeros('node',
                                               'water__inflow_specific_discharge')

        self.calc_ground_slope()

        self.sort_nodes_by_elevation()

    def calc_ground_slope(self):
        """Calculate and store square root of land-surface slope, and sign.
                
        Examples
        --------
        >>> from landlab.components.kinwave.kinwave_semi_implicit import KinwaveSemiImplicit
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((3, 4))
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = np.arange(rg.number_of_nodes)
        >>> kw = KinwaveSemiImplicit(rg)
        >>> kw.slope[rg.active_links]
        array([ 4.,  4.,  1.,  1.,  1.,  4.,  4.])
        >>> kw.sign_slope[rg.active_links]
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.])
        >>> kw.sqrt_slope[rg.active_links]
        array([ 2.,  2.,  1.,  1.,  1.,  2.,  2.])
        >>> kw.slope_sum_at_node
        array([ 0.,  2.,  2.,  0.,  1.,  6.,  6.,  1.,  0.,  2.,  2.,  0.])
        """
        self.slope[self._grid.active_links] = \
            self._grid.calc_grad_at_link(self.elev)[self._grid.active_links]
        self.sqrt_slope = np.sqrt( np.abs(self.slope) )
        self.sign_slope = np.sign( self.slope )
        
        # Calculate the sum of |S|^1/2 for all links at each node. Actually,
        # this should only be done for ACTIVE links. But we've already set it
        # up such that sqrt_slope is zero at inactive links, so we can safely
        # sum over all links. The one subtle catch is that nonexistent links
        # will have ID = -1, meaning we'll be adding in the sqrt_slope for
        # link -1 (i.e., the last node in the grid) in some places. As long as
        # the last numbered link is inactive, that's ok. If for some reason it
        # is active (unlikely), then we'll have problems. Hence the assertion.
        assert (self._grid.status_at_link[-1] == INACTIVE_LINK), \
            'Last numbered link must be inactive'
        self.slope_sum_at_node = np.sum(self.sqrt_slope[self._grid.links_at_node], 1)

    def sort_nodes_by_elevation(self):
        """Create an array of node IDs in descending order of elevation.

        Examples
        --------
        >>> from landlab.components.kinwave.kinwave_semi_implicit import KinwaveSemiImplicit
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((4, 5))
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = np.arange(rg.number_of_nodes)
        >>> kw = KinwaveSemiImplicit(rg)
        >>> kw.core_nodes_by_elev
        array([13, 12, 11,  8,  7,  6])
        """
        nodes_by_elevation = np.argsort(-self.elev)
        (core_nodes, ) = np.where(
            self._grid.status_at_node[nodes_by_elevation] == CORE_NODE)
        self.core_nodes_by_elev = nodes_by_elevation[core_nodes]

    def updated_boundary_conditions(self):
        """Call if boundary conditions are updated.
        """
        self.calc_ground_slope()

    def run_one_step(self, dt, precip_rate=1.0e-5, current_time=0.0, **kwds):
        """Calculate water flow for a time period `dt`.
        """

        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code
            
        # Zero out the inflow discharge for this time step
        self.q_in_at_node[:] = 0.0
        self.disch[:] = 0.0
        
        # TEMP: IF WORKS, MOVE TO FN
        link_mask = np.zeros((self._grid.number_of_nodes, 4))
        for n in self._grid.core_nodes:
            for k in range(link_mask.shape[1]):
                if self.elev[n] > self.elev[self._grid.neighbors_at_node[n,k]]:
                    link_mask[n, k] = 1

        # Loop from high to low elevation
        #
        # Note that for the moment "dx" is used, meaning this works only for
        # grids with uniform spacing. TODO!
        for n in self.core_nodes_by_elev:
            
            # Calculate inflow discharge
            self.q_in_at_node[n] = np.sum(self.disch[self._grid.links_at_node[n]])

            # Calculate coefficient for the denominator
            print n, self.slope_sum_at_node[n]
            coef = (dt * self.slope_sum_at_node[n]) / (self.vel_coef * self._grid.dx)
            
            # This is the denominator in the equation for new depth at this
            # node
            denom = 1.0 + (coef * (self.depth[n] ** 0.6667))
            
            # Now update the depth at this node
            self.depth[n] += ((dt * self.q_in_at_node[n] / self._grid.dx) + 
                              precip_rate * dt -
                              self.infilt * dt) / denom
            
            # Now, send the water downhill by giving it to links
            links = self._grid.links_at_node
            self.disch[links] += link_mask[n,:] * self.vel_coef * self.depth[n] ** (5./3.) * self.sqrt_slope[self._grid.links_at_node[n]]

            print n, self.q_in_at_node[n], coef, denom, self.depth, self.disch[self._grid.links_at_node[n]]


def main():
    """just for testing/debugging"""
    from landlab import RasterModelGrid

    grid = RasterModelGrid((4, 5))
    z = grid.add_zeros('node', 'topographic__elevation')
    z[:] = 5.0 - grid.x_of_node
    grid.set_closed_boundaries_at_grid_edges(False, True, True, True)
    kw = KinwaveSemiImplicit(grid)
    kw.run_one_step(1.0)

if __name__ == '__main__':
    main()