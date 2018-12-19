# -*- coding: utf-8 -*-
"""
Created on Nov 9, 2018

Base of code was copied from DepthDependentDiffuser, written by RC Glade
Also grabbed things from Space, written by CM Shobe

@author: Nicole Gasparini and Joel Johnson
"""

from landlab import Component
import numpy as np
from landlab import INACTIVE_LINK

class EvolvingThresholdFluvialEvolver(Component):

    """
    This component implements an evolving threshold for gravel entrainment, as 
    described by Joel Johnson, 2016, ESURF. 
    Implements Wong & Parker (2006) sediment transport equation:
        q_s* = 3.97(taustar - taustar_crit)^1.5 for taustar > taustar_crit
        (above is equation 1 in Johnson, 2016)
        and
        d(theta_s)/dt = function (depth of sediment deposited or eroded)
        d(theta_s)/dt = k Ar B(d(theta_s)/dx)^k_ent (for entrainment)
        d(theta_s)/dt = -k (1-B) (d(theta_s)/dx)^k_dep (for deposition)
        (above is equation 10 in Johnson, 2016)
        where 
        q_s* = dimensionless sediment transport rate
        taustar = dimensionless shear stress
        taustar_crit = dimensionless threshold shear stress
        theta_s = thickness of sediment eroded or deposited
        k = coefficient w/ units 1/time
        Ar = 1 (for a single grain size)
        B = (taustar_crit_max - taustar_cqs)/(taustar_crit_max - taustar_crit_min)
        taustar_crit_max = upper bound on taustar_dqs
        taustar_crit_min = lower bound on taustar_dqs
        taustar_cqs = evolving threshold shear stress
        k_ent, k_dep = dimensionless exponents
       
    Some notes:
        - In Wong & Parker (2006) we will assume that 
        taustar_crit = taustar_cqs
        - This model implements only one grain size.
        - At nodes that only drain themselves (i.e. upstream-most points), the 
        incoming sediment flux is zero (FOR NOW).
        
    Parameters
    ----------
    grid: ModelGrid
        Landlab ModelGrid object
    rho_s: float
        sediment density, (kg/m^3)
        default = 2600 kg/m^3
    rho_w: float
        water density, (kg/m^3)
        default = 1000 kg/m^3
    g: float
        gravity, (m/s^2)
        default = 9.81 m/s^2
    D: float
       grain diameter, (m)
        default = 0.05 m 
    taustar_crit_max: float
        upper bound on taustar_dqs, (dimensionless)
        default = 0.35
    taustar_crit_min: float
        lower bound on taustar_dqs, (dimensionless)
        default = 0.02
    k: float
        coefficient in crticial shear stress equation, (1/s)
        default = 2.8e-6 1/s
    lambda_p: float
        bed porosity, (dimensionless)
        default = 0.25
    k_ent: float
        exponent on calcuation of critical shear stress if entrainment, (dless)
        default = 0.4
    k_dep: float
        exponent on calcuation of critical shear stress if deposition, (dless)
        default = 0.2

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeatherer
    >>> from landlab.components import DepthDependentDiffuser
    >>> mg = RasterModelGrid((5, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(mg.at_node['soil_production__rate'][mg.core_nodes], 1.)
    True
    >>> DDdiff.soilflux(2.)
    >>> np.allclose(mg.at_node['topographic__elevation'][mg.core_nodes], 0.)
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'][mg.core_nodes], -2.)
    True
    >>> np.allclose(mg.at_node['soil__depth'][mg.core_nodes], 2.)
    True

    """

    _name = 'EvolvingThresholdFluvialEvolver'

    _input_var_names = (
        'flow__receiver_node',
        'flow__upstream_node_order',
        'topographic__steepest_slope',
        'drainage_area',
        'soil__depth'
    )

    _output_var_names = (
        'topographic__elevation'
        'soil__depth'
    )

    _var_units = {
        'flow__receiver_node': '-',
        'flow__upstream_node_order': '-',
        'topographic__steepest_slope': '-',
        'drainage_area': 'm**2',
        'soil__depth': 'm',
        'topographic__elevation': 'm',
    }

    _var_mapping = {
        'flow__receiver_node': 'node',
        'flow__upstream_node_order': 'node',
        'topographic__steepest_slope': 'node',
        'drainage_area': 'node',
        'soil__depth': 'node',
        'topographic__elevation': 'node',
    }

    _var_doc = {
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'topographic__steepest_slope':
            'Topographic slope at each node',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'soil__depth':
            'Depth of sediment above bedrock',
        'topographic__elevation':
            'Land surface topographic elevation',
    }
 
    # Johnson 2016 for now       
    _cite_as = """@Article{esurf-4-685-2-16,
                  AUTHOR = {Johnson, J. P. L.},
                  TITLE = {Gravel threshold of motion: a state function of sediment transport disequilibrium},
                  JOURNAL = {Esurf},
                  VOLUME = {4},
                  YEAR = {2016},
                  PAGES = {685--703},
                  URL = {www.earth-surf-dynam.net/4/685/2016/},
                  DOI = {10.5194/esurf-4-685-2016}
                  }"""        
        
        
    def __init__(self,grid, rho_s=2600, rho_w=1000, g=9.81, D=0.05, 
                 taustar_crit_max=0.35, taustar_crit_min=0.02,
                 k=2.8e-6, lambda_p=0.25, k_ent=0.4, k_dep=0.2, 
                 discharge_field='surface_water__discharge', **kwds):
        """Initialize the EvolvingThresholdFluvialEvolver.
        """
        if (grid.at_node['flow__receiver_node'].size != grid.size('node')):
            msg = ('A route-to-multiple flow director has been '
                   'run on this grid. This will not work with the ' 
                   'EvolvingThresholdFluvialEvolver.')
            raise NotImplementedError(msg)
            
        # Store grid and parameters
        self._grid = grid
        self.p_s = rho_s
        self.p = rho_w
        self.g = g
        self.D = D
        self.A_r = 1.
        self.ts_c_max = taustar_crit_max
        self.ts_c_min = taustar_crit_min
        self.k = k
        self.lambda_p = lambda_p
        self.k_ent = k_ent
        self.k_dep = d_dep

        # create fields
        # elevation
        if 'topographic__elevation' in self.grid.at_node:
            self.elev = self.grid.at_node['topographic__elevation']
        else:
            self.elev = self.grid.add_zeros('node', 'topographic__elevation')

        # slope
        if 'topographic__slope' in self.grid.at_link:
            self.slope = self.grid.at_link['topographic__slope']
        else:
            self.slope = self.grid.add_zeros('link', 'topographic__slope')

        # soil depth
        if 'soil__depth' in self.grid.at_node:
            self.depth = self.grid.at_node['soil__depth']
        else:
            self.depth = self.grid.add_zeros('node', 'soil__depth')

        # soil flux
        if 'soil__flux' in self.grid.at_link:
            self.flux = self.grid.at_link['soil__flux']
        else:
            self.flux=self.grid.add_zeros('link', 'soil__flux')

        # bedrock elevation
        if 'bedrock__elevation' in self.grid.at_node:
            self.bedrock = self.grid.at_node['bedrock__elevation']
        else:
            self.bedrock = self.grid.add_zeros('node', 'bedrock__elevation')
            self.bedrock__elevation[:] = self.topographic__elevation -\
                self.soil__depth
                
        # tau_crit_star as a function of q and s
        self.ts_cqs = np.zeros(grid.number_of_link)
        # tau_star
        self.ts = np.zeros(grid.number_of_link)
        # tau
        #self.t = np.zeros(grid.number_of_link)
        # qs_star
        self.qs_star = np.zeros(grid.number_of_link)
        # qs
        self.qs = np.zeros(grid.number_of_link)
        # dz/dt
        self.dzdt = np.zeros(grid.number_of_nodes)
        # dts_cqs/dt
        self.dts_cqsdt = np.zeros(grid.number_of_nodes)

    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.

        Parameters
        ----------

        dt: float (time)
            The imposed timestep.
        """

        #update soil thickness
        self.grid.at_node['soil__depth'][:] = (
            self.grid.at_node['topographic__elevation']
            - self.grid.at_node['bedrock__elevation'])

        #Calculate soil depth at links.
        H_link = self.grid.map_value_at_max_node_to_link(
            'topographic__elevation','soil__depth')

        #Calculate gradients
        slope = self.grid.calc_grad_at_link(self.elev)
        slope[self.grid.status_at_link == INACTIVE_LINK] = 0.

        #Calculate flux
        self.flux[:] = (-self.K
                        * slope
                        * (1.0 - np.exp(-H_link
                                        / self.soil_transport_decay_depth)))

        #Calculate flux divergence
        dqdx = self.grid.calc_flux_div_at_node(self.flux)

        #Calculate change in soil depth
        dhdt = self.soil_prod_rate - dqdx

        #Calculate soil depth at nodes
        self.depth[self.grid.core_nodes] += dhdt[self.grid.core_nodes] * dt

        #prevent negative soil thickness
        self.depth[self.depth < 0.0] = 0.0

        #Calculate bedrock elevation
        self.bedrock[self.grid.core_nodes] -= (
            self.soil_prod_rate[self.grid.core_nodes] * dt)

        #Update topography
        self.elev[self.grid.core_nodes] = (self.depth[self.grid.core_nodes]
                                         + self.bedrock[self.grid.core_nodes])

    def run_one_step(self, dt):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        self.soilflux(dt)
