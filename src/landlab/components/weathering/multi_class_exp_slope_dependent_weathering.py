#!/usr/bin/env python
"""

@author: YShmilovitz
"""

import numpy as np

from landlab import Component


class MultiClassExpSlopeDependentWeatherer(Component):
    """Calculate exponential and slope-dependent weathering of bedrock taking
    into account various sediment grains classes.

    Uses exponential soil production function in the style of Ahnert (1976).

    Consider that ``w0`` is the maximum soil production rate and
    that ``w_star`` is the characteristic soil production depth. The
    soil production rate ``w0`` is given as a function of the soil
    depth::

        soil_production =  w0 * exp(-soil__depth / w_star)

    Slope dependent weathering is calculate according to local topographic slope ``S`` and
    maximum weathering rate ``kappa``:

        slope_dependent_weathering = kappa * S
    
    References
    ----------
    Ahnert, F. (1976). Brief description of a comprehensive three-dimensional
    process-response model of landform development Z. Geomorphol. Suppl.  25,
    29 - 49.

    Armstrong, A. (1976). A three dimensional simulation of slope forms.
    Zeitschrift für Geomorphologie  25, 20 - 28.

    # Examples
    # --------
    # >>> from landlab.components.weathering import MultiClassExpSlopeDependentWeatherer
    # >>> from landlab import RasterModelGrid
    # >>> from landlab.components.soil_grading import SoilGrading
    # >>> from landlab.components import FlowAccumulator
    # >>> rows=3
    # >>> columns=3
    # >>> xy_spacing = 10
    # >>> porosity = 0.4
    # >>> sed_depth = 500
    # >>> rho_sed = 2650
    # >>> grain_sizes = 0.01
    # >>> grains_weight = 10
    # >>> p0=10**-2
    # >>> kappa=10**-4
    # >>> grid = RasterModelGrid((rows, columns), xy_spacing=xy_spacing)
    # >>> sg = SoilGrading(grid,
    # ...     meansizes=grain_sizes,
    # ...     grains_weight=grains_weight,
    # ...     phi=porosity,
    # ...     soil_density=rho_sed)
    # >>> weatherer = MultiClassExpSlopeDependentWeatherer(grid,
    # ...     phi=porosity,
    # ...     rho=rho_sed,
    # ...     soil_production_maximum_rate=p0,
    # ...     slope_dependent_weathering_maximum_rate=kappa)
    # >>> fa = FlowAccumulator(grid)
    # >>> fa.run_one_step()
    # >>> weatherer.run_one_step(dt=1)
    # >>> np.round(grid.at_node['soil__depth'][grid.core_nodes][0],4)
    np.float64(0.0162)
    """

    _name = "MultiClassExpSlopeDependentWeatherer"

    _unit_agnostic = True

    _cite_as = """
    """

    _info = {
        "soil__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "soil_production__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": "rate of soil production at nodes",
        },
    }

    def __init__(
        self, grid,
            phi=0.35,
            rho=2650,
            soil_production_maximum_rate=10**-2,
            slope_dependent_weathering_maximum_rate=10**-4,
            soil_production_decay_depth=1.0,
    ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        phi : float  (default 0.35)
            Bulk porosity of bed sediment
        rho : float (default 2650)
            Density of the sediment layer
        soil_production_maximum_rate : float, array of float
            Maximum weathering rate for bare bedrock with no dependecy on the
            topographic slope
        slope_dependent_weathering_maximum_rate : float, array of float
            Maximum slope-depedent weathering rate for bare bedrock
        soil_production_decay_depth : float, array of float
            Characteristic weathering depth
        """
        super().__init__(grid)

        self._phi = phi
        self._rho = rho
        self._p0 = soil_production_maximum_rate
        self._kappa = slope_dependent_weathering_maximum_rate
        self._depth_decay_scale = soil_production_decay_depth
        self._cores = self._grid.core_nodes

        # Create soil_production__rate field
        if "soil_production__rate" not in grid.at_node:
            grid.add_zeros("soil_production__rate", at="node")

    def calc_soil_prod_rate(self):
        """Calculate soil production rate."""

        slope = self._grid.at_node['topographic__steepest_slope']

        self.calc_rock_exposure_fraction()
        self.grid.at_node["soil_production__rate"][self._cores] = ((self._p0+self._kappa*slope[self._cores]) *
                         self._rock_exposure_fraction[self._cores])


    def calc_rock_exposure_fraction(self):
        """Update the bedrock exposure fraction.
        """
        self._rock_exposure_fraction = np.exp(-self._grid.at_node['soil__depth']/ self._depth_decay_scale)


    def _update_dz_and_mass(self,dt):

        # Pointers
        soil = self._grid.at_node['soil__depth']
        bedrock = self._grid.at_node['bedrock__elevation']
        grains__weight = self._grid.at_node['grains__weight']
        slope = self._grid.at_node['topographic__steepest_slope']
        topo = self._grid.at_node['topographic__elevation']

        # Convert dz/dt to dz
        dz_weathering = self.grid.at_node["soil_production__rate"][self._cores]*dt
        # Conevrt dz to mass per area [kg/m2]
        dmass_per_area = dz_weathering * self._rho * (1 - self._phi)

        # Update fields
        if np.ndim(self._grid.at_node['bed_grains__proportions']) > 1:
            grains__weight[self._cores, :] += dmass_per_area[:, np.newaxis] * \
                                                        self._grid.at_node['bed_grains__proportions'][
                                                            self._cores]
        else:
            grains__weight[self._cores] += dmass_per_area * self._grid.at_node['bed_grains__proportions'][
                self._cores]

        bedrock[self._cores] -= dz_weathering
        soil[self._cores] += dz_weathering
        topo[self._cores] = soil[self._cores] + bedrock[self._cores]

    def run_one_step(self,
                     dt=1):

        """
        Parameters
        ----------
        dt: float
            model time step (yr)
        """
        self.calc_soil_prod_rate()
        self._update_dz_and_mass(dt=dt)

    @property
    def maximum_soil_production_rate(self):
        """Maximum rate of weathering with no depdency on the topographic slope (m/yr)."""
        return self._p0


    @property
    def maximum_slope_dependent_weathering_rate(self):
        """Slope-dependent maximum weathering rate (m/yr)."""
        return self._kappa
