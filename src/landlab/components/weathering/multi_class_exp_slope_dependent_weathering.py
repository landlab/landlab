#!/usr/bin/env python
"""

@author: YShmilovitz
"""

import numpy as np
from landlab import Component, FieldError

class MultiClassExpSlopeDependentWeatherer(Component):
    """Calculate exponential and slope-dependent weathering of bedrock.
    This component spreads the weathered material into various sediment grain classes.

    Exponential soil production is calculated based on Ahnert (1976):
    Consider that ``w0`` is the maximum soil production rate and
    that ``w_star`` is the characteristic soil production depth. The
    soil production rate is given by:

        soil_production =  w0 * exp(-soil__depth / w_star)

    Slope dependent weathering is calculated according to local topographic slope ``S`` and
    maximum weathering rate ``kappa``:

        slope_dependent_weathering = kappa * S * exp(-soil__depth / w_star)

    References
    ----------
    Ahnert, F. (1976). Brief description of a comprehensive three-dimensional
    process-response model of landform development Z. Geomorphol. Suppl.  25,
    29 - 49.


    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.weathering import MultiClassExpSlopeDependentWeatherer
    >>> grid = RasterModelGrid((3, 4), xy_spacing=10.0)
    >>> phi, rho = 0.4, 2650.0
    >>> n_classes = 2
    >>> grid.at_node["grains__weight"] = np.full((grid.number_of_nodes, n_classes), 5.0)
    >>> grid.at_node["bed_grains__proportions"] = np.full(
    ...     (grid.number_of_nodes, n_classes), 1.0 / n_classes
    ... )

    Soil depth must match the grain mass: z = M / (rho * (1 - phi)).

    >>> soil = grid.add_field(
    ...     "soil__depth",
    ...     grid.at_node["grains__weight"].sum(axis=1) / (rho * (1 - phi)),
    ...     at="node",
    ... )
    >>> bedrock = grid.add_zeros("bedrock__elevation", at="node")
    >>> topo = grid.add_field("topographic__elevation", bedrock + soil, at="node")
    >>> slope = grid.add_zeros("topographic__steepest_slope", at="node")
    >>> slope[grid.core_nodes] = [0.0, 0.1]
    >>> weatherer = MultiClassExpSlopeDependentWeatherer(
    ...     grid,
    ...     phi=phi,
    ...     rho=rho,
    ...     soil_production_maximum_rate=1e-2,
    ...     slope_dependent_weathering_maximum_rate=1e-2,
    ...     soil_production_decay_depth=1.0,
    ... )
    >>> weatherer.run_one_step(dt=1.0)
    >>> np.round(grid.at_node["soil_production__rate"][grid.core_nodes], 5)
    array([0.00994, 0.01093])
    """

    _name = "MultiClassExpSlopeDependentWeatherer"

    _unit_agnostic = True

    _cite_as = """
    """

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Topographic elevation at node",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Bedrock elevation at node",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
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
        "grains__weight": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "kg/m2",
            "mapping": "node",
            "doc": "Mass per unit area of grains in each size class stored at node",
        },
        "bed_grains__proportions": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Proportional mass of each grain size class in the bed layer",
        },
    }

    def __init__(
        self, grid,
            phi=0.35,
            rho=2650,
            soil_production_maximum_rate=10**-2,
            slope_dependent_weathering_maximum_rate=10**-4,
            soil_production_decay_depth=1.0,
            slope_name = 'topographic__steepest_slope',
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
        soil_production_maximum_rate : float (default 10**-2), array of float [m/yr]
            Maximum weathering rate for bare bedrock with no dependecy on the
            topographic slope
        slope_dependent_weathering_maximum_rate : float (default 10**-4), array of float [m/yr]
            Maximum slope-depedent weathering rate for bare bedrock
        soil_production_decay_depth : float (default 1), array of float [m]
            Characteristic weathering depth
        slope_name : string (default 'topographic__steepest_slope')
            The name of the slope field used
        """
        super().__init__(grid)

        self._slope = self._get_slope_field(slope_name)
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

        self.calc_rock_exposure_fraction()
        self.grid.at_node["soil_production__rate"][self._cores] = ((self._p0+self._kappa*self._slope[self._cores]) *
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


    def _get_slope_field(self, slope):
        """Return the slope array at nodes, checking that it exists and is valid."""
        if isinstance(slope, str):
            if slope not in self._grid.at_node:
                raise FieldError(
                    f"Slope field '{slope}' not found at nodes. "
                )
            slope_array = self._grid.at_node[slope]
        else:
            slope_array = np.asarray(slope, dtype=float)
            if slope_array.shape != (self._grid.number_of_nodes,):
                raise ValueError(
                    "Slope array must have one value per node: expected shape "
                    f"({self._grid.number_of_nodes},), got {slope_array.shape}."
                )

        if not np.all(np.isfinite(slope_array[self._grid.core_nodes])):
            raise ValueError("Slope field contains NaN or inf values at core nodes.")

        return slope_array


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

