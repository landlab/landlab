import numpy as np

from landlab import Component

# =============================================================================
# Defining Classes
# =============================================================================


class WildfireGenerator(Component):
    """
    Simulate stochastic wildfires of different sizes and severity, and
    calculate post-fire vegetation recovery based on aridity.

    Landlab component that simulates landscape-scale fire
    activity driven by global climate-vegetation interactions, following
    the theory of Pausas and Paula (2012) and Pausas and Ribeiro (2013).

    **Fire ignition and spread**

    At each timestep, a Poisson-distributed number of ignition attempts is
    drawn from a mean rate ``potential_fires * dt``. For each attempt, a
    random core node is selected as the ignition point. Ignition succeeds if
    a uniform random number falls below a weighted combination of local fuel
    availability (``fuel_availability``) and aridity, where the relative
    weights are determined by the aridity bin (``_aridity_weights``). River
    channels, nodes where drainage area exceeds ``minimum_river_threshold``,
    act as firebreaks and block both ignition and spread. Setting
    ``minimum_river_threshold`` sufficiently high (e.g. ``np.inf``) effectively
    disables river firebreaks entirely, allowing fire to spread across the full
    grid regardless of drainage network structure. Note that river firebreaks
    are treated as absolute barriers; fire cannot cross nodes where drainage area
    exceeds ``minimum_river_threshold`` regardless of fire intensity or wind.
    This is a simplifying assumption that may not hold for large, wind-driven
    fires. Making firebreaks probabilistic could be tested in future versions.

    Fire spreads iteratively outward from the ignition node to orthogonal
    (D4) neighbors using the same fuel-aridity probability, modulated by a
    logistic slope factor following Rothermel-style upslope enhancement:

    .. math::

        f_{slope} = \\frac{1}{1 + e^{-\\alpha \\cdot S}}

    where :math:`S` is the slope between node and neighbour and :math:`\\alpha`
    is the ``upslope_preference`` parameter. Spread stops when no new nodes
    are added to the fire front.

    **Fire severity**

    Severity is calculated as a power-law function of a randomly perturbed
    aridity index following Grunig et al. (2023):

    .. math::

        \\text{severity} = \\text{aridity\\_index}^{\\gamma}

    where :math:`\\gamma` is ``severity_exponent`` and the aridity index is
    drawn from a normal distribution centred on ``aridity`` (std = 0.05),
    clamped to [0, 1]. Vegetation at burned nodes is reduced by the severity
    factor each fire event.

    **Vegetation recovery**

    After each timestep, vegetation recovers asymptotically toward
    ``max_vegetation`` following:

    .. math::

        \\Delta v = 1 - e^{-e \\cdot dt / T}

    where :math:`T` is the recovery timescale and :math:`e` is the recovery
    exponent, both determined by the aridity bin (``_regrowth_table``).

    **Fire logging**

    Every successful fire event is recorded in ``_fire_records`` and
    accessible via the ``fire_log`` property as a list of dicts containing
    year, ignition node, burned nodes, fire size, severity, aridity, and
    percentage of vegetation removed.


    Authors:
    Matheus de Almeida & Benjamin Campforts

    See the publication:

    de Almeida M., Shobe C.M., Roda-Boluda D.C., Gourbet L., Veraverbeke S.,
    Distelbrink A., Campforts B. (2026) FireLands 1.0: A landscape evolution
    model for simulating the effects of fires and post-fire erosion on sediment
    dynamics in evolving landscapes. Geosci Model Dev:


    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import WildfireGenerator

    >>> np.random.seed(5000)

    >>> dt = 1
    >>> mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> da = mg.add_zeros("drainage_area", at="node")
    >>> fuel = mg.add_zeros("fuel_availability", at="node")
    >>> fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75

    >>> wg = WildfireGenerator(mg)
    >>> wg.run_one_step(dt)
    >>> len(wg.fire_log)
    17


    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    de Almeida M., Shobe C.M., Roda-Boluda D.C., Gourbet L., Veraverbeke S.,
    Distelbrink A., Campforts B. (2026) FireLands 1.0: A landscape evolution
    model for simulating the effects of fires and post-fire erosion on sediment
    dynamics in evolving landscapes. Geosci Model Dev:


    **Additional References**
         Pausas, J. G., & Paula, S. (2012). Fuel shapes the fire-climate
         relationship: evidence from Mediterranean ecosystems. Global
         Ecology and Biogeography, 21(11), 1074-1082.
         https://doi.org/10.1111/j.1466-8238.2012.00769.x

         Pausas JG, Ribeiro E. 2013. The global fire-productivity
         relationship. Global Ecology and Biogeography 22: 728-736.
         https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.12043

         Grunig, M., Seidl, R., & Senf, C. (2023). Increasing aridity
         causes larger and more severe forest fires across Europe.
         Global Change Biology, 29(6), 1648-1659.
         https://doi.org/10.1111/gcb.16547

         Viedma, O., Melia, J., Segarra, D., & Garcia-Haro, J. (1997).
         Modeling rates of ecosystem recovery after fires by using
         landsat TM data. Remote Sensing of Environment, 61(3),
         383-398. https://doi.org/10.1016/S0034-4257(97)00048-5


    """

    _name = "WildfireGenerator"

    _unit_agnostic = True

    _cite_as = """
    @Article{firelandsmodel2026,
        AUTHOR = {de Almeida, M. and Shobe, C. M. and Roda-Boluda, D. C. and
                Gourbet, L. and Veraverbeke, S. and Distelbrink, A. and
                Campforts, B.},
        TITLE = {FireLands 1.0: A landscape evolution model for simulating
                the effects of fires and post-fire erosion on sediment
                dynamics in evolving landscapes},
        JOURNAL = {Geoscientific Model Development},
        YEAR = {2026},
        DOI = {}
    }"""

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "fuel_availability": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Fuel available above the land surface ",
        },
        "last_fire_time": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "y",
            "mapping": "node",
            "doc": "Time of the most recent fire at each node",
        },
    }

    def __init__(
        self,
        grid,
        potential_fires=100,
        max_vegetation=1,
        minimum_river_threshold=2e5,
        upslope_preference=0.3,
        aridity=0.5,
        severity_exponent=0.64,
        seed=5000,
    ):

        super().__init__(grid)

        if seed is not None:
            np.random.seed(seed)

        # --- aridity ---
        if aridity < 0.0:
            raise ValueError("aridity must be >= 0.0")
        if aridity > 1.0:
            raise ValueError("aridity must be <= 1.0")

        # --- potential_fires ---
        if potential_fires <= 0:
            raise ValueError("potential_fires must be > 0")

        # --- max_vegetation ---
        if max_vegetation <= 0.0:
            raise ValueError("max_vegetation must be > 0.0")

        # --- minimum_river_threshold ---
        if minimum_river_threshold < 0.0:
            raise ValueError("minimum_river_threshold must be >= 0.0")

        # --- upslope_preference ---
        if upslope_preference < 0.0:
            raise ValueError(
                "upslope_preference (alpha) must be >= 0.0; "
                "negative values would make downslope spread more likely than upslope"
            )

        # --- severity_exponent ---
        if severity_exponent <= 0.0:
            raise ValueError("severity_exponent must be > 0.0")
        if severity_exponent > 1.0:
            raise ValueError(
                "severity_exponent must be <= 1.0; "
                "values above 1.0 produce severity > aridity, which has no physical basis"
            )

        # --- fuel_availability field values ---
        fuel = grid.at_node["fuel_availability"]
        if np.any(fuel < 0.0):
            raise ValueError("fuel_availability cannot contain negative values")
        if np.any(fuel > max_vegetation):
            raise ValueError(
                "fuel_availability contains values exceeding max_vegetation "
                f"({max_vegetation}); initial fuel load cannot exceed the maximum"
            )

        # parameters (private)
        self._potential_fires = potential_fires
        self._max_vegetation = max_vegetation
        self._minimum_river_threshold = minimum_river_threshold
        self._alpha = upslope_preference
        self._aridity = aridity
        self._sev_exponent = severity_exponent
        self._current_time = None

        # grid field references
        self._vegetation = grid.at_node["fuel_availability"]
        self._drainage_area = grid.at_node["drainage_area"]
        self._topo = grid.at_node["topographic__elevation"]

        # internal state
        self._fire_records = []

        # initialize output fields
        self.initialize_output_fields()
        self._grid.at_node["last_fire_time"][
            :
        ] = -np.inf  # Initialize last_fire_time to -inf (no fire has occurred yet)
        self._last_fire_time = self._grid.at_node["last_fire_time"]

        self._regrowth_table = [
            (
                0.0,
                0.1,
                5.4,
                0.42,
            ),  # Trees + shrubs, data derived from Viedma et al. 1997
            (0.1, 0.3, 8.0, 0.29),  # Dense shrubs
            (0.3, 0.5, 12.4, 0.18),  # Dense shrubs
            (0.5, 0.7, 15.0, 0.15),  # Sparse shrubs
            (0.7, 1.0, 18.5, 0.12),  # Sparse shrubs
        ]

        self._aridity_weights = [
            (
                0.0,
                0.1,
                0.5,
                0.5,
            ),  # Trees + shrubs, data derived from Viedma et al. 1997
            (0.1, 0.3, 0.5, 0.5),  # Dense shrubs
            (0.3, 0.5, 0.5, 0.5),  # Dense shrubs
            (0.5, 0.7, 0.8, 0.2),  # Sparse shrubs
            (0.7, 1.0, 0.9, 0.1),  # Sparse shrubs
        ]

    @property
    def fire_sizes(self):
        """Fire sizes (km²) for all recorded fires."""
        return [r["fire_size (km2)"] for r in self._fire_records]

    @property
    def burned_nodes(self):
        """Burned node lists for all recorded fires."""
        return [r["burned_nodes"] for r in self._fire_records]
    
    @property
    def ignition_nodes(self):
        """lists node lists for all recorded fires."""
        return [r["ignition_node"] for r in self._fire_records]

    @property
    def severity(self):
        """Severity factors for all recorded fires."""
        return [r["severity_factor"] for r in self._fire_records]

    @property
    def _rivers(self):
        return self._drainage_area > self._minimum_river_threshold

    @property
    def fire_log(self):
        """Fire records as a list of dicts."""
        return self._fire_records

    @property
    def potential_fires(self):
        """Mean number of ignition attempts per unit time [-]."""
        return self._potential_fires

    @property
    def upslope_preference(self):
        """Slope factor exponent controlling upslope fire spread bias [-]."""
        return self._alpha

    @property
    def aridity(self):
        """Aridity index controlling fire severity and vegetation recovery [-]."""
        return self._aridity

    @property
    def last_step_fires(self):
        """Records of all fires ignited in the most recent timestep."""
        return [r for r in self._fire_records if r["year"] == self._current_time]

    @aridity.setter
    def aridity(self, new_val):
        if new_val < 0.0:
            raise ValueError("aridity must be >= 0.0")
        if new_val > 1.0:
            raise ValueError("aridity must be <= 1.0")
        self._aridity = new_val

    def _get_regrowth_params(self):
        for low, high, t, e in self._regrowth_table:
            if low <= self._aridity < high or (self._aridity == 1 and high == 1.0):
                return t, e
        raise ValueError(f"Aridity {self._aridity} is out of range [0,1].")

    def _get_aridity_weights(self):
        for low, high, f_weight, a_weight in self._aridity_weights:
            if low <= self._aridity < high or (self._aridity == 1 and high == 1.0):
                return f_weight, a_weight
        raise ValueError(f"Aridity {self._aridity} is out of range [0,1].")

    def _regrow_vegetation(self, dt):
        regrowth_time, regrowth_exponent = self._get_regrowth_params()

        delta = 1 - np.exp(-regrowth_exponent * dt / regrowth_time)

        # Update the vegetation grid (regrow)
        self._vegetation[:] = np.minimum(self._vegetation + delta, self._max_vegetation)

    def _get_neighbors(self, node):
        """Method to get the active neighbor nodes (for fire spreading)"""
        return self._grid.active_adjacent_nodes_at_node[node]

    def _calc_severity(self):
        """Calculate fire severity as a power-law function of aridity.
        Draws a random aridity index from a normal distribution centred on
        the component aridity (std=0.05), clamps it to [0, 1], then applies
        the power-law severity = aridity_index ** sev_exponent following
        Grunig et al. (2023). The result is clamped to a maximum of 1 and
        rounded to 2 decimal places.
        """
        aridity_index = np.clip(np.random.normal(self._aridity, 0.05), 0, 1)
        return np.round(np.minimum(np.power(aridity_index, self._sev_exponent), 1), 2)

    def _fire(self, dt):
        """Simulate fire ignition and spreading throughout the grid.
        Draws a Poisson number of ignition attempts. For each attempt, a random
        core node is selected; ignition succeeds if a uniform random number falls
        below a weighted combination of local fuel availability and aridity. Fire
        spreads iteratively to orthogonal neighbors using a logistic slope factor
        and the same fuel-aridity spread probability, stopping at river firebreaks
        (drainage area > minimum_river_threshold). After spreading, vegetation is
        reduced by the severity factor and the event is logged in _fire_records.
        """
        fire_ignitions = np.random.poisson(self._potential_fires * dt)
        if fire_ignitions == 0:
            return

        cell_area = self._grid.cell_area_at_node

        for _ in range(fire_ignitions):
            fuel_weight, aridity_weight = self._get_aridity_weights()
            center = np.random.choice(self._grid.core_nodes)

            if self._rivers[center]:
                continue

            ignition_prob = (
                self._vegetation[center] * fuel_weight + self._aridity * aridity_weight
            )
            if np.random.rand() > ignition_prob:
                continue

            severity_factor = self._calc_severity()
            fire_front = {center}
            burned = set()

            while fire_front:
                new_front = set()
                for node in fire_front:
                    if node in burned:
                        continue
                    burned.add(node)
                    for neighbor in self._get_neighbors(node):
                        if (
                            neighbor in burned
                            or neighbor == -1
                            or self._drainage_area[neighbor]
                            > self._minimum_river_threshold
                        ):
                            continue
                        dx = self._grid.x_of_node[neighbor] - self._grid.x_of_node[node]
                        dy = self._grid.y_of_node[neighbor] - self._grid.y_of_node[node]
                        slope = (self._topo[neighbor] - self._topo[node]) / np.sqrt(
                            dx**2 + dy**2
                        )
                        slope_factor = 1 / (1 + np.exp(-self._alpha * slope))
                        spread_prob = (
                            self._vegetation[neighbor] * fuel_weight
                            + self._aridity * aridity_weight
                        ) * slope_factor
                        if np.random.rand() < spread_prob:
                            new_front.add(neighbor)
                fire_front = new_front

            if not burned:
                continue

            changed_nodes = list(burned)
            veg_before = self._vegetation[changed_nodes].copy()
            self._vegetation[changed_nodes] *= 1 - severity_factor
            non_zero_veg = veg_before > 0
            veg_change = (
                1
                - np.mean(
                    self._vegetation[changed_nodes][non_zero_veg]
                    / veg_before[non_zero_veg]
                )
            ) * 100
            self._last_fire_time[changed_nodes] = self._current_time

            self._fire_records.append(
                {
                    "year": self._current_time,
                    "ignition_node": self._grid.xy_of_node[center],
                    "burned_nodes": changed_nodes,
                    "fire_size (km2)": cell_area[changed_nodes].sum() / 1e6,
                    "severity_factor": severity_factor,
                    "aridity": self._aridity,
                    "pct of vegetation removed": veg_change,
                }
            )

    def run_one_step(self, dt, current_time=None):

        if current_time is not None:
            self._current_time = current_time
        else:
            if self._current_time is None:
                self._current_time = dt
            else:
                self._current_time += dt
        self._fire(dt)
        self._regrow_vegetation(dt)
