import numpy as np
from requireit import require_between
from requireit import require_positive

from landlab import Component

from .fire_event_log import FireEvent

# Regrowth parameters by aridity bin, derived from Viedma et al. (1997).
# Columns: (aridity_low, aridity_high, regrowth_time, regrowth_exponent)
_REGROWTH_TABLE = (
    (0.0, 0.1, 5.4, 0.42),  # Trees + shrubs
    (0.1, 0.3, 8.0, 0.29),  # Dense shrubs
    (0.3, 0.5, 12.4, 0.18),  # Dense shrubs
    (0.5, 0.7, 15.0, 0.15),  # Sparse shrubs
    (0.7, 1.0, 18.5, 0.12),  # Sparse shrubs
)

# Ignition/spread weights by aridity bin, derived from Viedma et al. (1997).
# Columns: (aridity_low, aridity_high, fuel_weight, aridity_weight)
_ARIDITY_WEIGHTS = (
    (0.0, 0.1, 0.5, 0.5),  # Trees + shrubs
    (0.1, 0.3, 0.5, 0.5),  # Dense shrubs
    (0.3, 0.5, 0.5, 0.5),  # Dense shrubs
    (0.5, 0.7, 0.8, 0.2),  # Sparse shrubs
    (0.7, 1.0, 0.9, 0.1),  # Sparse shrubs
)


def _get_aridity_bin(aridity, table):
    """Return the two values in a table row for the bin containing aridity.

    Bins are half-open, ``[low, high)``, except the last bin, which also
    includes its upper edge so that the maximum aridity is valid.

    Examples
    --------
    >>> table = ((0.0, 0.5, 1.0, 2.0), (0.5, 1.0, 3.0, 4.0))
    >>> _get_aridity_bin(0.25, table)
    (1.0, 2.0)
    >>> _get_aridity_bin(0.5, table)
    (3.0, 4.0)
    >>> _get_aridity_bin(1.0, table)
    (3.0, 4.0)
    """
    upper_edge = table[-1][1]
    for low, high, a, b in table:
        if low <= aridity < high or aridity == high == upper_edge:
            return a, b
    raise ValueError(f"aridity {aridity} is outside the table range")


def get_regrowth_params(aridity, table=_REGROWTH_TABLE):
    """Return ``(regrowth_time, regrowth_exponent)`` for an aridity value.

    Examples
    --------
    >>> get_regrowth_params(0.2)
    (8.0, 0.29)
    >>> get_regrowth_params(1.0)
    (18.5, 0.12)
    """
    return _get_aridity_bin(aridity, table)


def get_aridity_weights(aridity, table=_ARIDITY_WEIGHTS):
    """Return ``(fuel_weight, aridity_weight)`` for an aridity value.

    Examples
    --------
    >>> get_aridity_weights(0.6)
    (0.8, 0.2)
    """
    return _get_aridity_bin(aridity, table)


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
    logistic slope factor following Rothermel-style upslope enhancement::



        slope_factor = 1 / (1 + exp(-upslope_preference * slope))

    where ``slope`` is the slope between node and neighbour and ``upslope_preference``
    is the ``upslope_preference`` parameter. Spread stops when no new nodes
    are added to the fire front.

    **Fire severity**

    Severity is calculated as a power-law function of a randomly perturbed
    aridity index following Grunig et al. (2023)::



        severity = aridity ** sev_exponent

    where the aridity index is drawn from a normal distribution centred on ``aridity`` (std = 0.05),
    clamped to [0, 1]. Vegetation at burned nodes is reduced by the severity
    factor each fire event.

    **Vegetation recovery**

    After each timestep, vegetation recovers asymptotically toward
    ``max_vegetation`` following::



        delta_vegetation = 1 - exp(-regrowth_exponent * dt / regrowth_time)

    where ``regrowth_time`` is the recovery timescale and ``regrowth_exponent`` is the recovery
    exponent, both determined by the aridity bin (``_regrowth_table``).

    **Fire logging**

    Logging is optional. If ``log`` is a callable, it is called with a
    :class:`FireEvent` for every fire. :class:`FireEventLog` keeps events in
    memory; any other callable, such as ``list.append`` or a CSV writer, also
    works. If ``log`` is ``None`` (the default), nothing is recorded.


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
    >>> from landlab.components.wildfire_generator import FireEventLog

    >>> np.random.seed(5000)

    >>> dt = 1
    >>> mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> da = mg.add_zeros("drainage_area", at="node")
    >>> fuel = mg.add_zeros("fuel_availability", at="node")
    >>> fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75


    >>> log = FireEventLog()
    >>> wg = WildfireGenerator(mg, log=log)
    >>> wg.run_one_step(dt)
    >>> len(log.records)
    13


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
        rng=5000,
        log=None,
    ):

        super().__init__(grid)

        self._rng = np.random.default_rng(rng)

        # --- aridity ---
        self._aridity = require_between(aridity, a_min=0.0, a_max=1.0, name="aridity")

        # --- potential_fires ---
        self._potential_fires = require_positive(
            potential_fires, name="potential_fires"
        )

        # --- max_vegetation ---
        self._max_vegetation = require_positive(max_vegetation, name="max_vegetation")
        # --- minimum_river_threshold ---
        self._minimum_river_threshold = require_positive(
            minimum_river_threshold, name="minimum_river_threshold"
        )
        # --- upslope_preference ---
        self._upslope_preference = require_positive(
            upslope_preference, name="upslope_preference"
        )

        # --- severity_exponent ---
        self._severity_exponent = require_between(
            severity_exponent, a_min=0.0, a_max=1.0, name="severity_exponent"
        )

        if log is not None and not callable(log):
            raise TypeError("log must be callable or None")
        self._log = log

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
        self._current_time = 0.0

        # grid field references
        self._vegetation = grid.at_node["fuel_availability"]
        self._drainage_area = grid.at_node["drainage_area"]
        self._topo = grid.at_node["topographic__elevation"]

        # initialize output fields
        self.initialize_output_fields()
        self._grid.at_node["last_fire_time"][
            :
        ] = -np.inf  # Initialize last_fire_time to -inf (no fire has occurred yet)
        self._last_fire_time = self._grid.at_node["last_fire_time"]

    @property
    def _rivers(self):
        return self._drainage_area > self._minimum_river_threshold

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

    @aridity.setter
    def aridity(self, new_val):
        if new_val < 0.0:
            raise ValueError("aridity must be >= 0.0")
        if new_val > 1.0:
            raise ValueError("aridity must be <= 1.0")
        self._aridity = new_val

    def _regrow_vegetation(self, dt):

        regrowth_time, regrowth_exponent = get_regrowth_params(self._aridity)

        self._vegetation[:] = self._max_vegetation - (
            self._max_vegetation - self._vegetation
        ) * np.exp(-regrowth_exponent * dt / regrowth_time)

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
        aridity_index = np.clip(self._rng.normal(self._aridity, 0.05), 0, 1)
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
        fire_ignitions = self._rng.poisson(self._potential_fires * dt)
        if fire_ignitions == 0:
            return

        cell_area = self._grid.cell_area_at_node

        for _ in range(fire_ignitions):
            fuel_weight, aridity_weight = get_aridity_weights(self._aridity)
            center = self._rng.choice(self._grid.core_nodes)

            if self._rivers[center]:
                continue

            ignition_prob = (
                self._vegetation[center] * fuel_weight + self._aridity * aridity_weight
            )
            if self._rng.random() > ignition_prob:
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
                        if self._rng.random() < spread_prob:
                            new_front.add(neighbor)
                fire_front = new_front

            changed_nodes = list(burned)
            veg_before = self._vegetation[changed_nodes].copy()
            self._vegetation[changed_nodes] *= 1 - severity_factor

            self._last_fire_time[changed_nodes] = self._current_time

            if self._log is not None:
                had_veg = veg_before > 0.0
                if had_veg.any():
                    removed_fraction = 1.0 - np.mean(
                        self._vegetation[changed_nodes][had_veg] / veg_before[had_veg]
                    )
                else:
                    removed_fraction = 0.0

                self._log(
                    FireEvent(
                        year=float(self._current_time),
                        xy_of_ignition=tuple(
                            float(v) for v in self._grid.xy_of_node[center]
                        ),
                        burned_nodes=tuple(int(n) for n in changed_nodes),
                        fire_size=float(cell_area[changed_nodes].sum()),
                        severity_factor=severity_factor,
                        aridity=float(self._aridity),
                        vegetation_removed_fraction=float(removed_fraction),
                    )
                )

    def run_one_step(self, dt):

        self._current_time += dt

        self._fire(dt)
        self._regrow_vegetation(dt)
