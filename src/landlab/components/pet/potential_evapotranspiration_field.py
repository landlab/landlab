import copy

import numpy as np

from landlab import Component
from landlab.grid.mappers import map_node_to_cell

_VALID_METHODS = {"PriestleyTaylor", "PenmanMonteith", "NetRadEqPE"}


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class PotentialEvapotranspiration(Component):
    """
    Potential Evapotranspiration Component calculates spatially distributed
    potential evapotranspiration based on input radiation factor (spatial
    distribution of incoming radiation) using chosen method such as constant
    or Priestley Taylor.
    ASCE-EWRI Task Committee Report Jan 2005 for 'PriestleyTaylor' method.
    Note: Calling 'PriestleyTaylor' method would generate/overwrite shortwave &
    longwave radiation fields.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu and Berkan Mertan

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.pet import PotentialEvapotranspiration

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    >>> grid.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [1.0, 1.0, 1.0, 1.0],
    ...     [2.0, 2.0, 2.0, 2.0],
    ...     [3.0, 4.0, 4.0, 3.0],
    ...     [4.0, 4.0, 4.0, 4.0],
    ... ]
    >>> PET = PotentialEvapotranspiration(grid)
    >>> PET.name
    'PotentialEvapotranspiration'
    >>> sorted(PET.output_var_names)
     ['surface__potential_evapotranspiration_rate']
    >>> sorted(PET.units)
    [('surface__potential_evapotranspiration_rate', 'mm')]
    >>> PET.grid.number_of_node_rows
    5
    >>> PET.grid.number_of_node_columns
    4
    >>> PET.grid is grid
    True
    >>> pet_rate = grid.at_cell["surface__potential_evapotranspiration_rate"]
    >>> np.allclose(pet_rate, 0.0)
    True
    >>> PET._current_time = 0.5
    >>> PET.update()
    >>> np.allclose(pet_rate, 0.0)
    False

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    ASCE-EWRI: The ASCE standardized reference evapotranspiration equation, in:
    Standardization of Reference Evapotranspiration Task Committee Final Report,
    edited by: Allen, R. G., Walter, I. A., Elliot, R. L., Howell, T. A.,
    Itenﬁsu, D., Jensen, M. E., and Snyder, R. L., Technical Committee report
    to the Environmental and Water Resources Institute of the American Society
    of Civil Engineers from the Task Committee on Standardization of Reference
    Evapotranspiration, Reston, VA, USA, 2005.

    Zhou, X., Istanbulluoglu, E., and Vivoni, E. R.: Modeling the
    ecohydrological role of aspect-controlled radiation on tree-grass-shrub
    coexistence in a semiarid climate, Water Resour. Res., 49, 2872– 2895,
    doi:10.1002/wrcr.20259, 2013.

    """

    _name = "PotentialEvapotranspiration"

    _unit_agnostic = False

    _info = {
        "surface__potential_evapotranspiration_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm",
            "mapping": "cell",
            "doc": "potential sum of evaporation and potential transpiration",
        }
    }

    def __init__(
        self,
        grid,
        method="PriestleyTaylor",
        priestley_taylor_const=1.26,
        relative_humidity=0.65,
        albedo=0.6,
        air_density=None,
        latent_heat_of_vaporization=28.34,
        psychometric_const=0.066,
        LAI=2.88,
        stefan_boltzmann_const=0.0000000567,
        solar_const=1366.67,
        latitude=34.0,
        elevation_of_measurement=300,
        adjustment_coeff=0.18,
        lt=0.0,
        nd=365.0,
        MeanTmaxF=12.0,
        delta_d=5.0,
        current_time=0.5,
        Tmin=10,
        Tmax=25,
        Tavg=17.5,
        Rl=100,
        Zm=2.0,
        Zd=0.084,
        Zo=0.012,
        Zveg=None,
        Vwind=3.12,
        temperatures=None,
        radiation=None,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        method: {'PriestleyTaylor', 'PenmanMonteith', 'NetRadEqPE'}, optional
            Priestley Taylor method will spit out radiation outputs too.
        priestley_taylor_constant: float, optional
            Alpha used in Priestley Taylor method.
        relative_humidity: float, optional
            Relative humidity factor used to determine real Priestley Taylor constant
        albedo: float, optional
            Albedo.
        air_density: float, field, optional
            Single value or spatially distributed air densities in kg/m^3.
        latent_heat_of_vaporization: float, optional
            Latent heat of vaporization for water Pwhv (Wd/(m*mm^2)).
        psychometric_const: float, optional
            Psychometric constant (kPa (deg C)^-1).
        LAI: float, field, optional
            LAI used for Penman Monteith equation
        stefan_boltzmann_const: float, optional
            Stefan Boltzmann's constant (W/(m^2K^-4)).
        solar_const: float, optional
            Solar constant (W/m^2).
        latitude: float, optional
            Latitude (radians).
        elevation_of_measurement: float, optional
            Elevation at which measurement was taken (m).
        adjustment_coeff: float, optional
            adjustment coeff to predict Rs from air temperature (deg C)^-0.5.
        lt: float, optional
            lag between peak TmaxF and solar forcing (days).
        nd: float, optional
            Number of days in year (days).
        MeanTmaxF: float, optional
            Mean annual rate of TmaxF (mm/d).
        delta_d: float, optional
            Calibrated difference between max & min daily TmaxF (mm/d).
        current_time: float, required only for 'Cosine' method
            Current time (Years)
        Tmin: float, required for 'Priestley Taylor' method
            Minimum temperature of the day (deg C)
        Tmax: float, required for 'Priestley Taylor' method
            Maximum temperature of the day (deg C)
        Tavg: float, required for 'Priestley Taylor' method
            methods
            Average temperature of the day (deg C)
        Rl: float, optional
            Stomatal wall resistance in s/m, set by default to 100s/m
        Zm: float, optional
            Elevation / height at which wind speed is recorded
        Zd: float, optional
            Zero plane displacement height
        Zo: float, optional
            Roughness length
        Zveg: float, optional
            Vegetation height, used if other Z variables are default
        Vwind: float, optional
            Wind speed / velocity
        temperatures: float, field, optional
            Spatially distributed temperatures (deg C) to apply to PET grid/fields as a whole
        radiation float, optional
            User-provided net radiation field, if not given net radiation will be
            computed internally (W/m^2)
        """
        super().__init__(grid)

        # Grid copy is used for node to cell mapping operations
        self._gridCopy = copy.deepcopy(self._grid)

        self._current_time = current_time

        self._zm = Zm
        self._zd = Zd
        self._zo = Zo
        self._zveg = Zveg
        self._vz = Vwind
        self._rl = Rl

        self._user_radiation = radiation
        self._method = method
        # For Priestley Taylor
        self._alpha = priestley_taylor_const
        self._a = albedo
        self._pa = air_density
        self._pwhv = latent_heat_of_vaporization
        self._y = psychometric_const
        self._sigma = stefan_boltzmann_const
        self._Gsc = solar_const
        self._latitude = latitude
        self._Krs = adjustment_coeff
        self._LT = lt
        self._LAI = self._validate_lai(LAI)
        self._ND = nd
        self._TmaxF_mean = MeanTmaxF
        self._DeltaD = delta_d
        self._relative_humidity = relative_humidity
        self._temperatures = temperatures
        _assert_method_is_valid(self._method)

        self.initialize_output_fields()

        self.Tmin, self.Tmax = self._validate_temperature_range(Tmin, Tmax)

        # If user provides Tmin/Tmax and not Tavg, calculate Tavg
        if not isinstance(Tavg, np.ndarray) and Tavg == 17.5:
            self.Tavg = (self._Tmin + self._Tmax) / 2
        else:
            self.Tavg = Tavg

        # Import Radiation to instantiate an internal radiation component
        from landlab.components import Radiation

        if "topographic__elevation" not in self._gridCopy.at_node.keys():
            self._gridCopy.add_zeros("topographic__elevation", at="node")

        # Internal radiation component
        self._etpRad = Radiation(
            self._gridCopy,
            method="Grid",
            latitude=self._latitude,
            current_time=self._current_time,
            albedo=self._a,
            max_daily_temp=self._Tmax,
            min_daily_temp=self._Tmin,
        )

        self._cell_values = self._grid["cell"]

        # Internal closed elevations field enables the component to disregard closed nodes when calculating PET
        self._gridCopy = copy.deepcopy(self._grid)
        self._gridCopy.add_field(
            "pet_status_at_node", self._grid.status_at_node, at="node"
        )
        self._cellular_status = map_node_to_cell(self._gridCopy, "pet_status_at_node")
        self._closed_elevations = (
            self._cellular_status == self._gridCopy.BC_NODE_IS_CLOSED
        )

        self._UpdateRad()

    @property
    def radiation(self):
        """User-provided net radiation field (W/m^2)

        radiation float, optional user-provided net radiation field.
        """
        return self._user_radiation

    @radiation.setter
    def radiation(self, radiation):
        self._user_radiation = radiation

    @property
    def Tmin(self):
        """Minimum temperature of the day (deg C)

        Tmin: float, required for 'Priestley Taylor' method.
        """
        return self._Tmin

    @Tmin.setter
    def Tmin(self, Tmin):
        self._Tmin = Tmin

    @property
    def Tmax(self):
        """Maximum temperature of the day (deg C)

        Tmax: float, required for 'Priestley Taylor' method.
        """
        return self._Tmax

    @Tmax.setter
    def Tmax(self, Tmax):
        self._Tmax = Tmax

    @property
    def Tavg(self):
        """Average temperature of the day (deg C)
        Tavg: float, required for 'Priestley Taylor'
        methods.
        """
        return self._Tavg

    @Tavg.setter
    def Tavg(self, Tavg):
        self._Tavg = Tavg

    @property
    def grid(self):
        return self._grid

    @grid.setter
    def grid(self, grid):
        self._grid = grid
        self._gridCopy = grid
        self._etpRad._grid = self._gridCopy

    def _process_field(self, field, field_name):
        if isinstance(field, np.ndarray) and np.shape(field) == np.shape(
            self._grid.at_node["topographic__elevation"]
        ):
            self._gridCopy.add_field(field_name, field, at="node")
            return map_node_to_cell(self._gridCopy, field_name)

        return field

    def _validate_temperature_range(self, min_temp, max_temp):
        if np.any(min_temp is None) or np.any(max_temp is None):
            raise ValueError("Tmin and Tmax are required fields")
        if np.any(min_temp > max_temp):
            raise ValueError("Tmin must be less than Tmax")

        min_temp = self._process_field(min_temp, "min_temperature")
        max_temp = self._process_field(max_temp, "max_temperature")

        return min_temp, max_temp

    def _validate_lai(self, lai):
        return self._process_field(lai, "leaf_area_index")

    # Function used to adjust field/variable values that would raise errors
    def _fix_values(self, field, error_value, fixed_value):
        if isinstance(field, np.ndarray):
            if np.any(field == error_value):
                field[field == error_value] = fixed_value

        elif field == error_value:
            field = fixed_value

    def update(self):
        """Update fields with current conditions.

        If the 'PenmanMonteith' method is used, this method looks to the properties of
        ``Tmin``, ``Tavg``, ``Zveg``, and all the other Penman / vegetation factors.

        If the 'PriestleyTaylor' method is used, this method looks to the
        values of the ``Tmin``, ``Tmax``, and ``Tavg`` properties.

        If the 'NetRadEqPE' method is used, this method looks to the values of
        the ``radiation`` (considered radiation field) and
        ``latent_heat_of_vaporization`` properties.
        """

        self._UpdateRad()

        if self._method == "PriestleyTaylor":
            self._PET_value = self._PriestleyTaylor()

        elif self._method == "PenmanMonteith":
            self._PET_value = self._PenmanMonteith()

        elif self._method == "NetRadEqPE":
            self._PET_value = self._NetRadEqPE()

        if isinstance(self._PET_value, np.ndarray):
            self._PET_value[self._closed_elevations] = 0.0

        self._cell_values["surface__potential_evapotranspiration_rate"][
            :
        ] = self._PET_value

    def _PriestleyTaylor(self):
        self._ETp = np.maximum(
            self._inp_alpha
            * (self._delta / (self._delta + self._y))
            * (self._input_rad / self._pwhv),
            0,
        )
        return self._ETp

    # Penman Monteith method
    def _PenmanMonteith(self):
        self._deltaTerm = self._delta * (self._input_rad - self._G)
        # LAIa is 50% of the LAI
        # rl is stomatal resistance of the cell wall
        self._LAIa = self._LAI * 0.5
        self._rs = self._rl / self._LAIa

        self._fix_values(self._vz, 0, 1.0)

        # Ensure that wind measurements are done above the vegetation canopy
        if self._is_zveg_defined and np.any((self._zm - self._zveg) <= 0):
            raise ValueError(
                f"""Elevation of wind speed observation, Zm ({self._zm}) must exceed
                the vegetation height, Zdveg ({self._zveg}).
                if wind was measured at a local station, use a logarithmic wind profile to estimate wind speed above canopy height."""
            )
        # If user didn't provide Zveg, ensure that wind measurements are done above the zero-plane displacement height
        elif np.any((self._zm - self._zd) <= 0):
            raise ValueError(
                f"""Elevation of wind speed observation, Zm ({self._zm}) must exceed
                the zero-plane displacement height, ({self._zd})."""
            )

        # Aerodynamic resistance
        self._ra = (np.log((self._zm - self._zd) / self._zo)) ** 2 / (
            self._k**2 * self._vz
        )

        self._vaporTerm = self._ca * self._pa * (self._es - self._ea) / self._ra

        self._denom = self._pwhv * (self._delta + self._y * (1 + self._rs / self._ra))
        self._ETp = (self._deltaTerm + self._vaporTerm) / self._denom

        return self._ETp

    # Net Rad equivalent PE method
    def _NetRadEqPE(self):
        self._E = self._input_rad / self._pwhv

        return self._E

    # Called before every PET update, ensures internal
    # calculations use user-updated variables
    def _UpdateRad(self):
        self._Tmin, self._Tmax = self._validate_temperature_range(
            self._Tmin, self._Tmax
        )

        self._LAI = self._validate_lai(self._LAI)

        # Update internal radiation component variables
        self._etpRad._current_time = self._current_time
        self._etpRad._Tmax = self._Tmax
        self._etpRad._Tmin = self._Tmin
        self._etpRad._a = self._a

        self._etpRad._latitude = self._latitude

        self._etpRad.update()

        self._input_rad = (
            self._etpRad._cell_values["radiation__net_flux"]
            if np.all(self._user_radiation is None)
            else self._user_radiation
        )

        # Ground heat flux, 10% of net rad
        self._G = self._input_rad * 0.1

        self._inp_alpha = self._alpha

        if not isinstance(self._Tavg, np.ndarray) and self._Tavg == 17.5:
            self._Tavg = (self._Tmax + self._Tmin) / 2

        self._inp_temp = self._Tavg
        if self._temperatures is not None:
            self._inp_temp = self._temperatures

        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)
        self._es = 0.6108 * np.exp((17.27 * self._inp_temp) / (237.7 + self._inp_temp))

        # Actual Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 8, (38)
        self._ea = (
            0.6108
            * np.exp((17.27 * self._Tmin) / (237.7 + self._Tmin))
            * self._relative_humidity
        )

        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es) / ((237.3 + self._inp_temp) ** 2.0)

        # Checks if user-provided vegetation variables are spatially distributed
        is_zveg_field = isinstance(self._zveg, np.ndarray)
        is_zo_field = isinstance(self._zo, np.ndarray)
        is_zd_field = isinstance(self._zd, np.ndarray)

        # Used in Penman method
        self._is_zveg_defined = (self._zveg is not None) if not is_zveg_field else True
        is_zo_default = (self._zo == 0.012) if not is_zo_field else False
        is_zd_default = (self._zd == 0.084) if not is_zd_field else False

        cell_sz = self._grid["cell"].size

        # Assert proper dimensions for these fields, all sizes must match grid cell dimensions
        if is_zveg_field and self._zveg.size != cell_sz:
            raise ValueError(
                f"""Zveg field dimensions must match grid dimensions, Zveg size was
                {self._zveg.size} while grid size was {cell_sz}"""
            )

        if is_zo_field and self._zo.size != cell_sz:
            raise ValueError(
                f"""Zo field dimensions must match grid dimensions, Zo size was
                {self._zo.size} while grid size was {cell_sz}"""
            )

        if is_zd_field and self._zd.size != cell_sz:
            raise ValueError(
                f"""Zd field dimensions must match grid dimensions, Zd size was
                {self._zd.size} while grid size was {cell_sz}"""
            )

        if self._is_zveg_defined and is_zo_default and is_zd_default:
            self._zd = self._zveg * 0.7
            self._zo = self._zveg * 0.1

        # von karman constant for air and clear water,
        # should be smaller for sediment-laden flows
        self._k = 0.4

        if self._pa is None:
            self._pa = 3.47 * self._etpRad._P / (273.3 + self._inp_temp)
            self._pa *= 1000

        self._ca = 1.22
