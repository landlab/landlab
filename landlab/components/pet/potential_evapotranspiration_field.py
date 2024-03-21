import numpy as np

from landlab import Component

_VALID_METHODS = {"Constant", "PriestleyTaylor", "MeasuredRadiationPT", "Cosine"}


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class PotentialEvapotranspiration(Component):
    """
    Potential Evapotranspiration Component calculates spatially distributed
    potential evapotranspiration based on input radiation factor (spatial
    distribution of incoming radiation) using chosen method such as constant
    or Priestley Taylor. Ref: Xiaochi et. al. 2013 for 'Cosine' method and
    ASCE-EWRI Task Committee Report Jan 2005 for 'PriestleyTaylor' method.
    Note: Calling 'PriestleyTaylor' method would generate/overwrite shortwave &
    longwave radiation fields.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.pet import PotentialEvapotranspiration

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    >>> grid["cell"]["radiation__ratio_to_flat_surface"] = np.array(
    ...     [0.38488566, 0.38488566, 0.33309785, 0.33309785, 0.37381705, 0.37381705]
    ... )
    >>> PET = PotentialEvapotranspiration(grid)
    >>> PET.name
    'PotentialEvapotranspiration'
    >>> PET.input_var_names
    ('radiation__ratio_to_flat_surface',)
    >>> sorted(PET.output_var_names)
    ['radiation__incoming_shortwave_flux',
     'radiation__net_flux',
     'radiation__net_longwave_flux',
     'radiation__net_shortwave_flux',
     'surface__potential_evapotranspiration_rate']
    >>> sorted(PET.units)
    [('radiation__incoming_shortwave_flux', 'W/m^2'),
     ('radiation__net_flux', 'W/m^2'),
     ('radiation__net_longwave_flux', 'W/m^2'),
     ('radiation__net_shortwave_flux', 'W/m^2'),
     ('radiation__ratio_to_flat_surface', 'None'),
     ('surface__potential_evapotranspiration_rate', 'mm')]
    >>> PET.grid.number_of_cell_rows
    3
    >>> PET.grid.number_of_cell_columns
    2
    >>> PET.grid is grid
    True
    >>> pet_rate = grid.at_cell["surface__potential_evapotranspiration_rate"]
    >>> np.allclose(pet_rate, 0.0)
    True
    >>> PET.current_time = 0.5
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
        "radiation__incoming_shortwave_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "incident shortwave radiation",
        },
        "radiation__net_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "net radiation",
        },
        "radiation__net_longwave_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "net incident longwave radiation",
        },
        "radiation__net_shortwave_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "net incident shortwave radiation",
        },
        "radiation__ratio_to_flat_surface": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": (
                "ratio of incident shortwave radiation on sloped "
                "surface to flat surface"
            ),
        },
        "surface__potential_evapotranspiration_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm",
            "mapping": "cell",
            "doc": "potential sum of evaporation and potential transpiration",
        },
    }

    def __init__(
        self,
        grid,
        method="Cosine",
        priestley_taylor_const=1.26,
        albedo=0.6,
        latent_heat_of_vaporization=28.34,
        psychometric_const=0.066,
        stefan_boltzmann_const=0.0000000567,
        solar_const=1366.67,
        latitude=34.0,
        elevation_of_measurement=300,
        adjustment_coeff=0.18,
        lt=0.0,
        nd=365.0,
        MeanTmaxF=12.0,
        delta_d=5.0,
        current_time=None,
        const_potential_evapotranspiration=12.0,
        Tmin=0.0,
        Tmax=1.0,
        Tavg=0.5,
        obs_radiation=350.0,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        method: {'Constant', 'PriestleyTaylor', 'MeasuredRadiationPT', 'Cosine'}, optional
            Priestley Taylor method will spit out radiation outputs too.
        priestley_taylor_constant: float, optional
            Alpha used in Priestley Taylor method.
        albedo: float, optional
            Albedo.
        latent_heat_of_vaporization: float, optional
            Latent heat of vaporization for water Pwhv (Wd/(m*mm^2)).
        psychometric_const: float, optional
            Psychometric constant (kPa (deg C)^-1).
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
        const_potential_evapotranspiration: float, optional for
            'Constant' method
            Constant PET value to be spatially distributed.
        Tmin: float, required for 'Priestley Taylor' method
            Minimum temperature of the day (deg C)
        Tmax: float, required for 'Priestley Taylor' method
            Maximum temperature of the day (deg C)
        Tavg: float, required for 'Priestley Taylor' and 'MeasuredRadiationPT'
            methods
            Average temperature of the day (deg C)
        obs_radiation float, required for 'MeasuredRadiationPT' method
            Observed radiation (W/m^2)
        """
        super().__init__(grid)

        self.current_time = current_time
        self.const_potential_evapotranspiration = const_potential_evapotranspiration
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Tavg = Tavg
        self.obs_radiation = obs_radiation

        self._method = method
        # For Priestley Taylor
        self._alpha = priestley_taylor_const
        self._a = albedo
        self._pwhv = latent_heat_of_vaporization
        self._y = psychometric_const
        self._sigma = stefan_boltzmann_const
        self._Gsc = solar_const
        self._phi = (np.pi / 180.0) * latitude
        self._z = elevation_of_measurement
        self._Krs = adjustment_coeff
        self._LT = lt
        self._ND = nd
        self._TmaxF_mean = MeanTmaxF
        self._DeltaD = delta_d
        _assert_method_is_valid(self._method)

        self.initialize_output_fields()

        self._cell_values = self._grid["cell"]

    @property
    def const_potential_evapotranspiration(self):
        """Constant PET value to be spatially distributed.

        Used by 'Constant' method.
        """
        return self._const_potential_evapotranspiration

    @const_potential_evapotranspiration.setter
    def const_potential_evapotranspiration(self, const_potential_evapotranspiration):
        self._const_potential_evapotranspiration = const_potential_evapotranspiration

    @property
    def obs_radiation(self):
        """Observed radiation (W/m^2)

        obs_radiation float, required for 'MeasuredRadiationPT' method.
        """
        return self._obs_radiation

    @obs_radiation.setter
    def obs_radiation(self, obs_radiation):
        self._obs_radiation = obs_radiation

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

        Tavg: float, required for 'Priestley Taylor' and 'MeasuredRadiationPT'
        methods.
        """
        return self._Tavg

    @Tavg.setter
    def Tavg(self, Tavg):
        self._Tavg = Tavg

    def update(self):
        """Update fields with current conditions.

        If the 'Constant' method is used, this method looks to the value of
        the ``const_potential_evapotranspiration`` property.

        If the 'PriestleyTaylor' method is used, this method looks to the
        values of the ``Tmin``, ``Tmax``, and ``Tavg`` properties.

        If the 'MeasuredRadiationPT' method is use this method looks to the
        values of the ``Tavg`` and ``obs_radiation`` property.
        """

        if self._method == "Constant":
            self._PET_value = self._const_potential_evapotranspiration
        elif self._method == "PriestleyTaylor":
            self._PET_value = self._PriestleyTaylor(
                self._current_time, self._Tmax, self._Tmin, self._Tavg
            )
            self._cell_values["radiation__incoming_shortwave_flux"] = (
                self._Rs * self._cell_values["radiation__ratio_to_flat_surface"]
            )
            self._cell_values["radiation__net_shortwave_flux"] = (
                self._Rns * self._cell_values["radiation__ratio_to_flat_surface"]
            )
            self._cell_values["radiation__net_longwave_flux"] = (
                self._Rnl * self._cell_values["radiation__ratio_to_flat_surface"]
            )
            self._cell_values["radiation__net_flux"] = (
                self._Rn * self._cell_values["radiation__ratio_to_flat_surface"]
            )
        elif self._method == "MeasuredRadiationPT":
            Robs = self._obs_radiation
            self._PET_value = self._MeasuredRadPT(self._Tavg, (1 - self._a) * Robs)
        elif self._method == "Cosine":
            self._J = np.floor(
                (self._current_time - np.floor(self._current_time)) * 365.0
            )
            self._PET_value = max(
                (
                    self._TmaxF_mean
                    + self._DeltaD
                    / 2.0
                    * np.cos(
                        (2 * np.pi) * (self._J - self._LT - self._ND / 2) / self._ND
                    )
                ),
                0.0,
            )

        self._PET = (
            self._PET_value * self._cell_values["radiation__ratio_to_flat_surface"]
        )
        self._cell_values["surface__potential_evapotranspiration_rate"][:] = self._PET

    def _PriestleyTaylor(self, current_time, Tmax, Tmin, Tavg):
        # Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn 25, (52)
        self._J = np.floor((current_time - np.floor(current_time)) * 365)
        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)
        self._es = 0.6108 * np.exp((17.27 * Tavg) / (237.7 + Tavg))

        # Actual Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 8, (38)
        self._ea = 0.6108 * np.exp((17.27 * Tmin) / (237.7 + Tmin))

        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es) / ((237.3 + Tavg) ** 2.0)

        # Solar Declination Angle - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 24,(51)
        self._sdecl = 0.409 * np.sin(((np.pi / 180.0) * self._J) - 1.39)

        # Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 23,(50)
        self._dr = 1 + (0.033 * np.cos(np.pi / 180.0 * self._J))

        # To calculate ws - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 29,(61)
        self._x = 1.0 - (((np.tan(self._phi)) ** 2.0) * (np.tan(self._sdecl) ** 2.0))
        if self._x <= 0:
            self._x = 0.00001
            # Sunset Hour Angle - ASCE-EWRI Task Committee Report,
            # Jan-2005 - Eqn 28,(60)
        self._ws = (np.pi / 2.0) - np.arctan(
            (-1 * np.tan(self._phi) * np.tan(self._sdecl)) / (self._x**2.0)
        )

        # Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 21, (48)
        # 11.57 converts 1 MJ/m^2/day to W/m^2
        self._Ra = (
            11.57
            * (24.0 / np.pi)
            * 4.92
            * self._dr
            * (
                (self._ws * np.sin(self._phi) * np.sin(self._sdecl))
                + (np.cos(self._phi) * np.cos(self._sdecl) * (np.sin(self._ws)))
            )
        )

        # Clear-sky Solar Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 19, (47)
        self._Rso = (0.75 + ((2.0 * (10 ** (-5.0))) * self._z)) * self._Ra
        self._Rs = min(self._Krs * self._Ra * np.sqrt(Tmax - Tmin), self._Rso)

        # Net Short Wave Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 16, (43)
        self._Rns = self._Rs * (1 - self._a)

        # Relative Cloudiness - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Page 20,35
        if self._Rso > 0:
            self._u = self._Rs / self._Rso
        else:
            self._u = 0

        if self._u < 0.3:
            self._u = 0.3
        elif self._u > 1:
            self._u = 1.0

        # Cloudiness Function - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 18, (45)
        self._fcd = (1.35 * self._u) - 0.35

        # Net Long Wave Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 17, (44)
        self._Rnl = (
            self._sigma
            * self._fcd
            * (
                0.34
                - (0.14 * np.sqrt(self._ea))
                * (((Tmax + 273.16) ** 4.0 + (Tmin + 273.16) ** 4.0) / 2.0)
            )
        )

        # Net Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 15, (42)
        self._Rn = self._Rns - self._Rnl

        self._ETp = max(
            self._alpha
            * (self._delta / (self._delta + self._y))
            * (self._Rn / self._pwhv),
            0,
        )

        return self._ETp

    def _MeasuredRadPT(self, Tavg, Rnobs):
        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)
        self._es = 0.6108 * np.exp((17.27 * Tavg) / (237.7 + Tavg))

        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es) / ((237.3 + Tavg) ** 2.0)
        self._ETp = max(
            self._alpha
            * (self._delta / (self._delta + self._y))
            * (Rnobs / self._pwhv),
            0,
        )
        return self._ETp
