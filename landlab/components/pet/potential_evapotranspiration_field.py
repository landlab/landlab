import numpy as np
from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import math

_VALID_METHODS = set(
    ["Constant", "Cosine", "MeasuredRadiationPT", "PenmanMonteith", "PriestleyTaylor"]
)


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class PotentialEvapotranspiration(Component):
    """
    Potential Evapotranspiration Component calculates spatially distributed
    potential evapotranspiration based on input radiation factor (spatial
    distribution of incoming radiation) using one of the five ways
    to calculate potential evapotranspiration. The methods
    are 'Constant','Cosine', 'MeasuredRadiationPT',
    'PenmanMonteith', and 'PriestleyTaylor'.

    Note: Calling 'PriestleyTaylor' or 'PenmanMonteith' methods
    will overwrite the following fields:
        "radiation__incoming_shortwave_flux",
        "radiation__net_shortwave_flux",
        "radiation__net_longwave_flux", and
        "radiation__net_flux".

    As mentioned above, this component can be configured to use
    one of the five ways to calculate potential evapotranspiration.
    The methods are 'Constant','Cosine', 'MeasuredRadiationPT',
    'PenmanMonteith', and 'PriestleyTaylor'.

    'Constant' method: Receives the user-input const_potential_evapotranspiration
    during a call to the 'update' method. This value is
    then distributed spatially using the 'radiation__ratio_to_flat_surface'
    field.

    'Cosine' method: This method calculates the flat surface potential
    evapotranspiration value as a sinusoidal function of the day of the year.
    This method uses the inputs of MeanTmaxF, delta_d, lt, and nd
    obtained during instantiation of the component or during a call
    to the 'initialize' method, and current_time during a call to the
    'update' method. This value is then distributed spatially using the
    'radiation__ratio_to_flat_surface' field. [Equation 28 - Zhou et al. 2013]

    'MeasuredRadiationPET' method: This method calculates the flat surface
    potential evapotranspiration value using the inputs of
    obs_radiation, albedo, and Tavg obtained during the call to
    'update' method. This value is then distributed spatially
    using the 'radiation__ratio_to_flat_surface' field.
    [ASCE-EWRI Task Committee Report Jan 2005]

    'PenmanMonteith' method: This method calculates the flat surface
    potential evapotranspiration values using Penman equation.
    This method uses the inputs of zm, zh, zveg, rl,
    air_density, latent_heat_of_vaporization,
    stefan_boltzmann_const, and albedo obtained during
    instantiation of the component or during a call
    to the 'initialize' method, and Tavg, radiation_sw,
    wind_speed, relative_humidity, and co2_concentration
    during a call to the 'update' method. This value is then
    distributed spatially using the
    'radiation__ratio_to_flat_surface' field.[Penman 1948]

    'PriestleyTaylor' method: This method calculates the flat surface
    potential evapotranspiration values using Priestley-Taylor
    equation. This method uses the inputs of
    stefan_boltzmann_const, latitude, elevation_of_measurement,
    adjustment_coeff, priestley_taylor_constant, albedo,
    psychometric_const, and latent_heat_of_vaporization obtained during
    instantiation of the component or during a call
    to the 'initialize' method, and current_time, Tmax,
    Tmin, and Tavg during a call to the 'update' method.
    This value is then distributed spatially using the
    'radiation__ratio_to_flat_surface' field.
    [Priestley and Taylor 1972; ASCE-EWRI Task Committee Report Jan 2005]


    Suggested References:

    Allen, R. G., Walter, I. A., Elliot, R., Howell, T., Itenfisu,
    D., Jensen, M., & Snyder, R. (2005). The ASCE standardized
    reference evapotranspiration equation. ASCE-EWRI task committee
    final report.

    Penman, H. L. (1948). Natural evaporation from open water, bare
    soil and grass. Proceedings of the Royal Society of London.
    Series A. Mathematical and Physical Sciences, 193(1032), 120-145.

    Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of
    surface heat flux and evaporation using large-scale parameters.
    Monthly weather review, 100(2), 81-92.

    Yang, Y., Roderick, M. L., Zhang, S., McVicar, T. R., & Donohue,
    R. J. (2019). Hydrologic implications of vegetation response to
    elevated CO 2 in climate projections. Nature Climate Change, 9(1), 44.

    Zhou, X., Istanbulluoglu, E., & Vivoni, E. R. (2013). Modeling the
    ecohydrological role of aspect controlled radiation on
    tree grass shrub coexistence in a semiarid climate.
    Water Resources Research, 49(5), 2872-2895.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import PotentialEvapotranspiration

    Let's create a raster grid with 5 rows and 4 columns
    of nodes (therefore, 6 cells).

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))

    To run the processes represented by PotentialEvapotranspiration
    class, we need to instantiate the class.

    >>> pet_obj = PotentialEvapotranspiration(grid)

    Note that the default method for PotentialEvapotranspiration
    is 'Cosine'.

    >>> pet_obj.name
    'Potential Evapotranspiration'
    >>> pet_obj.input_var_names
    ('radiation__ratio_to_flat_surface',)
    >>> sorted(pet_obj.output_var_names)
    ['radiation__incoming_shortwave_flux',
     'radiation__net_flux',
     'radiation__net_longwave_flux',
     'radiation__net_shortwave_flux',
     'surface__potential_evapotranspiration_rate']
    >>> sorted(pet_obj.units) # doctest: +NORMALIZE_WHITESPACE
    [('radiation__incoming_shortwave_flux', 'W/m^2'),
     ('radiation__net_flux', 'W/m^2'),
     ('radiation__net_longwave_flux', 'W/m^2'),
     ('radiation__net_shortwave_flux', 'W/m^2'),
     ('radiation__ratio_to_flat_surface', 'None'),
     ('surface__potential_evapotranspiration_rate', 'mm/d')]
    >>> pet_obj.grid.number_of_cell_rows
    3
    >>> pet_obj.grid.number_of_cell_columns
    2
    >>> pet_obj.grid is grid
    True
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, 0.)
    True
    >>> grid.at_cell['radiation__ratio_to_flat_surface'] = np.array([
    ...       0.38488566, 0.38488566,
    ...       0.33309785, 0.33309785,
    ...       0.37381705, 0.37381705])
    >>> current_time = 0.5
    >>> pet_obj.update(current_time)
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, np.array([
    ...       5.58, 5.58, 4.83,
    ...       4.83, 5.42, 5.42]),
    ...       rtol=1e-02)
    True

    A reminder here that potential evapotranspiration is calculated
    on a flat surface and is distributed over the sloped surface
    using a grid.at_cell['radiation__ratio_to_flat_surface'] that
    is passed as an input.

    If daily temperature data is available, we can
    use Priestly Taylor's method to calculate the
    potential evapotranspiration. Let us look at an
    example where:
    minimum daily temperature = 13.36 C,
    maximum daily temperature = 25.6 C,
    Day = 9th June (current_time = 0.438), and
    Location = Ainsworth, NE (latitude = 42.545 N),
    elevation = 786 m,
    albedo = 0.05, and
    priestly_taylor_constant = 1.26.
    let's re-initialize the 'surface__potential_evapotranspiration_rate'
    to zeros.

    >>> grid.at_cell['surface__potential_evapotranspiration_rate'] = np.zeros(
    ...     grid.number_of_cells)
    >>> pet_pt = PotentialEvapotranspiration(
    ...     grid, method='PriestleyTaylor',
    ...     albedo=0.05, latitude=42.545,
    ...     elevation_of_measurement=786,
    ...     priestley_taylor_const=1.26)
    >>> pet_pt.update(
    ...     current_time=0.438, Tmax=25.6,
    ...     Tmin=13.36, Tavg=19.48)
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, np.array([
    ...     3.99, 3.99, 3.46,
    ...     3.46, 3.88, 3.88]),
    ...     rtol=2)
    True

    If we have observed radiation at this location,
    we can use 'MeasuredRadiationPT' method to calculate
    Potential Evapotranspiration using Priestly Taylor equation.
    Let us assume, observed radiation = 361.17 W/m^2.

    >>> pet_mpt = PotentialEvapotranspiration(
    ...     grid, method='MeasuredRadiationPT',
    ...     albedo=0.05, latitude=42.545,
    ...     elevation_of_measurement=786,
    ...     priestley_taylor_const=1.26)
    >>> pet_mpt.update(
    ...     Tavg=19.48, obs_radiation=361.17)
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, np.array([
    ...     3.99, 3.99, 3.46,
    ...     3.46, 3.88, 3.88]),
    ...     rtol=2)
    True

    Now, let's look at an example using 'PenmanMonteith' method
    to calculate potential evapotranspiration. To reuse the
    grid, let's re-initialize the 'surface__potential_evapotranspiration_rate'
    to zeros.

    >>> grid.at_cell['surface__potential_evapotranspiration_rate'] = np.zeros(
    ...    grid.number_of_cells)

    For this example, let us consider the following:
    albedo = 0.05,
    average air temperature = 19.48 C,
    observed shortwave radiation = 307 W/m^2,
    relative humidity = 33% (measured at a height of 2 m),
    wind speed = 3.12 m/s (measured at a height of 2 m),
    height of vegetation = 0.12 m,
    leaf area index of reference crop = 2.88, and
    stomatal resistance of well illuminated leaf = 70 m/s.

    >>> pet_pm = PotentialEvapotranspiration(
    ...     grid, method='PenmanMonteith',
    ...     zm=2, zveg=0.12, lai=2.88, zh=2,
    ...     rl=70, albedo=0.05)
    >>> pet_pm.update(Tavg=19.48, obs_radiation=307,
    ...               wind_speed=3.12, relative_humidity=33)
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, np.array([
    ...     3.56, 3.56, 3.08,
    ...     3.08, 3.46, 3.46]),
    ...     rtol=2)
    True

    If potential evapotranspiration value on a flat
    surface is already available, we can distribute
    it on a sloped surface using the 'Constant' method.

    >>> pet_c = PotentialEvapotranspiration(
    ... grid, method='Constant')
    >>> pet_c.update(const_potential_evapotranspiration=5)
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, np.array([
    ...     1.92, 1.92, 1.67,
    ...     1.67, 1.87, 1.87]),
    ...     rtol=2)
    True
    """

    _name = "Potential Evapotranspiration"

    _input_var_names = ("radiation__ratio_to_flat_surface",)

    _output_var_names = (
        "surface__potential_evapotranspiration_rate",
        "radiation__incoming_shortwave_flux",
        "radiation__net_shortwave_flux",
        "radiation__net_longwave_flux",
        "radiation__net_flux",
    )

    _var_units = {
        "radiation__ratio_to_flat_surface": "None",
        "surface__potential_evapotranspiration_rate": "mm/d",
        "radiation__incoming_shortwave_flux": "W/m^2",
        "radiation__net_shortwave_flux": "W/m^2",
        "radiation__net_longwave_flux": "W/m^2",
        "radiation__net_flux": "W/m^2",
    }

    _var_mapping = {
        "radiation__ratio_to_flat_surface": "cell",
        "surface__potential_evapotranspiration_rate": "cell",
        "radiation__incoming_shortwave_flux": "cell",
        "radiation__net_shortwave_flux": "cell",
        "radiation__net_longwave_flux": "cell",
        "radiation__net_flux": "cell",
    }

    _var_doc = {
        "radiation__ratio_to_flat_surface": "ratio of total incident shortwave radiation on sloped surface \
             to flat surface",
        "surface__potential_evapotranspiration_rate": "potential sum of evaporation and potential transpiration",
        "radiation__incoming_shortwave_flux": "total incident shortwave radiation over the time step",
        "radiation__net_shortwave_flux": "net incident shortwave radiation over the time step",
        "radiation__net_longwave_flux": "net incident longwave radiation over the time step",
        "radiation__net_flux": "net total radiation over the time step",
    }

    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        method="Cosine",
        priestley_taylor_const=1.26,
        albedo=0.12,
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
        rl=130,
        zveg=0.3,
        lai=2.0,
        zm=3.3,
        zh=2,
        air_density=1.22,
        **kwds,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        method: {'Constant', 'PriestleyTaylor', 'MeasuredRadiationPT',
                 'Cosine', 'PenmanMonteith'}, optional
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
        delta_d: float, required only for method(s): Cosine
            Calibrated difference between max & min daily TmaxF (mm/d).
        zm: float, required only for method(s): PenmanMonteith
            Wind speed anemometer height (m).
        zh: float, required only for method(s): PenmanMonteith (with
            correction for Relative Humidity measurement elevation)
            Relative Humidity measurement height (m).
        rl: float, required only for method(s): PenmanMonteith
            Stomatal resistance of a well-illuminated leaf (s/m)
        lai: float, required only for method(s): PenmanMonteith
            Leaf Area Index of reference crop (-).
        air_density: float, optional
            Density of air (Kg/m^3).
        """

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
        self._zm = zm  # (m) wind speed anemometer height
        self._zh = zh  # (m) Relative Humidity probe height
        self._lai = lai  # Leaf Area Index
        self._zveg = zveg  # (m) Vegetation height
        self._rl = rl  # (sec/m) reverse of conductance
        self._rho_w = 1000.0  # (Kg/m^3) Density of water
        self._rho_a = air_density  # (Kg/m^3) Density of Air
        self._ca = 1000.0  # (J/Kg/C) or (Ws/Kg/C) Specific heat of air
        self._von_karman = 0.41  # Von Karman Constant
        _assert_method_is_valid(self._method)

        super(PotentialEvapotranspiration, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros("cell", name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros("cell", name, units=self._var_units[name])

        self._cell_values = self.grid["cell"]

    def update(
        self,
        current_time=None,
        const_potential_evapotranspiration=12.0,
        Tmin=None,
        Tmax=None,
        Tavg=None,
        obs_radiation=None,
        relative_humidity=None,
        wind_speed=None,
        co2_concentration=300.0,
        ground_heat_flux=0.0,
        **kwds,
    ):
        """Update fields with current conditions.

        Parameters
        ----------
        current_time: float, required only for method(s): Cosine
            Current time (Years)
        constant_potential_evapotranspiration: float, required only for
        method(s): Constant
            Constant PET value to be spatially distributed.
        Tmin: float, required for method(s): Priestley Taylor
            Minimum temperature of the day (deg C)
        Tmax: float, required for method(s): Priestley Taylor
            Maximum temperature of the day (deg C)
        Tavg: float, required for method(s): Priestley Taylor,
        MeasuredRadiationPT, and PenmanMonteith
            Average temperature of the day (deg C)
        obs_radiation: float, required for method(s): MeasuredRadiationPT, and
        PenmanMonteith
            Observed radiation (W/m^2)
        relative_humidity: float, required for method(s): PenmanMonteith
            Observed relative humidity (%)
        wind_speed: float, required for method(s): PenmanMonteith
            Observed wind speed (m/s)
        co2_concentration: float (default=300.),
            CO2 concentration (ppm)
        ground_heat_flux: float (default=0)
            Ground heat flux (W/m^2)
        """
        if self._method in ["PriestleyTaylor", "MeasuredRadiationPT", "PenmanMonteith"]:
            if Tavg is None:
                Tavg = (Tmax + Tmin) / 2.0
        if self._method == "Constant":
            self._PET_value = const_potential_evapotranspiration
        elif self._method == "PriestleyTaylor":
            self._PET_value = self._PriestleyTaylor(current_time, Tmax, Tmin, Tavg)
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
            Robs = obs_radiation
            self._PET_value = self._MeasuredRadPT(Tavg, (1 - self._a) * Robs)
        elif self._method == "Cosine":
            self._J = np.floor((current_time - np.floor(current_time)) * 365.0)
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
        elif self._method == "PenmanMonteith":
            self._PET_value = self._PenmanMonteith(
                Tavg,
                obs_radiation,
                wind_speed,
                relative_humidity,
                co2_concentration,
                ground_heat_flux,
            )
            if math.isnan(self._PET_value):
                self._PET_value = 0.0
            if self._PET_value < 0.0:
                self._PET_value = 0.0
            self._cell_values["radiation__incoming_shortwave_flux"] = (
                obs_radiation * self._cell_values["radiation__ratio_to_flat_surface"]
            )
            self._cell_values["radiation__net_shortwave_flux"] = (
                obs_radiation
                * (1 - self._a)
                * self._cell_values["radiation__ratio_to_flat_surface"]
            )
            self._cell_values["radiation__net_longwave_flux"] = (
                self._radiation_lw
                * self._cell_values["radiation__ratio_to_flat_surface"]
            )
            self._cell_values["radiation__net_flux"] = (
                self._net_radiation
                * self._cell_values["radiation__ratio_to_flat_surface"]
            )

        # Spatially distributing PET
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
            (-1 * np.tan(self._phi) * np.tan(self._sdecl)) / (self._x ** 2.0)
        )

        # Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 21, (48)
        # 11.57 converts 1 MJ/m^2/day to W/m^2
        self._Ra = (11.57 * (24.0 / np.pi) * 4.92 * self._dr) * (
            (self._ws * np.sin(self._phi) * np.sin(self._sdecl))
            + (np.cos(self._phi) * np.cos(self._sdecl) * (np.sin(self._ws)))
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

    def _PenmanMonteith(
        self,
        Tavg,
        radiation_sw,
        wind_speed,
        relative_humidity,
        co2_concentration,
        ground_heat_flux,
    ):
        zm = self._zm
        zh = self._zh
        zd = 0.7 * self._zveg  # (m) zero-plane displacement height
        z0 = 0.123 * self._zveg
        z0h = 0.1 * z0
        kv = self._von_karman

        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)
        self._es = 0.6108 * np.exp((17.27 * Tavg) / (237.7 + Tavg))
        # Actual Vapor Pressure
        self._ea = self._es * relative_humidity * 0.01
        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es) / ((237.3 + Tavg) ** 2.0)
        self._ra = ((np.log((zm - zd) / z0)) ** 2) / ((kv ** 2) * wind_speed)
        if self._zh is not None:
            self._ra = ((np.log((zm - zd) / z0)) * (np.log((zh - zd) / z0h))) / (
                kv ** 2 * wind_speed
            )  # Correction for RH

        self._Kat = 2.158 / (self._rho_w * (273.3 + Tavg))
        # Evaporation
        self._E = (
            self._Kat * (1.0 / self._ra) * (self._es - self._ea) * 86400.0 * 1000.0
        )
        # Potential Evapotranspiration
        self._E2 = radiation_sw / self._pwhv
        self._Ts = Tavg + 273.15 - 0.825 * np.exp(3.54 * (10.0 ** (-3)) * radiation_sw)
        self._radiation_lw = max(
            self._sigma * (self._Ts ** 4 - (Tavg + 273.15) ** 4), 0.0
        )
        self._net_radiation = max(
            ((1 - self._a) * radiation_sw + self._radiation_lw - ground_heat_flux), 0.0,
        )
        self._penman_numerator = (self._delta * self._net_radiation) + (
            self._rho_a * self._ca * (self._es - self._ea) / self._ra
        )
        self._rs = self._rl / (0.5 * self._lai) + 0.05 * (
            co2_concentration - 300.0
        )  # Yang et al. 2019
        self._penman_denominator = self._pwhv * (
            self._delta + self._y * (1 + self._rs / self._ra)
        )
        self._ETp = self._penman_numerator / self._penman_denominator
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
