import numpy as np

from landlab import Component

_VALID_METHODS = {"Grid"}


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class Radiation(Component):

    """Compute 1D and 2D total incident shortwave radiation.

    Landlab component that computes 1D and 2D total incident shortwave
    radiation. This code also computes relative incidence shortwave radiation
    compared to a flat surface. Ref: Bras, Rafael L. Hydrology: an
    introduction to hydrologic science. Addison Wesley Publishing Company,
    1990.

    .. codeauthor:: Sai Nudurupati & Erkan Istanbulluoglu & Berkan Mertan

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Radiation
    >>> import numpy as np
    >>> from numpy.testing import assert_array_almost_equal

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    >>> z = grid.add_zeros("node", "topographic__elevation")
    >>> rad = Radiation(grid)
    >>> rad.name
    'Radiation'
    >>> rad.input_var_names
    ('topographic__elevation',)
    >>> sorted(rad.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['radiation__incoming_shortwave_flux',
     'radiation__net_flux',
     'radiation__net_longwave_flux',
     'radiation__net_shortwave_flux',
     'radiation__ratio_to_flat_surface']
    >>> sorted(rad.units) # doctest: +NORMALIZE_WHITESPACE
    [('radiation__incoming_shortwave_flux', 'W/m^2'),
     ('radiation__net_flux', 'W/m^2'),
     ('radiation__net_longwave_flux', 'W/m^2'),
     ('radiation__net_shortwave_flux', 'W/m^2'),
     ('radiation__ratio_to_flat_surface', 'None'),
     ('topographic__elevation', 'm')]

    >>> rad.grid.number_of_node_rows
    5
    >>> rad.grid.number_of_node_columns
    4
    >>> rad.grid is grid
    True
    >>> np.all(grid.at_cell['radiation__ratio_to_flat_surface'] == 0.)
    True
    >>> np.all(grid.at_cell['radiation__incoming_shortwave_flux'] == 0.)
    True
    >>> np.all(grid.at_cell['radiation__net_shortwave_flux'] == 0.)
    True
    >>> np.all(grid.at_cell['radiation__net_longwave_flux'] == 0.)
    True
    >>> np.all(grid.at_cell['radiation__net_flux'] == 0.)
    True
    >>> np.all(grid.at_node['topographic__elevation'] == 0.)
    True

    >>> grid['node']['topographic__elevation'] = np.array([
    ...       0., 0., 0., 0.,
    ...       1., 1., 1., 1.,
    ...       2., 2., 2., 2.,
    ...       3., 4., 4., 3.,
    ...       4., 4., 4., 4.])
    >>> rad.current_time = 0.5
    >>> rad.update()
    >>> np.all(grid.at_cell['radiation__ratio_to_flat_surface'] == 0.)
    False
    >>> np.all(grid.at_cell['radiation__incoming_shortwave_flux'] == 0.)
    False
    >>> np.all(grid.at_cell['radiation__net_shortwave_flux'] == 0.)
    False
    >>> np.all(grid.at_cell['radiation__net_longwave_flux'] == 0.)
    False
    >>> np.all(grid.at_cell['radiation__net_flux'] == 0.)
    False
    >>> grid['node']['topographic__elevation'] = np.array([
    ...       0., 0., 0., 0.,
    ...       100., 100., 100., 100.,
    ...       200., 200., 200., 200.,
    ...       300., 400., 400., 300.,
    ...       400., 400., 400., 400.])
    >>> calc_rad = Radiation(grid, current_time=0.0)
    >>> calc_rad.update()
    >>> proven_net_shortwave_field = [189.608, 189.608, 188.839,
    ...                                 188.686, 184.152, 183.741]
    >>> nsflux = grid.at_cell['radiation__net_shortwave_flux']
    >>> print(assert_array_almost_equal(proven_net_shortwave_field, nsflux, decimal=3))
    None

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Bras, R. L.: Hydrology: an introduction to hydrologic science, Addison
    Wesley Publishing Company, Boston, Mass., USA, 643 pp., 1990.

    """

    _name = "Radiation"

    _unit_agnostic = False

    # Setting a really arbitrary value to check
    # if variable was custom initialized during instantiation
    _default_assertion_val = -999.999

    _info = {
        "radiation__incoming_shortwave_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "total incident shortwave radiation over the time step",
        },
        "radiation__net_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "net total radiation over the time step",
        },
        "radiation__net_longwave_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "net incident longwave radiation over the time step",
        },
        "radiation__net_shortwave_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "net incident shortwave radiation over the time step",
        },
        "radiation__ratio_to_flat_surface": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": (
                "ratio of total incident shortwave radiation on sloped "
                "surface to flat surface"
            ),
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        method="Grid",
        cloudiness=0.2,
        latitude=34.0,
        albedo=0.2,
        solar_constant=1366.67,
        clearsky_turbidity=_default_assertion_val,
        opt_airmass=_default_assertion_val,
        current_time=0.0,
        hour=12.0,
        average_daily_temp=17.5,
        max_daily_temp=25.0,
        min_daily_temp=10.0,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        method: {'Grid'}, optional
            Currently, only default is available.
        cloudiness: float, optional
            Cloudiness.
        latitude: float, optional
            Latitude (radians).
        albedo: float, optional
            Albedo.
        solar_constant: float, optional
            Solar Constant (W/m^2).
        clearsky_turbidity: float, optional
            Clear sky turbidity.
        opt_airmass: float, optional
            Optical air mass.
        average_daily_temp: float, optional
            Average daily temperature (Celsius)
        max_daily_temp: float, optional
            Maximum daily temperature (Celsius)
        min_daily_Temp: float, optional
            Minimum daily temperature (Celsius)
        current_time: float
              Current time (years).
        hour: float, optional
              Hour of the day. Default is 12 (solar noon)
        """
        super().__init__(grid)

        self.current_time = current_time
        self.hour = hour

        self._method = method
        self._N = cloudiness
        self._latitude = latitude
        self._A = albedo
        self._Io = solar_constant
        self._n = clearsky_turbidity
        self._m = opt_airmass

        # For computations requiring temperature
        self._Tavg = average_daily_temp
        self._Tmax = max_daily_temp
        self._Tmin = min_daily_temp

        _assert_method_is_valid(self._method)

        self.initialize_output_fields()

        if "Slope" not in self._grid.at_cell:
            self._grid.add_zeros("Slope", at="cell", units="radians")

        if "Aspect" not in self._grid.at_cell:
            self._grid.add_zeros("Aspect", at="cell", units="radians")

        self._nodal_values = self._grid["node"]
        self._cell_values = self._grid["cell"]

        # The elevation used in calculations is taken from
        # topographic elevation field
        self._max_elev = np.max(self._nodal_values["topographic__elevation"])
        self._avg_elev = np.average(self._nodal_values["topographic__elevation"])

        self._slope, self._aspect = grid.calculate_slope_aspect_at_nodes_burrough(
            vals="topographic__elevation"
        )

        self._cell_values["Slope"] = self._slope
        self._cell_values["Aspect"] = self._aspect

    @property
    def hour(self):
        """Hour of the day.

        Default is 12 (solar noon).
        """
        return self._hour

    @hour.setter
    def hour(self, hour):
        assert hour >= 0.0
        assert hour <= 24.0
        self._hour = hour

    def update(self):
        """Update fields with current loading conditions.

        This method looks to the properties ``current_time`` and
        ``hour`` and uses their values in updating fields.
        """

        self._t = self._hour
        self._radf = self._cell_values["radiation__ratio_to_flat_surface"]
        self._Rs = self._cell_values["radiation__incoming_shortwave_flux"]
        self._Rnet = self._cell_values["radiation__net_shortwave_flux"]
        # Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn 25, (52)
        self._julian = np.floor(
            ((self.current_time) - np.floor(self.current_time)) * 365
        )
        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)

        self._es = 0.6108 * np.exp((17.27 * self._Tavg) / (237.7 + self._Tavg))

        # Actual Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 8, (38)
        self._ea = 0.6108 * np.exp((17.27 * self._Tmin) / (237.7 + self._Tmin))

        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es) / ((237.3 + self._Tavg) ** 2.0)

        # Solar Declination Angle - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 24,(51)
        self._sdecl = 0.409 * np.sin(((np.pi / 180.0) * self._julian) - 1.39)

        # Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 23,(50)
        self._dr = 1 + (0.033 * np.cos(np.pi / 180.0 * self._julian))

        """ The following section computes the ratio to flat surface radiation field """
        ####################################################################

        # Apply distance factor to solar constant
        # Correction of solar constant using distance,
        # Src: https://www.lehman.edu/academics/eggs/documents/IJOC-2001-Becker.pdf
        self._Io *= self._dr

        # Convert latitude to radians and store trig calculations
        self._phi = np.radians(self._latitude)  # Latitude in Radians
        self._sinLat = np.sin(self._phi)
        self._cosLat = np.cos(self._phi)

        # Get the hour angle using time of day
        self._tau = (self._t + 12.0) * np.pi / 12.0  # Hour angle

        # Calculate solar altitude using declination angle, hour angle, and azimuth
        self._alpha = np.arcsin(
            np.sin(self._sdecl) * self._sinLat
            + (np.cos(self._sdecl) * self._cosLat * np.cos(self._tau))
        )  # Solar Altitude/Angle

        if self._alpha <= 0.25 * np.pi / 180.0:  # If altitude is -ve,
            self._alpha = 0.25 * np.pi / 180.0  # sun is beyond the horizon

        # For less constant calculation, keeping solar altitude trig in fixed variables
        self._cosSA = np.cos(self._alpha)
        self._sinSA = np.sin(self._alpha)

        # Avoid div by zero
        if self._cosSA <= 0.0001:
            self._cosSA = 0.0001

        # self._hAboveSea = 300
        self._hAboveSea = self._max_elev
        self._temp = self._Tavg  # in C

        # To calculate ws - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 29,(61)
        self._x = 1.0 - (((np.tan(self._phi)) ** 2.0) * (np.tan(self._sdecl) ** 2.0))
        if self._x <= 0:
            self._x = 0.00001

        # Sunset Hour Angle - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 28,(60)
        self._ws = (np.pi / 2.0) - np.arctan(
            (-1 * np.tan(self._phi) * np.tan(self._sdecl)) / (self._x**0.5)
        )

        # DIRECT SOLAR RADIATION (I)
        # Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 21, (48)
        self._Rgl = (
            11.57
            * (24.0 / np.pi)
            * 4.92
            * self._dr
            * (
                (self._ws * np.sin(self._phi) * np.sin(self._sdecl))
                + (np.cos(self._phi) * np.cos(self._sdecl) * (np.sin(self._ws)))
            )
        )

        # Sun's Azimuth calculation code
        self._F = np.tan(self._alpha) * np.tan(self._phi) - (
            np.sin(self._sdecl) / (self._cosSA * self._cosLat)
        )

        # Rounding
        if self._F > 0.99999 or self._F < -0.99999:
            self._F = 0.99999

        #
        if self._t < 12.0:
            self._phisun = np.pi - np.arccos(self._F)
        else:
            self._phisun = np.pi + np.arccos(self._F)

        self._flat = np.cos(np.arctan(0)) * np.sin(self._alpha) + np.sin(
            np.arctan(0)
        ) * np.cos(self._alpha) * np.cos(
            self._phisun - 0
        )  # flat surface reference

        # flat surface total incoming shortwave radiation
        self._Rsflat = self._Rgl * self._flat

        # flat surface Net incoming shortwave radiation
        self._Rnetflat = (1 - self._A) * (1 - 0.65 * (self._N**2)) * self._Rsflat

        self._sloped = np.cos(self._slope) * self._sinSA + np.sin(
            self._slope
        ) * self._cosSA * np.cos(self._phisun - self._aspect)

        self._sloped[self._sloped <= 0.0] = 0.0

        # Ratio of cosine of solar incidence angle of sloped surface to that
        # of a flat surface
        self._radf = self._sloped / self._flat

        self._radf[self._radf <= 0.0] = 0.0
        self._radf[self._radf > 6.0] = 6.0

        ####################################################################
        """ Radiation flux field computations follow below: """

        # Use the same direct radiation calculation from before, copy variable values
        self._Rext = self._Rgl

        # Clear-sky Radiation Rcs - Lecture_Radiation_Oct_13_18, Page 37
        # will only be computed if turbidity and optical air mass are both user defined
        # Optical airmass must be >0 or Rcs2 is disregarded (log of zero or / 0 error)

        self._rcs2valid = True
        if self._m > 0:
            self._Rcs2 = self._Rext * np.exp(
                -self._n * self._m * (0.128 - 0.054 * np.log10(self._m))
            )
        else:
            self._rcs2valid = False

        # Clear-sky Solar Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 19, (47)
        # This one is more ideal, and does not need as many parameters
        self._Rcs = self._Rext * (0.75 + 2 * (10**-5) * self._hAboveSea)

        # Rcfactor, set to Rcs2 for empirical method, Rcs for accurate method
        # IF optical air mass and turbidity are both appropriately defined,
        # then use the more parameter-intensive formula. Otherwise, use the
        # regular convenient empirical method.
        if (
            self._n != self._default_assertion_val
            and self._m != self._default_assertion_val
            and self._rcs2valid
        ):
            # print("Using turbidity + airmass clearsky factor...")
            self._Rcfactor = self._Rcs2
        else:
            # print("Using empirical clearsky factor...")
            self._Rcfactor = self._Rcs

        # Assuming 2500 foot elevation
        # Equation 34
        self._P = 101.3 * ((293 - 0.0065 * self._avg_elev) / 293) ** 5.26

        # If it does work, set _Po to a pressure according to the current month's
        # average monthly pressure
        self._Po = 101.325

        # non-constant KT (adjustment coefficient) - Handout_4_Radiation_Balance, Page 3
        self._KT = 0.2 * (self._P / self._Po) ** 0.5

        # Incoming shortwave radiation - Lecture_Radiation_Oct_13_18, Page 42
        # Shortwave radiation should ideally be below clearsky radiation
        # (Rcfactor), if there is a case where the standard shortwave rad
        # formula yields a result greater than Rcfactor, set shortwave
        # radiation to the clearsky radiation itself.
        self._Rs = min(
            self._KT * self._Rext * np.sqrt(self._Tmax - self._Tmin), self._Rcfactor
        )

        # Net shortwave Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn (43)
        self._Rns = self._Rs * (1 - self._A)

        # Relative Cloudiness - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Page 20,35
        self._u = 1.35 * (self._Rs / self._Rcfactor) - 0.35

        # As according to ASCE-EWRI Task Committee Report, Page 19
        # Cloudiness should be within 0.05 and 1 so as to not
        # nullify or incorrectly calculate the net longwave radiation
        if self._u < 0.05:
            self._u = 0.05
        elif self._u > 1:
            self._u = 1.0

        # Net Longwave Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn (17)
        self._Rnl = (
            5.67
            * (10**-8)
            * ((self._Tmax + 273.15) ** 4 + (self._Tmin + 273.15) ** 4)
            / 2
            * (0.34 - (0.14 * np.sqrt(self._ea)))
            * self._u
        )

        # Net Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Page 17, Eqn 15
        self._Rn = self._Rns - self._Rnl

        self._Rnet = self._Rns * self._cell_values["radiation__ratio_to_flat_surface"]
        # Apply shortwave, longwave, net radiation calculations to
        # the grid slopes.
        self._cell_values["radiation__ratio_to_flat_surface"] = self._radf
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
