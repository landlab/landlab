import copy

import numpy as np

from landlab import Component
from landlab.grid.mappers import map_node_to_cell

_VALID_METHODS = {"Grid"}


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class Radiation(Component):
    """Compute 1D and 2D daily incident shortwave radiation.

    Landlab component that computes 1D and 2D daily extraterrestiral, clear-sky,
    incident shortwave, net shortwave, longwave, and net radiation. This code also
    computes relative incidence shortwave radiation compared to a flat surface
    calculated at noon.

    **References**

    Bras, R. L.: Hydrology: an introduction to hydrologic science, Addison
    Wesley Publishing Company, Boston, Mass., USA, 643 pp., 1990.

    ASCE-EWRI (2005) The ASCE standardized reference evapotranspiration equation.
    In: Allen RG, Walter IA, Elliot RL et al (eds) Environmental and Water Resources
    Institute (EWRI) of the American Society of Civil Engineers, ASCE,
    Standardization of Reference Evapotranspiration Task Committee final report.
    American Society of Civil Engineers (ASCE), Reston

    Allen, R.G., 1996. Assessing integrity of weather data for  reference
    evapotranspiration estimation. J. Irrig. Drain.  Eng., ASCE 122 (2), 97-106.

    Flores-Cervantes, J.H., E. Istanbulluoglu, E.R. Vivoni, and R.L. Bras (2014).
    A geomorphic perspective on terrain-modulate organization of vegetation productivity:
    Analysis in two semiarid grassland ecosystems in Southwestern United States.
    Ecohydrol., 7: 242-257. doi: 10.1002/eco.1333.

    .. codeauthor:: Sai Nudurupati & Erkan Istanbulluoglu & Berkan Mertan


    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Radiation
    >>> import numpy as np

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    >>> z = grid.add_zeros("node", "topographic__elevation")
    >>> rad = Radiation(grid)

    >>> grid.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [1.0, 1.0, 1.0, 1.0],
    ...     [2.0, 2.0, 2.0, 2.0],
    ...     [3.0, 4.0, 4.0, 3.0],
    ...     [4.0, 4.0, 4.0, 4.0],
    ... ]
    >>> rad.current_time = 0.5
    >>> rad.update()

    >>> grid.at_cell["radiation__net_shortwave_flux"].reshape((3, 2))
    array([[251.63813643, 251.63813643],
           [251.62345462, 251.62345462],
           [251.59409258, 251.59409258]])
    >>> grid.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [100.0, 100.0, 100.0, 100.0],
    ...     [200.0, 200.0, 200.0, 200.0],
    ...     [300.0, 400.0, 400.0, 300.0],
    ...     [400.0, 400.0, 400.0, 400.0],
    ... ]
    >>> calc_rad = Radiation(grid, current_time=0.0, kt=0.2)
    >>> calc_rad.update()

    >>> grid.at_cell["radiation__net_shortwave_flux"].reshape((3, 2))
    array([[188.10745478, 188.10745478],
           [187.84329564, 187.69076199],
           [183.82445291, 183.41439585]])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed


    """

    _name = "Radiation"

    _unit_agnostic = False

    _info = {
        "radiation__extraterrestrial_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "extraterrestrial radiation",
        },
        "radiation__clearsky_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W/m^2",
            "mapping": "cell",
            "doc": "clearsky radiation",
        },
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
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": (
                "ratio of incident shortwave radiation on sloped "
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
        kt=0.17,
        clearsky_turbidity=None,
        opt_airmass=None,
        current_time=0.0,
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
        kt: float, optional
            Regional coefficient applied to actual KT coefficient (0.15-0.2). Default is 0.15
        clearsky_turbidity: float, optional
            Clear sky turbidity.
        opt_airmass: float, optional
            Optical air mass.
        max_daily_temp: float, optional
            Maximum daily temperature (Celsius)
        min_daily_Temp: float, optional
            Minimum daily temperature (Celsius)
        current_time: float
              Current time (years).
        """
        super().__init__(grid)

        self.current_time = current_time
        self._hour = 12

        self._method = method
        self._N = cloudiness
        self._latitude = self._validate_latitude(latitude)
        self._A = self._validate_albedo(albedo)

        # note that kt provided by the user is just a
        # 0.15-0.2 range value meant to indicate the type of region
        # where the model is run, e.g 0.2 for coastal regions,
        # 0.17 for interior regions, where the default is 0.17.
        # this parameter can be used as a calibration coefficient
        self._kt = kt

        self._n = clearsky_turbidity
        self._m = opt_airmass

        # For computations requiring temperature
        self._Tmin, self._Tmax = self._validate_temperature_range(
            min_daily_temp, max_daily_temp
        )

        _assert_method_is_valid(self._method)

        self.initialize_output_fields()

        if "Slope" not in self._grid.at_cell:
            self._grid.add_zeros("Slope", at="cell", units="radians")

        if "Aspect" not in self._grid.at_cell:
            self._grid.add_zeros("Aspect", at="cell", units="radians")

        self._nodal_values = self._grid["node"]
        self._cell_values = self._grid["cell"]

        self._slope, self._aspect = grid.calculate_slope_aspect_at_nodes_burrough(
            vals="topographic__elevation"
        )

        self._cell_values["Slope"] = self._slope
        self._cell_values["Aspect"] = self._aspect

        # Create a 'status' nodal field to store the status_at_node values
        # then map it to a cell-based field the same statuses. Use this to
        # generate a "closed elevations" field of boolean indexing.self
        self._gridCopy = copy.deepcopy(self._grid)
        self._gridCopy.add_field(
            "radiation_status_at_node", self._grid.status_at_node, at="node"
        )
        self._cellular_status = map_node_to_cell(
            self._gridCopy, "radiation_status_at_node"
        )
        self._closed_elevations = (
            self._cellular_status == self._gridCopy.BC_NODE_IS_CLOSED
        )

    def run_one_step(self, dt=None):
        if dt is None:
            dt = 1.0 / 365.0
        self.current_time += dt
        self.update()

    def _validate_latitude(self, latitude):
        if latitude < -90.0 or latitude > 90.0:
            raise ValueError("latitude must be between -90 and 90 degrees")
        return latitude

    def _validate_albedo(self, albedo):
        if albedo < 0.0 or albedo > 1.0:
            raise ValueError("albedo must be between 0 and 1")
        return albedo

    def _validate_temperature_range(self, min_temp, max_temp):
        if min_temp > max_temp:
            raise ValueError(
                f"minimum temperature ({min_temp}) must be less than maximum ({max_temp})"
            )
        return min_temp, max_temp

    @property
    def day_of_year(self):
        return (self.current_time - np.floor(self.current_time)) * 365

    @property
    def solar_declination(self):
        return 0.409 * np.sin(2.0 * np.pi / 365.0 * self.day_of_year - 1.39)

    @property
    def relative_distance_factor(self):
        return 1 + (0.033 * np.cos(2.0 * np.pi / 365.0 * self.day_of_year))

    @property
    def actual_vapor_pressure(self):
        return 0.6108 * np.exp((17.27 * self._Tmin) / (237.7 + self._Tmin))

    def update(self):
        """Update fields with current loading conditions.

        This method looks to the properties ``current_time`` and
        ``hour`` and uses their values in updating fields.
        """
        self._t = self._hour

        # Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn 25, (52)
        self._julian = (self.current_time - np.floor(self.current_time)) * 365

        # Actual Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 8, (38)
        self._ea = 0.6108 * np.exp((17.27 * self._Tmin) / (237.7 + self._Tmin))

        # Solar Declination Angle - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 24,(51)
        self._sdecl = 0.409 * np.sin(((np.pi / 180.0) * self._julian) - 1.39)

        # Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 23,(50)
        self._dr = 1 + (0.033 * np.cos(np.pi / 180.0 * self._julian))

        # Generate spatially distributed field of flat surface to
        # sloped surface radiation incidence ratios
        self._ratio_flat_surface_calc()

        # Extraterrestrial radmodel.docx - ASCE-EWRI (2005), Eqn (21)
        self._Rext = (
            11.57
            * (24.0 / np.pi)
            * 4.92
            * self._dr
            * (
                (self._ws * np.sin(self._phi) * np.sin(self._sdecl))
                + (np.cos(self._phi) * np.cos(self._sdecl) * (np.sin(self._ws)))
            )
        )

        # Clear-sky Radiation, Rcs, ASCE-EWRI (2005) as default method
        # if (optionally) turbidity and optical air mass are both user defined,
        # an exponential decay model Bras (2.25) is used.
        # Optical airmass must be >0 or Rcs2 is disregarded (log of zero or / 0 error)
        self._rcs2valid = True
        if self._m is not None:
            self._Rcs2 = self._Rext * np.exp(
                -self._n * self._m * (0.128 - 0.054 * np.log10(self._m))
            )
        else:
            self._rcs2valid = False

        # Clear-sky Solar Radiation - ASCE-EWRI (2005), Eqn 19
        self._Rcs1 = self._Rext * (0.75 + 2 * (10**-5) * self._elevation)

        # Rcs1, set to Rcs2 for empirical method, Rcs for accurate method
        # Using optical air mass and turbidity is optional when watershed
        # is relatively flat. n and m are not required fields at all
        if self._n is not None and self._m is not None and self._rcs2valid:
            self._Rc = self._Rcs2
        else:
            self._Rc = self._Rcs1

        # KT adjustment factor for incoming short-wave calculations
        # this uses the ratio of local to sealevel barometric
        # pressure following Allen (1996)
        self._Po = 101.325

        # Local atmospheric pressure ASCE-EWRI (2005), Eqn (34)
        self._P = self._Po * ((293 - 0.0065 * self._elevation) / 293) ** 5.26

        # non-constant KT (adjustment coefficient) - Based on Allen (1995)
        self._KT = self._kt * (self._P / self._Po) ** 0.5

        # Incoming shortwave radiation cannot exceed clear-sky radiation
        # Shortwave radiation should ideally be below clearsky radiation
        # (Rc), if there is a case where the standard shortwave rad
        # formula yields a result greater than Rc, set shortwave
        # radiation to the clearsky radiation itself.
        self._Rs = np.minimum(
            self._KT * self._Rext * np.sqrt(self._Tmax - self._Tmin), self._Rc
        )

        # Net shortwave Radiation - ASCE-EWRI (2005), Eqn (43)
        self._Rns = self._Rs * (1 - self._A)

        # Relative Cloudiness - ASCE-EWRI (2005), Eqn (18)
        self._u = 1.35 * (self._Rs / self._Rc) - 0.35

        # Cloudiness should be within 0.05 and 1 so as to not
        # nullify or incorrectly calculate the net longwave radiation
        self._u = np.clip(self._u, 0.05, 1.0)

        # Net Longwave Radiation - ASCE-EWRI (2005), Eqn (17) in W/M^2
        self._Rnl = (
            5.67
            * (10**-8)
            * ((self._Tmax + 273.15) ** 4 + (self._Tmin + 273.15) ** 4)
            / 2
            * (0.34 - (0.14 * np.sqrt(self._ea)))
            * self._u
        )

        # Load spatially distributed ratio to flat surface function values
        self._cell_values["radiation__ratio_to_flat_surface"] = self._radf

        # If sun does not rise, accounted by an invalid trig function
        # input, set all NaN grid values to 0.0 instead.
        self._radf[np.isnan(self._radf)] = 0.0
        self._Rext = 0.0 if np.isnan(self._Rext) else self._Rext

        # Net Radiation - ASCE-EWRI (2005), Eqn 15
        # Apply the ratio to flat surface to net shortwave
        # radiation first, then use spatially distributed
        # shortwave radiation to calculate net radiation
        np.multiply(
            self._Rext,
            self._cell_values["radiation__ratio_to_flat_surface"],
            out=self._cell_values["radiation__extraterrestrial_flux"],
        )

        # Clearsky flux
        np.multiply(
            self._Rc,
            self._cell_values["radiation__ratio_to_flat_surface"],
            out=self._cell_values["radiation__clearsky_flux"],
        )

        # Incoming shortwave flux
        np.multiply(
            self._Rs,
            self._cell_values["radiation__ratio_to_flat_surface"],
            out=self._cell_values["radiation__incoming_shortwave_flux"],
        )

        # Net shortwave flux
        np.multiply(
            self._Rns,
            self._cell_values["radiation__ratio_to_flat_surface"],
            out=self._cell_values["radiation__net_shortwave_flux"],
        )

        # Net longwave flux
        np.multiply(
            self._Rnl,
            self._cell_values["radiation__ratio_to_flat_surface"],
            out=self._cell_values["radiation__net_longwave_flux"],
        )

        # Net radiation flux is net shortwave - net longwave
        np.subtract(
            self._cell_values["radiation__net_shortwave_flux"],
            self._cell_values["radiation__net_longwave_flux"],
            out=self._cell_values["radiation__net_flux"],
        )

    def _ratio_flat_surface_calc(self):
        """generate radiation__ratio_to_flat_surface field

        This method looks to the slope, aspect values across
        the grid provided to the component, then runs calculations
        with solar altitude, latitude, and elevation to create a
        spatially distributed field of flat surface to sloped surface
        ratios / factors.
        """

        # self._elevation, elevation of surface above sea level
        self._elevation = self._nodal_values["topographic__elevation"]

        # Handle invalid values (closed nodes are not invalid)
        if np.any(
            (self._elevation < 0.0)
            & (self._grid.status_at_node != self._grid.BC_NODE_IS_CLOSED)
        ):
            raise ValueError(
                "No negative (< 0.0) values allowed in an above sea level elevation field."
            )

        # Map nodal elevation values to cells to calculate clearsky incidence
        # across a spatially distributed field
        self._elevation = map_node_to_cell(self._grid, "topographic__elevation")

        # Convert latitude to radians and store trig calculations
        self._phi = np.radians(self._latitude)  # Latitude in Radians
        self._sinLat = np.sin(self._phi)
        self._cosLat = np.cos(self._phi)

        # Get the hour angle using time of day
        self._tau = (self._t + 12.0) * np.pi / 12.0  # Hour angle

        # Calculate solar altitude using declination angle, hour angle, and latitude
        self._alpha = np.arcsin(
            np.sin(self._sdecl) * self._sinLat
            + (np.cos(self._sdecl) * self._cosLat * np.cos(self._tau))
        )  # Solar Altitude/Angle

        if self._alpha <= 0.25:  # If altitude is -ve,
            self._alpha = 0.25  # sun is beyond the horizon

        # For less constant calculation, keeping solar altitude trig in fixed variables
        self._cosSA = np.cos(self._alpha)
        self._sinSA = np.sin(self._alpha)

        # Avoid div by zero
        if self._cosSA <= 0.0001:
            self._cosSA = 0.0001

        # Sunset Hour Angle - ASCE-EWRI (2005), Eqn (59)
        self._ws = np.arccos(-np.tan(self._sdecl) * np.tan(self._phi))

        # Sun's Azimuth calculation code
        F = np.tan(self._alpha) * np.tan(self._phi) - (
            np.sin(self._sdecl) / (self._cosSA * self._cosLat)
        )

        # Clip azimuth within these bounds
        F = np.clip(F, -0.99999, 0.99999)

        if self._t < 12.0:
            self._phisun = np.pi - np.arccos(F)
        else:
            self._phisun = np.pi + np.arccos(F)

        self._flat = np.sin(self._alpha)

        # solar angle of incidence, Multiplying this with incoming radiation
        # gives radiation on sloped surface see Flores-Cervantes, J.H. (2012)
        self._sloped = np.cos(self._slope) * self._sinSA + np.sin(
            self._slope
        ) * self._cosSA * np.cos(self._phisun - self._aspect)

        self._sloped[self._sloped <= 0.0] = 0.0

        # Ratio of cosine of solar incidence angle of sloped surface to that
        # of a flat surface
        self._radf = self._sloped / self._flat

        self._radf[self._radf <= 0.0] = 0.0
        self._radf[self._radf > 6.0] = 6.0

        # Closed nodes will be omitted from spatially distributed ratio calculations
        self._radf[self._closed_elevations] = 0.0
