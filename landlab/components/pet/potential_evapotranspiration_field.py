"""Landlab component that simulates potential evapotranspiration rate.

Potential Evapotranspiration Component calculates spatially distributed
potential evapotranspiration based on input radiation factor (spatial
distribution of incoming radiation) using chosen method such as constant 
or Priestly Taylor. Ref: ASCE-EWRI Task Committee Report Jan 2005.

.. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.pet import PotentialEvapotranspiration

Create a grid on which to calculate incident shortwave radiation

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> PotentialEvapotranspiration.input_var_names
('radiation__ratio_to_flat_surface',)

Check the units for the fields.

>>> PotentialEvapotranspiration.var_units('radiation__ratio_to_flat_surface')
'None'

Create the input fields.

>>> grid['cell']['radiation__ratio_to_flat_surface'] = np.array([
...       0.38488566, 0.38488566,
...       0.33309785, 0.33309785,
...       0.37381705, 0.37381705])

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> PotentialEvapotranspiration.var_help('radiation__ratio_to_flat_surface')
name: radiation__ratio_to_flat_surface
description:
  ratio of total incident shortwave radiation on sloped surface
  to flat surface
units: None
at: cell
intent: in

Check the output variable names

>>> sorted(PotentialEvapotranspiration.output_var_names)
['radiation__incoming_shortwave',
 'radiation__net',
 'radiation__net_longwave',
 'radiation__net_shortwave',
 'surface__potential_evapotranspiration_rate']

Instantiate the 'PotentialEvapotranspiration' component to work on this grid,
and run it.

>>> PET = PotentialEvapotranspiration(grid, method='PriestlyTaylor')

Run the *update* method to update output variables with current time
>>> current_time = 0.5
>>> PET.update(current_time)

>>> PET.grid.at_cell['radiation__incoming_shortwave']
array([ 33.09968448,  33.09968448,  28.64599771,  28.64599771,
        32.14779789,  32.14779789])

>>> PET.grid.at_cell['radiation__net']
array([ 13.9764353 ,  13.9764353 ,  12.09585347,  12.09585347,
        13.57449849,  13.57449849])

>>> PET.grid.at_cell['radiation__net_shortwave']
array([ 13.23987379,  13.23987379,  11.45839908,  11.45839908,
        12.85911915,  12.85911915])

>>> PET.grid.at_cell['surface__potential_evapotranspiration_rate']
array([ 0.25488065,  0.25488065,  0.22058551,  0.22058551,  0.24755075,
        0.24755075])
"""
from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import numpy as np


_VALID_METHODS = set(['Constant', 'PriestlyTaylor', 'MeasuredRadiationPT',
                      'Cosine'])


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class PotentialEvapotranspiration(Component):

    """Calculate spatially distributed potential evapotranspiration.

    Potential Evapotranspiration Component calculates spatially distributed
    potential evapotranspiration based on input radiation factor (spatial
    distribution of incoming radiation) using chosen method such as constant 
    or Priestly Taylor. Ref: Xiaochi et. al. 2013 for 'Cosine' method and
    ASCE-EWRI Task Committee Report Jan 2005 for 'PriestlyTaylor' method.
    Note: Calling 'PriestlyTaylor' method would generate/overwrite shortwave &
    longwave radiation fields.

    Construction::

        PotentialEvapotranspiration(grid, method='Cosine',
            priestly_taylor_const=1.26, albedo=0.6,
            latent_heat_of_vaporization=28.34, psychometric_const=0.066,
            stefan_boltzmann_const=0.0000000567, solar_const=1366.67,
            latitude=34., elevation_of_measurement=300, adjustment_coeff=0.18,
            lt=0., nd=365., MeanTmaxF=12., delta_d=5., **kwds)
        
    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    method: {'Constant', 'PriestlyTaylor', 'MeasuredRadiationPT', 'Cosine'}, optional
        Priestly Taylor method will spit out radiation outputs too.
    priestly_taylor_constant: float, optional
        Alpha used in Priestly Taylor method.
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

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.pet import PotentialEvapotranspiration

    >>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    >>> PET = PotentialEvapotranspiration(grid)
    >>> PET.name
    'Potential Evapotranspiration'
    >>> PET.input_var_names
    ('radiation__ratio_to_flat_surface',)
    >>> sorted(PET.output_var_names)
    ['radiation__incoming_shortwave',
     'radiation__net',
     'radiation__net_longwave',
     'radiation__net_shortwave',
     'surface__potential_evapotranspiration_rate']
    >>> sorted(PET.units) # doctest: +NORMALIZE_WHITESPACE
    [('radiation__incoming_shortwave', 'W/m^2'),
     ('radiation__net', 'W/m^2'),
     ('radiation__net_longwave', 'W/m^2'),
     ('radiation__net_shortwave', 'W/m^2'),
     ('radiation__ratio_to_flat_surface', 'None'),
     ('surface__potential_evapotranspiration_rate', 'mm')]
    >>> PET.grid.number_of_cell_rows
    3
    >>> PET.grid.number_of_cell_columns
    2
    >>> PET.grid is grid
    True
    >>> pet_rate = grid.at_cell['surface__potential_evapotranspiration_rate']
    >>> np.allclose(pet_rate, 0.)
    True
    >>> grid['cell']['radiation__ratio_to_flat_surface'] = np.array([
    ...       0.38488566, 0.38488566,
    ...       0.33309785, 0.33309785,
    ...       0.37381705, 0.37381705])
    >>> current_time = 0.5
    >>> PET.update(current_time)
    >>> np.allclose(pet_rate, 0.)
    False
    """

    _name = 'Potential Evapotranspiration'

    _input_var_names = (
        'radiation__ratio_to_flat_surface',
    )

    _output_var_names = (
        'surface__potential_evapotranspiration_rate',
        'radiation__incoming_shortwave',
        'radiation__net_shortwave',
        'radiation__net_longwave',
        'radiation__net',
    )

    _var_units = {
        'radiation__ratio_to_flat_surface': 'None',
        'surface__potential_evapotranspiration_rate': 'mm',
        'radiation__incoming_shortwave': 'W/m^2',
        'radiation__net_shortwave': 'W/m^2',
        'radiation__net_longwave': 'W/m^2',
        'radiation__net': 'W/m^2',
    }

    _var_mapping = {
        'radiation__ratio_to_flat_surface': 'cell',
        'surface__potential_evapotranspiration_rate': 'cell',
        'radiation__incoming_shortwave': 'cell',
        'radiation__net_shortwave': 'cell',
        'radiation__net_longwave': 'cell',
        'radiation__net': 'cell',
    }

    _var_doc = {
        'radiation__ratio_to_flat_surface':
            'ratio of total incident shortwave radiation on sloped surface \
             to flat surface',
        'radiation__incoming_shortwave':
            'total incident shortwave radiation over the time step',
        'radiation__net_shortwave':
            'net incident shortwave radiation over the time step',
        'radiation__net_longwave':
            'net incident longwave radiation over the time step',
        'radiation__net':
            'net total radiation over the time step',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, method='Cosine', priestly_taylor_const=1.26,
                 albedo=0.6, latent_heat_of_vaporization=28.34,
                 psychometric_const=0.066, stefan_boltzmann_const=0.0000000567,
                 solar_const=1366.67, latitude=34.,
                 elevation_of_measurement=300, adjustment_coeff=0.18,
                 lt=0., nd=365., MeanTmaxF=12., delta_d=5., **kwds):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        method: {'Constant', 'PriestlyTaylor', 'MeasuredRadiationPT',
                 'Cosine'}, optional
            Priestly Taylor method will spit out radiation outputs too.
        priestly_taylor_constant: float, optional
            Alpha used in Priestly Taylor method.
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
        """

        self._method = method
        # For Priestly Taylor
        self._alpha = priestly_taylor_const
        self._a = albedo
        self._pwhv = latent_heat_of_vaporization
        self._y = psychometric_const
        self._sigma = stefan_boltzmann_const
        self._Gsc = solar_const
        self._phi = (np.pi / 180.) * latitude
        self._z = elevation_of_measurement
        self._Krs = adjustment_coeff
        self._LT = lt
        self._ND = nd
        self._TmaxF_mean = MeanTmaxF
        self._DeltaD = delta_d

        _assert_method_is_valid(self._method)

        super(PotentialEvapotranspiration, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        self._cell_values = self.grid['cell']

    def update(self, current_time, const_potential_evapotranspiration=12.,
               Tmin=0., Tmax=1., Tavg=0.5, obs_radiation=350., **kwds):
        """Update fields with current conditions."""

        if self._method == 'Constant':
            self._PET_value = const_potential_evapotranspiration
        elif self._method == 'PriestlyTaylor':
            Tmin = Tmin
            Tmax = Tmax
            Tavg = Tavg
            self._PET_value = (
                self._PriestlyTaylor(current_time, Tmax, Tmin, Tavg))
            self._cell_values['radiation__incoming_shortwave'] = (
                self._Rs *
                self._cell_values['radiation__ratio_to_flat_surface'])
            self._cell_values['radiation__net_shortwave'] = (
                self._Rns *
                self._cell_values['radiation__ratio_to_flat_surface'])
            self._cell_values['radiation__net_longwave'] = (
                self._Rnl *
                self._cell_values['radiation__ratio_to_flat_surface'])
            self._cell_values['radiation__net'] = (
                self._Rn *
                self._cell_values['radiation__ratio_to_flat_surface'])
        elif self._method == 'MeasuredRadiationPT':
            Tavg = Tavg
            Robs = obs_radiation
            self._PET_value = self._MeasuredRadPT(Tavg, (1-self._a)*Robs)
        elif self._method == 'Cosine':
            self._J = np.floor((current_time - np.floor(current_time)) * 365.)
            self._PET_value = (
                max((self._TmaxF_mean + self._DeltaD / 2. *
                     np.cos((2 * np.pi) * (self._J - self._LT - self._ND / 2) /
                            self._ND)), 0.0))

        self._PET = (
            self._PET_value *
            self._cell_values['radiation__ratio_to_flat_surface'])
        self._cell_values['surface__potential_evapotranspiration_rate'][:] = (
            self._PET)

    def _PriestlyTaylor(self, current_time, Tmax, Tmin, Tavg):

        # Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn 25, (52)
        self._J = np.floor((current_time - np.floor(current_time)) * 365)

        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)
        self._es = 0.6108 * np.exp((17.27 * Tavg)/(237.7 + Tavg))

        # Actual Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 8, (38)
        self._ea = 0.6108 * np.exp((17.27 * Tmin)/(237.7 + Tmin))

        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es)/((237.3 + Tavg) ** 2.0)

        # Solar Declination Angle - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 24,(51)
        self._sdecl = 0.409 * np.sin(((np.pi / 180.0) * self._J) - 1.39)

        # Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 23,(50)
        self._dr = 1 + (0.033 * np.cos(np.pi / 180.0 * self._J))

        # To calculate ws - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 29,(61)
        self._x = 1.0 - (((np.tan(self._phi)) ** 2.0) *
                         (np.tan(self._sdecl) ** 2.0))
        if self._x <= 0:
            self._x = 0.00001
            # Sunset Hour Angle - ASCE-EWRI Task Committee Report,
            # Jan-2005 - Eqn 28,(60)
        self._ws = ((np.pi / 2.0) -
                    np.arctan((- 1 * np.tan(self._phi) *
                               np.tan(self._sdecl)) / (self._x ** 2.0)))

        # Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 21, (48)
        self._Ra = (
                    11.57 * (24.0 / np.pi) * 4.92 * self._dr *
                    ((self._ws * np.sin(self._phi) * np.sin(self._sdecl)) +
                     (np.cos(self._phi) * np.cos(self._sdecl) *
                     (np.sin(self._ws)))))

        # Clear-sky Solar Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 19, (47)
        self._Rso = (0.75 + ((2.0 * (10 ** (- 5.0))) * self._z)) * self._Ra
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
        self._Rnl = (self._sigma *
                     self._fcd * (0.34 - (0.14 * np.sqrt(self._ea)) *
                                  (((Tmax + 273.16) ** 4.0 +
                                    (Tmin + 273.16) ** 4.0) / 2.0)))

        # Net Radiation - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 15, (42)
        self._Rn = self._Rns - self._Rnl

        self._ETp = max(self._alpha * (self._delta / (self._delta + self._y)) *
                        (self._Rn / self._pwhv), 0)

        return self._ETp

    def _MeasuredRadPT(self, Tavg, Rnobs):
        # Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 6, (37)
        self._es = 0.6108 * np.exp((17.27 * Tavg)/(237.7 + Tavg))

        # Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report,
        # Jan-2005 - Eqn 5, (36)
        self._delta = (4098.0 * self._es) / ((237.3 + Tavg) ** 2.0)
        self._ETp = max(self._alpha * (self._delta / (self._delta + self._y)) *
                        (Rnobs / self._pwhv), 0)

        return self._ETp
