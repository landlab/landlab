#!/usr/bin/env python

import numpy as np

from landlab import Component


class SoilInfiltrationGreenAmpt(Component):
    """Infiltrate surface water into a soil following the Green-Ampt method.

    This component calculates the infiltation of surface water into the soil,
    using the Green-Ampt method. The component tracks the depth of infiltrated
    water over time, in the field soil_water_infiltration__depth. It also
    modifies the depth of surface water (surface_water__depth) as surface water
    progressively infiltrates into the soil below.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import SoilInfiltrationGreenAmpt
    >>> mg = RasterModelGrid((4, 3), xy_spacing=10.0)
    >>> hydraulic_conductivity = mg.ones("node") * 1.0e-6
    >>> hydraulic_conductivity.reshape((4, 3))[0:2, :] *= 10000.0
    >>> h = mg.add_ones("surface_water__depth", at="node")
    >>> h *= 0.01
    >>> d = mg.add_ones("soil_water_infiltration__depth", at="node", dtype=float)
    >>> d *= 0.2
    >>> SI = SoilInfiltrationGreenAmpt(
    ...     mg, hydraulic_conductivity=hydraulic_conductivity
    ... )
    >>> for i in range(10):  # 100s total
    ...     SI.run_one_step(10.0)
    ...
    >>> mg.at_node["surface_water__depth"]
    array([1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           9.88530416e-03, 9.88530416e-03, 9.88530416e-03,
           9.88530416e-03, 9.88530416e-03, 9.88530416e-03])
    >>> mg.at_node["soil_water_infiltration__depth"]
    array([0.20999999, 0.20999999, 0.20999999, 0.20999999, 0.20999999,
           0.20999999, 0.2001147 , 0.2001147 , 0.2001147 , 0.2001147 ,
           0.2001147 , 0.2001147 ])

    Notes
    -----
    This code is based on an overland flow model by Francis Rengers and
    colleagues, after Julien et al., 1995. The infiltration scheme follows the
    Green and Ampt equation. It was implemented in Landlab by DEJH, March 2016.

    **Where to learn more**

    A description of the Green-Ampt infiltration equation can be found in many
    hydrology texts, as well as online resources. The original theory was
    published by Green and Ampt (1911).

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Rengers, F. K., McGuire, L. A., Kean, J. W., Staley, D. M., and Hobley, D.:
    Model simulations of flood and debris flow timing in steep catchments after
    wildfire, Water Resour. Res., 52, 6041–6061, doi:10.1002/2015WR018176, 2016.

    **Additional References**

    Julien, P. Y., Saghaﬁan, B., and Ogden, F. L.: Raster-based hydrologic
    modeling of spatially-varied surface runoff, J. Am. Water Resour. As., 31,
    523–536, doi:10.1111/j.17521688.1995.tb04039.x, 1995.

    Green, W. H., & Ampt, G. A. (1911). Studies on Soil Phyics. The Journal of
    Agricultural Science, 4(1), 1-24.
    """

    _name = "SoilInfiltrationGreenAmpt"

    _unit_agnostic = False

    _cite_as = """
    @article{rengers2016model,
      author = {Rengers, F K and McGuire, L A and Kean, J W and Staley, D M
                and Hobley, D E J},
      title = {{Model simulations of flood and debris flow timing in steep
                catchments after wildfire}},
      doi = {10.1002/2015wr018176},
      pages = {6041 -- 6061},
      number = {8},
      volume = {52},
      journal = {Water Resources Research},
      year = {2016},
    }
    """

    _info = {
        "soil_water_infiltration__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "Water column height above the surface previously absorbed "
                "into the soil. Note that this is NOT the actual depth of "
                "the wetted front, which also depends on the porosity."
            ),
        },
        "surface_water__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
    }

    # This follows mean values from Rawls et al., 1992; lambda then h_b
    SOIL_PROPS = {
        "sand": (0.694, 0.0726),
        "loamy sand": (0.553, 0.0869),
        "sandy loam": (0.378, 0.1466),
        "loam": (0.252, 0.1115),
        "silt loam": (0.234, 0.2076),
        "sandy clay loam": (0.319, 0.2808),
        "clay loam": (0.242, 0.2589),
        "silty clay loam": (0.177, 0.3256),
        "sandy clay": (0.223, 0.2917),
        "silty clay": (0.150, 0.3419),
        "clay": (0.165, 0.3730),
    }

    def __init__(
        self,
        grid,
        hydraulic_conductivity=0.005,
        soil_bulk_density=1590.0,
        rock_density=2650.0,
        initial_soil_moisture_content=0.15,
        soil_type="sandy loam",
        volume_fraction_coarse_fragments=0.2,
        coarse_sed_flag=False,
        surface_water_minimum_depth=1.0e-8,
        soil_pore_size_distribution_index=None,
        soil_bubbling_pressure=None,
        wetting_front_capillary_pressure_head=None,
    ):
        """
        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        hydraulic_conductivity : float, array, or field name (m/s)
            The soil effective hydraulic conductivity.
        soil_bulk_density : float (kg/m**3)
            The dry bulk density of the soil.
        rock_density : float (kg/m**3)
            The density of the soil constituent material (i.e., lacking porosity).
        initial_soil_moisture_content : float (m**3/m**3, 0. to 1.)
            The fraction of the initial pore space filled with water.
        soil_type : str
            A soil type to automatically set soil_pore_size_distribution_index
            and soil_bubbling_pressure, using mean values from Rawls et al.,
            1992. The following options are supported: 'sand', loamy sand',
            'sandy loam', 'loam', 'silt loam', 'sandy clay loam', 'clay loam',
            'silty clay loam', 'sandy clay', 'silty clay', or 'clay'.
        volume_fraction_coarse_fragments : float (m**3/m**3, 0. to 1.)
            The fraction of the soil made up of rocky fragments with very
            little porosity, with diameter > 2 mm.
        coarse_sed_flag : boolean, optional
            If this flag is set to true, the fraction of coarse material in the
            soil column with be used as a correction for phi, the porosity factor.
        surface_water_minimum_depth : float (m), optional
            A minimum water depth to stabilize the solutions for surface flood
            modelling. Leave as the default in most normal use cases.
        soil_pore_size_distribution_index : float, optional
            An index describing the distribution of pore sizes in the soil,
            and controlling effective hydraulic conductivity at varying water
            contents, following Brooks and Corey (1964). Can be set by
            soil_type. Typically denoted "lambda".
        soil_bubbling_pressure : float (m), optional
            The bubbling capillary pressure of the soil, controlling effective
            hydraulic conductivity at varying water contents, following Brooks
            and Corey (1964). Can be set by soil_type. Typically denoted "h_b".
        wetting_front_capillary_pressure_head : float (m), optional
            The effective head at the wetting front in the soil driven by
            capillary pressure in the soil pores. If not set, will be
            calculated by the component from the pore size distribution and
            bubbling pressure, following Brooks and Corey.

        """
        super().__init__(grid)

        self._min_water = surface_water_minimum_depth
        self._hydraulic_conductivity = hydraulic_conductivity

        if not coarse_sed_flag:
            volume_fraction_coarse_fragments = 0.0

        self._moisture_deficit = self.calc_moisture_deficit(
            soil_bulk_density=soil_bulk_density,
            rock_density=rock_density,
            volume_fraction_coarse_fragments=volume_fraction_coarse_fragments,
            soil_moisture_content=initial_soil_moisture_content,
        )

        if wetting_front_capillary_pressure_head is None:
            self._capillary_pressure = self.calc_soil_pressure(
                soil_type=soil_type,
                soil_pore_size_distribution_index=soil_pore_size_distribution_index,
                soil_bubbling_pressure=soil_bubbling_pressure,
            )
        else:
            self._capillary_pressure = wetting_front_capillary_pressure_head

    @staticmethod
    def calc_soil_pressure(
        soil_type=None,
        soil_pore_size_distribution_index=1.0,
        soil_bubbling_pressure=0.0,
    ):
        """Calculate capillary pressure in a soil type.

        Parameters
        ----------
        soil_type : str, optional
            The name of a soil type.
        soil_pore_size_distribution_index : float
            Pore-size distribution index [-].
        soil_bubbling_pressure : float
            Bubbling pressure [m].
        """
        if soil_type is None:
            soil_props = (soil_pore_size_distribution_index, soil_bubbling_pressure)
        else:
            try:
                soil_props = SoilInfiltrationGreenAmpt.SOIL_PROPS[soil_type]
            except KeyError as exc:
                raise ValueError(f"{soil_type}: unknown soil type") from exc

        return SoilInfiltrationGreenAmpt.calc_pressure_head(*soil_props)

    @staticmethod
    def calc_pressure_head(lam, h_b):
        """Calculate pressure head.

        Pressure head is set using *lambda* and *h_b*, using an
        equation after Brooks-Corey (1964), following Rawls et al., 1992.

        Parameters
        ----------
        lam : float, optional
            Pore-size distribution index. Exponent that describes the
            distribution of pore sizes in the soil, and controls
            effective hydraulic conductivity at varying water
            contents, following Brooks and Corey (1964) [-].
        h_b : float (m), optional
            Bubbling pressure. Capillary pressure of the soil,
            controlling effective hydraulic conductivity at varying
            water contents, following Brooks and Corey (1964) [m]
        """
        return (2.0 + 3.0 * lam) / (1.0 + 3.0 * lam) * h_b * 0.5

    @staticmethod
    def calc_moisture_deficit(
        soil_bulk_density=1590.0,
        rock_density=2650.0,
        volume_fraction_coarse_fragments=0.0,
        soil_moisture_content=0.0,
    ):
        """Calculate the moisture deficit in a soil.

        Parameters
        ----------
        soil_bulk_density : float or array of float
            Bulk density of the soil [kg / m3].
        rock_density : float or array of float
            Density of rock [kg / m3].
        volume_fraction_coarse_fragments : float or array of float
            Volume fraction of sediment made up of coarse grains [-].
        soil_moisture_content : float or array of float
            Fraction of soil filled with water [-].

        Returns
        -------
        float or array of float
            Moisture deficit.
        """
        if np.any(soil_bulk_density <= 0.0):
            raise ValueError("non-positive soil bulk density")
        if np.any(rock_density < soil_bulk_density):
            raise ValueError("soil bulk density greater than rock density")
        if np.any(volume_fraction_coarse_fragments < 0.0):
            raise ValueError("negative volume fraction of coarse grains")
        if np.any(volume_fraction_coarse_fragments > 1.0):
            raise ValueError("volume fraction of coarse grains")

        soil_porosity = 1.0 - np.true_divide(soil_bulk_density, rock_density)
        soil_porosity *= 1.0 - volume_fraction_coarse_fragments

        if np.any(soil_moisture_content > soil_porosity):
            raise ValueError("soil moisture greater than porosity")

        return soil_porosity - soil_moisture_content

    @property
    def min_water(self):
        """Minimum surface water depth."""
        return self._min_water

    @min_water.setter
    def min_water(self, new_value):
        if np.any(new_value <= 0.0):
            raise ValueError("minimum water depth must be positive")
        self._min_water = new_value

    @property
    def hydraulic_conductivity(self):
        """Hydraulic conductivity of soil."""
        return self._hydraulic_conductivity

    @hydraulic_conductivity.setter
    def hydraulic_conductivity(self, new_value):
        if isinstance(new_value, str):
            new_value = self._grid.at_node[new_value]
        if np.any(new_value < 0.0):
            raise ValueError("hydraulic conductivity must be positive")
        self._hydraulic_conductivity = new_value

    @property
    def moisture_deficit(self):
        """Moisture deficit of soil."""
        return self._moisture_deficit

    @moisture_deficit.setter
    def moisture_deficit(self, new_value):
        if np.any(new_value < 0.0):
            raise ValueError("negative moisture deficit")
        self._moisture_deficit = new_value

    @property
    def capillary_pressure(self):
        """Capillary pressure of soil."""
        return self._capillary_pressure

    @capillary_pressure.setter
    def capillary_pressure(self, new_value):
        if np.any(new_value < 0.0):
            raise ValueError("negative capillary pressure")
        self._capillary_pressure = new_value

    def run_one_step(self, dt):
        """Update fields with current hydrologic conditions.

        Parameters
        ----------
        dt : float (s)
            The imposed timestep for the model.
        """
        water_depth = self._grid.at_node["surface_water__depth"]
        infiltration_depth = self._grid.at_node["soil_water_infiltration__depth"]

        assert np.all(infiltration_depth >= 0.0)

        wettingfront_depth = infiltration_depth / self._moisture_deficit

        potential_infilt = (
            dt
            * self._hydraulic_conductivity
            * (
                (wettingfront_depth + self._capillary_pressure + water_depth)
                / wettingfront_depth
            )
        )
        np.clip(potential_infilt, 0.0, None, out=potential_infilt)

        available_water = water_depth - self._min_water
        np.clip(available_water, 0.0, None, out=available_water)

        actual_infiltration = np.choose(
            potential_infilt > available_water, (potential_infilt, available_water)
        )

        water_depth -= actual_infiltration
        infiltration_depth += actual_infiltration
