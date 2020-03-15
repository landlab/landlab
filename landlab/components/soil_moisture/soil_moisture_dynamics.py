
import numpy as np
import six
from landlab import Component

from ...utils.decorators import use_file_name_or_kwds


_VALID_METHODS = set(["Grid"])


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class SoilMoisture(Component):
    """
    Landlab component that simulates root-zone average soil moisture at each
    cell using inputs of potential evapotranspiration, live leaf area index,
    and vegetation cover.

    This component uses a single soil moisture layer and models soil moisture
    loss through transpiration by plants, evaporation by bare soil, and
    leakage. The solution of water balance is based on Laio et. al 2001. The
    component requires fields of initial soil moisture, rainfall input (if
    any), time to the next storm and potential transpiration.

    Ref: Laio, F., Porporato, A., Ridolfi, L., & Rodriguez-Iturbe, I. (2001).
    Plants in water-controlled ecosystems: active role in hydrologic processes
    and response to water stress: II. Probabilistic soil moisture dynamics.
    Advances in Water Resources, 24(7), 707-723.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.soil_moisture import SoilMoisture
    >>> grid = RasterModelGrid((5, 5), xy_spacing=(0.2, 0.2))
    >>> SoilMoisture.name
    'Soil Moisture'
    >>> sorted(SoilMoisture.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['soil_moisture__root_zone_leakage',
     'soil_moisture__saturation_fraction',
     'surface__evapotranspiration',
     'surface__runoff',
     'surface__runon',
     'vegetation__water_stress']
    >>> sorted(SoilMoisture.units) # doctest: +NORMALIZE_WHITESPACE
    [('rainfall__daily_depth', 'mm'),
     ('soil_moisture__initial_saturation_fraction', 'None'),
     ('soil_moisture__root_zone_leakage', 'mm'),
     ('soil_moisture__saturation_fraction', 'None'),
     ('surface__evapotranspiration', 'mm/d'),
     ('surface__potential_evapotranspiration_rate', 'mm/d'),
     ('surface__potential_evapotranspiration_rate__grass', 'mm/d'),
     ('surface__runoff', 'mm'),
     ('surface__runon', 'mm'),
     ('vegetation__cover_fraction', 'None'),
     ('vegetation__live_leaf_area_index', 'None'),
     ('vegetation__plant_functional_type', 'None'),
     ('vegetation__water_stress', 'None')]
    >>> grid.at_cell["vegetation__plant_functional_type"] = (
    ...            np.zeros(grid.number_of_cells, dtype=int))

    Let us look at an example where we don't consider runoff from
    upstream cells (runon_switch = 0 - default condition).

    >>> sm = SoilMoisture(grid)
    >>> sm.grid.number_of_cell_rows
    3
    >>> sm.grid.number_of_cell_columns
    3
    >>> sm.grid is grid
    True
    >>> import numpy as np
    >>> np.allclose(grid.at_cell["soil_moisture__saturation_fraction"], 0.)
    True
    >>> grid.at_cell["surface__potential_evapotranspiration_rate"]= np.array([
    ...        7.5, 2., 3.5, 0., 1., 5., 3., 1.8, 0.])
    >>> grid.at_cell["surface__potential_evapotranspiration_rate__grass"]= np.array([
    ...        7.5, 2., 3.5, 0., 1., 5., 3., 1.8, 0.])
    >>> grid.at_cell["soil_moisture__initial_saturation_fraction"]= (
    ...        np.full(grid.number_of_cells, 0.4))
    >>> grid.at_cell["vegetation__live_leaf_area_index"]= (
    ...        np.full(grid.number_of_cells, 0.8))
    >>> grid.at_cell["vegetation__cover_fraction"]= (
    ...        np.ones(grid.number_of_cells))
    >>> current_time = 0.5
    >>> grid.at_cell["rainfall__daily_depth"] = (
    ...        60. * np.ones(grid.number_of_cells))
    >>> current_time = sm.update(current_time)
    >>> np.allclose(
    ...        grid.at_cell["soil_moisture__saturation_fraction"],
    ...        np.array([ 0.546,  0.575,  0.567,
    ...                   0.579,  0.579,  0.559,
    ...                   0.570,  0.576,  0.579]), rtol=1e-02)
    True

    Now, let's look at an example where we consider runoff
    from upstream cells (runon_switch=1). Please note that
    you will need two set of extra inputs to use this method:
    (1) a numpy array of 'ordered_cells'; and
    (2) a nodal field of flow receivers, "flow__receiver_node".
    If you have a DEM or elevation data at each node of the
    domain, you can use the helpful function
    'get_ordered_cells_for_soil_moisture' from the file
    'funcs_for_including_runon_in_soil_moisture.py', that
    takes the grid and optionally the 'outlet_id' as inputs,
    and outputs the ordered_cells and the grid. This function
    also creates the required nodal field "flow__receiver_node"
    on the grid. Note that the above mentioned function expects the
    input elevation data to be passed as a nodal field,
    "topographic__elevation", of the grid.

    >>> grid_r = RasterModelGrid((5, 5), spacing=(10, 10))
    >>> grid_r.at_node["topographic__elevation"] = np.array([
    ...         5., 5., 5., 5., 5.,
    ...         5., 4., 3., 2., 5.,
    ...         5., 4.5, 3.5, 1.5, 5.,
    ...         5., 4., 3., 1., 5.,
    ...         5., 5., 5., 5., 0.])
    >>> from landlab.components.soil_moisture.funcs_for_including_runon_in_soil_moisture import get_ordered_cells_for_soil_moisture
    >>> ordered_cells, grid_r = (
    ...     get_ordered_cells_for_soil_moisture(grid_r, outlet_id=24))
    >>> grid_r.at_cell["vegetation__plant_functional_type"] = (
    ...            np.zeros(grid_r.number_of_cells, dtype=int))
    >>> grid_r.at_cell["surface__potential_evapotranspiration_rate"]= np.array([
    ...        7.5, 2., 3.5, 0., 1., 5., 3., 1.8, 0.])
    >>> grid_r.at_cell["surface__potential_evapotranspiration_rate__grass"]= np.array([
    ...        7.5, 2., 3.5, 0., 1., 5., 3., 1.8, 0.])
    >>> grid_r.at_cell["soil_moisture__initial_saturation_fraction"]= (
    ...        np.full(grid_r.number_of_cells, 0.4))
    >>> grid_r.at_cell["vegetation__live_leaf_area_index"]= (
    ...        np.full(grid_r.number_of_cells, 0.8))
    >>> grid_r.at_cell["vegetation__cover_fraction"]= (
    ...        np.ones(grid_r.number_of_cells))
    >>> current_time = 0.5
    >>> grid_r.at_cell["rainfall__daily_depth"] = (
    ...        60. * np.ones(grid_r.number_of_cells))
    >>> sm_r = SoilMoisture(grid_r, runon_switch=1, ordered_cells=ordered_cells)
    >>> current_time = sm_r.update(current_time)
    >>> np.allclose(
    ...        grid_r.at_cell["soil_moisture__saturation_fraction"],
    ...        np.array([ 0.546,  0.575,  0.567,
    ...                   0.579,  0.579,  0.559,
    ...                   0.570,  0.576,  0.579]), rtol=1e-02)
    True
    >>> np.all(grid_r.at_cell["surface__runon"] == 0)
    False

    We can also specify the duration of the storm (Tr), and
    the duration between the storms (Tb). Remember that the
    units for Tb and Tr are hours. Let's consider a storm
    with a depth of 20 mm, Tr = 2 hours, Tb = 480 hrs (20 days).

    >>> grid_r.at_cell["rainfall__daily_depth"] = (
    ...        20. * np.ones(grid_r.number_of_cells))
    >>> current_time = sm_r.update(current_time, Tb=480, Tr=2)
    >>> np.allclose(
    ...        grid_r.at_cell["soil_moisture__saturation_fraction"],
    ...        np.array([ 0.195,  0.503,  0.371,
    ...                   0.560,  0.560,  0.280,
    ...                   0.411,  0.521,  0.560]), rtol=1e-02)
    True

    Adding more examples for doctesting.

    >>> grid_r.at_cell["soil_moisture__initial_saturation_fraction"]= (
    ...        np.full(grid_r.number_of_cells, 0.9))
    >>> grid_r.at_cell["rainfall__daily_depth"] = (
    ...        250. * np.ones(grid_r.number_of_cells))
    >>> current_time = sm_r.update(current_time, Tb=1, Tr=1)
    >>> np.allclose(
    ...        grid_r.at_cell["soil_moisture__saturation_fraction"],
    ...        np.array([ 0.875,  0.876,  0.876,
    ...                   0.876,  0.876,  0.876,
    ...                   0.875,  0.876,  0.876]), rtol=1e-02)
    True
    >>> grid_r.at_cell["soil_moisture__initial_saturation_fraction"]= (
    ...        np.full(grid_r.number_of_cells, 0.5))
    >>> grid_r.at_cell["rainfall__daily_depth"] = (
    ...        50. * np.ones(grid_r.number_of_cells))
    >>> current_time = sm_r.update(current_time, Tb=2400, Tr=1)
    >>> np.allclose(
    ...        grid_r.at_cell["soil_moisture__saturation_fraction"],
    ...        np.array([ 0.106,  0.255,  0.123,
    ...                   0.559,  0.559,  0.111,
    ...                   0.137,  0.313,  0.559]), rtol=1e-02)
    True
    >>> grid_r.at_cell["soil_moisture__initial_saturation_fraction"]= (
    ...        np.array([ 0.5, 0.5, 0.33,
    ...                   0.2, 0.2, 0.2,
    ...                   0.05, 0.05, 0.5]))
    >>> grid_r.at_cell["rainfall__daily_depth"] = (
    ...        5. * np.ones(grid_r.number_of_cells))
    >>> current_time = sm_r.update(current_time, Tb=2400, Tr=1)
    >>> np.allclose(
    ...        grid_r.at_cell["soil_moisture__saturation_fraction"],
    ...        np.array([ 0.106,  0.230,  0.113,
    ...                   0.168,  0.168,  0.105,
    ...                   0.098,  0.098,  0.530]), rtol=1e-02)
    True
    """

    _name = "Soil Moisture"

    _input_var_names = (
        "vegetation__cover_fraction",
        "vegetation__live_leaf_area_index",
        "surface__potential_evapotranspiration_rate",
        "surface__potential_evapotranspiration_rate__grass",
        "soil_moisture__initial_saturation_fraction",
        "vegetation__plant_functional_type",
        "rainfall__daily_depth",
    )

    _output_var_names = (
        "vegetation__water_stress",
        "soil_moisture__saturation_fraction",
        "soil_moisture__root_zone_leakage",
        "surface__runoff",
        "surface__runon",
        "surface__evapotranspiration",
    )

    _var_units = {
        "vegetation__cover_fraction": "None",
        "vegetation__live_leaf_area_index": "None",
        "surface__potential_evapotranspiration_rate": "mm/d",
        "surface__potential_evapotranspiration_rate__grass": "mm/d",
        "vegetation__plant_functional_type": "None",
        "vegetation__water_stress": "None",
        "soil_moisture__saturation_fraction": "None",
        "soil_moisture__initial_saturation_fraction": "None",
        "soil_moisture__root_zone_leakage": "mm",
        "surface__runoff": "mm",
        "surface__runon": "mm",
        "surface__evapotranspiration": "mm/d",
        "rainfall__daily_depth": "mm",
    }

    _var_mapping = {
        "vegetation__cover_fraction": "cell",
        "vegetation__live_leaf_area_index": "cell",
        "surface__potential_evapotranspiration_rate": "cell",
        "surface__potential_evapotranspiration_rate__grass": "cell",
        "vegetation__plant_functional_type": "cell",
        "vegetation__water_stress": "cell",
        "soil_moisture__saturation_fraction": "cell",
        "soil_moisture__initial_saturation_fraction": "cell",
        "soil_moisture__root_zone_leakage": "cell",
        "surface__runoff": "cell",
        "surface__runon": "cell",
        "surface__evapotranspiration": "cell",
        "rainfall__daily_depth": "cell",
    }

    _var_doc = {
        "vegetation__cover_fraction":
            "fraction of land covered by vegetation",
        "vegetation__live_leaf_area_index":
            "one-sided green leaf area per unit ground surface area",
        "surface__potential_evapotranspiration_rate":
            "potential sum of evaporation and plant transpiration",
        "surface__potential_evapotranspiration_rate__grass":
            "potential sum of evaporation and grass transpiration, \
             for partitioning bare soil evapotranspiration rate",
        "vegetation__plant_functional_type":
            "classification of plants (int), grass=0, shrub=1, tree=2, \
             bare=3, shrub_seedling=4, tree_seedling=5",
        "vegetation__water_stress":
            "parameter that represents nonlinear effects of water deficit \
             on plants",
        "soil_moisture__saturation_fraction":
            "relative volumetric water content (theta) - limits=[0,1]",
        "soil_moisture__initial_saturation_fraction":
            "initial soil_moisture__saturation_fraction",
        "soil_moisture__root_zone_leakage":
            "leakage of water into deeper portions of the soil not accessible \
             to the plant",
        "surface__runoff":
            "infiltration excess runoff from ground surface",
        "surface__runon":
            "infiltration excess runon",
        "surface__evapotranspiration":
            "actual sum of evaporation and plant transpiration",
        "rainfall__daily_depth":
            "Rain in (mm) as a field, allowing spatio-temporal soil moisture \
             saturation analysis.",
    }

    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        ordered_cells=None,
        runon_switch=0,
        f_bare=0.7,
        soil_ew=0.1,
        intercept_cap_grass=1.,
        zr_grass=0.3,
        I_B_grass=20.,
        I_V_grass=24.,
        K_s_grass=42.,
        pc_grass=0.43,
        fc_grass=0.56,
        sc_grass=0.33,
        wp_grass=0.13,
        hgw_grass=0.1,
        beta_grass=13.8,
        LAI_max_grass=2.,
        LAIR_max_grass=2.88,
        intercept_cap_shrub=1.5,
        zr_shrub=0.5,
        I_B_shrub=20.,
        I_V_shrub=40.,
        K_s_shrub=42.,
        pc_shrub=0.43,
        fc_shrub=0.56,
        sc_shrub=0.24,
        wp_shrub=0.13,
        hgw_shrub=0.1,
        beta_shrub=13.8,
        LAI_max_shrub=2.,
        LAIR_max_shrub=2.,
        intercept_cap_tree=2.,
        zr_tree=1.3,
        I_B_tree=20.,
        I_V_tree=40.,
        K_s_tree=42.,
        pc_tree=0.43,
        fc_tree=0.56,
        sc_tree=0.22,
        wp_tree=0.15,
        hgw_tree=0.1,
        beta_tree=13.8,
        LAI_max_tree=4.,
        LAIR_max_tree=4.,
        intercept_cap_bare=1.,
        zr_bare=0.15,
        I_B_bare=20.,
        I_V_bare=20.,
        K_s_bare=42.,
        pc_bare=0.43,
        fc_bare=0.56,
        sc_bare=0.33,
        wp_bare=0.13,
        hgw_bare=0.1,
        beta_bare=13.8,
        LAI_max_bare=0.01,
        LAIR_max_bare=0.01,
        **kwds,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        runon_switch: int, optional
            To indicate whether runon needs to be considered (mm);
            0 - No runon, 1 - runon.
        ordered_cells: numpy.array, required if runon_switch = 1
            ordered_cells has the grid cells sorted in an order of descending
            channel length in a delineated watershed
         f_bare: float, optional
            Fraction to partition PET for bare soil (None).
        soil_ew: float, optional
            Residual Evaporation after wilting (mm/day).
        intercept_cap: float, optional
            Plant Functional Type (PFT) specific full canopy interception
            capacity.
        zr: float, optional
            Root depth (m).
        I_B: float, optional
            Infiltration capacity of bare soil (mm/h).
        I_V: float, optional
            Infiltration capacity of vegetated soil (mm/h).
        K_s: float, optional
            Hydraulic conductivity of soil (mm/h).
        pc: float, optional
            Soil porosity (None).
        fc: float, optional
            Soil saturation degree at field capacity (None).
        sc: float, optional
            Soil saturation degree at stomatal closure (None).
        wp: float, optional
            Soil saturation degree at wilting point (None).
        hgw: float, optional
            Soil saturation degree at hygroscopic point (None).
        beta: float, optional
            Deep percolation constant = 2*b+4 where b is
            water retention (None).
        LAI_max: float, optional
            Maximum leaf area index (m^2/m^2).
        LAIR_max: float, optional
            Reference leaf area index (m^2/m^2).
        """
        self._method = kwds.pop("method", "Grid")

        assert_method_is_valid(self._method)

        super(SoilMoisture, self).__init__(grid, **kwds)

        self.initialize(
            ordered_cells=ordered_cells,
            runon_switch=runon_switch,
            f_bare=f_bare,
            soil_ew=soil_ew,
            intercept_cap_grass=intercept_cap_grass,
            zr_grass=zr_grass,
            I_B_grass=I_B_grass,
            I_V_grass=I_V_grass,
            K_s_grass=K_s_grass,
            pc_grass=pc_grass,
            fc_grass=fc_grass,
            sc_grass=sc_grass,
            wp_grass=wp_grass,
            hgw_grass=hgw_grass,
            beta_grass=beta_grass,
            LAI_max_grass=LAI_max_grass,
            LAIR_max_grass=LAIR_max_grass,
            intercept_cap_shrub=intercept_cap_shrub,
            zr_shrub=zr_shrub,
            I_B_shrub=I_B_shrub,
            I_V_shrub=I_V_shrub,
            K_s_shrub=K_s_shrub,
            pc_shrub=pc_shrub,
            fc_shrub=fc_shrub,
            sc_shrub=sc_shrub,
            wp_shrub=wp_shrub,
            hgw_shrub=hgw_shrub,
            beta_shrub=beta_shrub,
            LAI_max_shrub=LAI_max_shrub,
            LAIR_max_shrub=LAIR_max_shrub,
            intercept_cap_tree=intercept_cap_tree,
            zr_tree=zr_tree,
            I_B_tree=I_B_tree,
            I_V_tree=I_V_tree,
            K_s_tree=K_s_tree,
            pc_tree=pc_tree,
            fc_tree=fc_tree,
            sc_tree=sc_tree,
            wp_tree=wp_tree,
            hgw_tree=hgw_tree,
            beta_tree=beta_tree,
            LAI_max_tree=LAI_max_tree,
            LAIR_max_tree=LAIR_max_tree,
            intercept_cap_bare=intercept_cap_bare,
            zr_bare=zr_bare,
            I_B_bare=I_B_bare,
            I_V_bare=I_V_bare,
            K_s_bare=K_s_bare,
            pc_bare=pc_bare,
            fc_bare=fc_bare,
            sc_bare=sc_bare,
            wp_bare=wp_bare,
            hgw_bare=hgw_bare,
            beta_bare=beta_bare,
            LAI_max_bare=LAI_max_bare,
            LAIR_max_bare=LAIR_max_bare,
            **kwds
        )
        for name in self._input_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros("cell", name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros("cell", name, units=self._var_units[name])

        self._nodal_values = self.grid["node"]

        self._cell_values = self.grid["cell"]

    def initialize(
        self,
        ordered_cells=None,
        runon_switch=0,
        f_bare=0.7,
        soil_ew=0.1,
        intercept_cap_grass=1.,
        zr_grass=0.3,
        I_B_grass=20.,
        I_V_grass=24.,
        K_s_grass=42.,
        pc_grass=0.43,
        fc_grass=0.56,
        sc_grass=0.33,
        wp_grass=0.13,
        hgw_grass=0.1,
        beta_grass=13.8,
        LAI_max_grass=2.,
        LAIR_max_grass=2.88,
        intercept_cap_shrub=1.5,
        zr_shrub=0.5,
        I_B_shrub=20.,
        I_V_shrub=40.,
        K_s_shrub=42.,
        pc_shrub=0.43,
        fc_shrub=0.56,
        sc_shrub=0.24,
        wp_shrub=0.13,
        hgw_shrub=0.1,
        beta_shrub=13.8,
        LAI_max_shrub=2.,
        LAIR_max_shrub=2.,
        intercept_cap_tree=2.,
        zr_tree=1.3,
        I_B_tree=20.,
        I_V_tree=40.,
        K_s_tree=42.,
        pc_tree=0.43,
        fc_tree=0.56,
        sc_tree=0.22,
        wp_tree=0.15,
        hgw_tree=0.1,
        beta_tree=13.8,
        LAI_max_tree=4.,
        LAIR_max_tree=4.,
        intercept_cap_bare=1.,
        zr_bare=0.15,
        I_B_bare=20.,
        I_V_bare=20.,
        K_s_bare=42.,
        pc_bare=0.43,
        fc_bare=0.56,
        sc_bare=0.33,
        wp_bare=0.13,
        hgw_bare=0.1,
        beta_bare=13.8,
        LAI_max_bare=0.01,
        LAIR_max_bare=0.01,
        **kwds
    ):
        # GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
        # SHRUBSEEDLING = 4; TREESEEDLING = 5
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        runon_switch: int, optional
            To indicate whether runon needs to considered (mm);
            0 - No runon, 1 - runon.
        ordered_cells: numpy.array, required if runon_switch = 1
            ordered_cells has the grid cells sorted in an order of descending
            channel length in a delineated watershed
        f_bare: float, optional
            Fraction to partition PET for bare soil (None).
        soil_ew: float, optional
            Residual Evaporation after wilting (mm/day).
        intercept_cap: float, optional
            Plant Functional Type (PFT) specific full canopy interception
            capacity.
        zr: float, optional
            Root depth (m).
        I_B: float, optional
            Infiltration capacity of bare soil (mm/h).
        I_V: float, optional
            Infiltration capacity of vegetated soil (mm/h).
        K_s: float, optional
            Hydraulic conductivity of soil (mm/h).
        pc: float, optional
            Soil porosity (None).
        fc: float, optional
            Soil saturation degree at field capacity (None).
        sc: float, optional
            Soil saturation degree at stomatal closure (None).
        wp: float, optional
            Soil saturation degree at wilting point (None).
        hgw: float, optional
            Soil saturation degree at hygroscopic point (None).
        beta: float, optional
            Deep percolation constant = 2*b+4 where b is
            water retention (None).
        parameter (None)
        LAI_max: float, optional
            Maximum leaf area index (m^2/m^2).
        LAIR_max: float, optional
            Reference leaf area index (m^2/m^2).
        """

        self._vegtype = self.grid.at_cell["vegetation__plant_functional_type"]
        self._runon_switch = runon_switch
        self._fbare = f_bare
        self.ordered_cells = ordered_cells
        self._interception_cap = np.choose(
            self._vegtype,
            [
                intercept_cap_grass,
                intercept_cap_shrub,
                intercept_cap_tree,
                intercept_cap_bare,
                intercept_cap_shrub,
                intercept_cap_tree,
            ],
        )
        self._zr = np.choose(
            self._vegtype,
            [zr_grass, zr_shrub, zr_tree, zr_bare, zr_shrub, zr_tree],
        )

        self._soil_Ib = np.choose(
            self._vegtype,
            [I_B_grass, I_B_shrub, I_B_tree, I_B_bare, I_B_shrub, I_B_tree],
        )

        self._soil_Iv = np.choose(
            self._vegtype,
            [I_V_grass, I_V_shrub, I_V_tree, I_V_bare, I_V_shrub, I_V_tree],
        )

        self._soil_Ks = np.choose(
            self._vegtype,
            [K_s_grass, K_s_shrub, K_s_tree, K_s_bare, K_s_shrub, K_s_tree],
        )

        self._soil_Ew = soil_ew

        self._soil_pc = np.choose(
            self._vegtype,
            [pc_grass, pc_shrub, pc_tree, pc_bare, pc_shrub, pc_tree],
        )

        self._soil_fc = np.choose(
            self._vegtype,
            [fc_grass, fc_shrub, fc_tree, fc_bare, fc_shrub, fc_tree],
        )

        self._soil_sc = np.choose(
            self._vegtype,
            [sc_grass, sc_shrub, sc_tree, sc_bare, sc_shrub, sc_tree],
        )

        self._soil_wp = np.choose(
            self._vegtype,
            [wp_grass, wp_shrub, wp_tree, wp_bare, wp_shrub, wp_tree],
        )

        self._soil_hgw = np.choose(
            self._vegtype,
            [hgw_grass, hgw_shrub, hgw_tree, hgw_bare, hgw_shrub, hgw_tree],
        )

        self._soil_beta = np.choose(
            self._vegtype,
            [
                beta_grass,
                beta_shrub,
                beta_tree,
                beta_bare,
                beta_shrub,
                beta_tree,
            ],
        )

        self._LAI_max = np.choose(
            self._vegtype,
            [
                LAI_max_grass,
                LAI_max_shrub,
                LAI_max_tree,
                LAI_max_bare,
                LAI_max_shrub,
                LAI_max_tree,
            ],
        )

        self._LAIR_max = np.choose(
            self._vegtype,
            [
                LAIR_max_grass,
                LAIR_max_shrub,
                LAIR_max_tree,
                LAIR_max_bare,
                LAIR_max_shrub,
                LAIR_max_tree,
            ],
        )

    def update(self, current_time, Tb=24., Tr=0., **kwds):
        """
        Update fields with current loading conditions.

        Parameters
        ----------
        current_time: float
            Current time (years).
        Tr: float, optional
            Storm duration (hours).
        Tb: float, optional
            Inter-storm duration (hours).
        """
        P_ = self._cell_values["rainfall__daily_depth"]
        self._PET = (
            self._cell_values["surface__potential_evapotranspiration_rate"])
        self._pet_g = (
            self._cell_values["surface__potential_evapotranspiration_rate__grass"])
        self._SO = (
            self._cell_values["soil_moisture__initial_saturation_fraction"])
        self._vegcover = self._cell_values["vegetation__cover_fraction"]
        self._water_stress = self._cell_values["vegetation__water_stress"]
        self._S = self._cell_values["soil_moisture__saturation_fraction"]
        self._D = self._cell_values["soil_moisture__root_zone_leakage"]
        self._ETA = self._cell_values["surface__evapotranspiration"]
        self._fr = (
            self._cell_values["vegetation__live_leaf_area_index"]
            / self._LAIR_max
        )
        self._runoff = self._cell_values["surface__runoff"]
        self._runon = self._cell_values["surface__runon"]
        self._runoff[:] = 0.   # Initializing runoff to zero
        self._runon[:] = 0.    # Initializing runon to zero

        # LAIl = self._cell_values["vegetation__live_leaf_area_index"]
        # LAIt = LAIl+self._cell_values["DeadLeafAreaIndex"]
        # if LAIt.all() == 0.:
        #     self._fr = np.zeros(self.grid.number_of_cells)
        # else:
        #     self._fr = (self._vegcover[0]*LAIl/LAIt)
        self._fr[self._fr > 1.0] = 1.0
        self._Sini = np.zeros(self._SO.shape)
        self._ETmax = np.zeros(self._SO.shape)
        self._ts = np.zeros(self._SO.shape)  # record Time to Saturation
        self._precip_int = np.zeros(self._SO.shape)
        # Adding routine to add runon & runoff
        if self._runon_switch:
            # Make sure that flow_router has been called before
            r_cell = (self.grid.cell_at_node[
                self.grid.at_node["flow__receiver_node"]])
            r_cell = r_cell[r_cell != -1]
            # finding cells that aren't part of
            # ordered cells (delineated watershed)
            diff_cells = (np.setdiff1d(
                range(0, self.grid.number_of_cells), self.ordered_cells))
            # trvrsl_order has ordered cells computed first and then the rest
            # of the cells
            trvrsl_order = np.concatenate((self.ordered_cells, diff_cells),
                                          axis=0)
        else:
            # trvrsl_order is just regular 'all cells' if no runon is computed
            trvrsl_order = range(0, self.grid.number_of_cells)

        for cell in trvrsl_order:  # Routine to calculate runon
            if self._runon_switch:
                if cell in self.ordered_cells:
                    donors = []
                    donors = list(np.where(r_cell == cell)[0])
                    if len(donors) != 0:
                        for k in range(0, len(donors)):
                            self._runon[cell] += self._runoff[donors[k]]
            P = P_[cell]
            runon = self._runon[cell]
            if runon < 0:
                six._print('Runon < 0!')
            # print cell
            s = self._SO[cell]
            fbare = self._fbare
            ZR = self._zr[cell]
            pc = self._soil_pc[cell]
            fc = self._soil_fc[cell]
            scc = self._soil_sc[cell]
            wp = self._soil_wp[cell]
            hgw = self._soil_hgw[cell]
            beta = self._soil_beta[cell]
            Ks = self._soil_Ks[cell]
            if self._vegtype[cell] == 0:   # 0 - GRASS
                sc = scc * self._fr[cell] + (1 - self._fr[cell]) * fc
            else:
                sc = scc

            # Infiltration capacity
            Inf_cap = (
                self._soil_Ib[cell] * (1 - self._vegcover[cell])
                + self._soil_Iv[cell] * self._vegcover[cell]
            )
            # Interception capacity
            Int_cap = min(
                self._vegcover[cell] * self._interception_cap[cell],
                P * self._vegcover[cell]
            )
            # Effective precipitation depth
            Peff = max((P + max(runon, 0.) - Int_cap), 0.)
            mu = (Ks / 1000.0) / (pc * ZR * (np.exp(beta * (1. - fc)) - 1.))

            if self._vegtype[cell] == 3:
                Ep = max((fbare * self._pet_g[cell]), 0.0001)
            else:
                Ep = max(
                    (
                        self._PET[cell] * self._fr[cell]
                        + fbare * self._pet_g[cell] * (1. - self._fr[cell])
                    )
                    - Int_cap,
                    0.0001,
                )  # mm/d
            self._ETmax[cell] = Ep
            nu = ((Ep / 24.) / 1000.) / (pc * ZR)   # Loss function parameter
            # Loss function parameter
            nuw = ((self._soil_Ew / 24.) / 1000.) / (pc * ZR)
            # Precipitation Intensity
            if Tr == 0:
                precip_int = 0.
            else:
                precip_int = Peff / Tr
            self._precip_int[cell] = precip_int
            # Time to saturation Ts
            if precip_int <= 0.:
                Ts = np.inf
            else:
                Ts = (
                    (1 - self._SO[cell]) * (pc * ZR * 1000.)
                    / (precip_int * (1 - np.exp((-1) * Inf_cap / precip_int)))
                )
            self._ts[cell] = Ts
            # Computing runoff
            # If using Poisson storms with Tr = 0, (Precip_int * Tr = Precip)
            if Tr == 0.:
                self._runoff[cell] = max((Peff - Inf_cap), 0.)
                sini = min(
                    self._SO[cell]
                    + ((Peff - self._runoff[cell]) / (pc * ZR * 1000.)),
                    1.,
                )
            # If using regular storms with (Tr != 0.)
            elif Tr < Ts:
                self._runoff[cell] = max((precip_int - Inf_cap) * Tr, 0.)
                sini = min(
                    self._SO[cell]
                    + (
                        (precip_int * Tr - self._runoff[cell])
                        / (pc * ZR * 1000.)
                    ),
                    1.,
                )
            else:
                sini = 1
                self._runoff[cell] = max(
                    (precip_int - Inf_cap) * Ts + (precip_int * (Tr - Ts)),
                    0.,
                )
            if sini >= fc:
                tfc = (1.0 / (beta * (mu - nu))) * (
                    beta * (fc - sini)
                    + np.log((nu - mu + mu * np.exp(beta * (sini - fc))) / nu)
                )
                tsc = ((fc - sc) / nu) + tfc
                twp = ((sc - wp) / (nu - nuw)) * np.log(nu / nuw) + tsc

                if Tb < tfc:
                    s = abs(
                        sini
                        - (1.0 / beta)
                        * np.log(
                            (
                                (nu - mu + mu * np.exp(beta * (sini - fc)))
                                * np.exp(beta * (nu - mu) * Tb)
                                - mu * np.exp(beta * (sini - fc))
                            )
                            / (nu - mu)
                        )
                    )

                    self._D[cell] = ((pc * ZR * 1000.0) * (sini - s)) - (
                        Tb * (Ep / 24.0)
                    )
                    self._ETA[cell] = Tb * (Ep / 24.0)

                elif tfc <= Tb < tsc:
                    s = fc - (nu * (Tb - tfc))
                    self._D[cell] = ((pc * ZR * 1000.0) * (sini - fc)) - (
                        (tfc) * (Ep / 24.0)
                    )
                    self._ETA[cell] = Tb * (Ep / 24.0)

                elif tsc <= Tb < twp:
                    s = wp + (sc - wp) * (
                        (nu / (nu - nuw))
                        * np.exp((-1) * ((nu - nuw) / (sc - wp)) * (Tb - tsc))
                        - (nuw / (nu - nuw))
                    )
                    self._D[cell] = ((pc * ZR * 1000.0) * (sini - fc)) - (
                        tfc * Ep / 24.0
                    )
                    self._ETA[cell] = (1000.0 * ZR * pc * (sini - s)) - self._D[cell]

                else:
                    s = hgw + (wp - hgw) * np.exp(
                        (-1) * (nuw / (wp - hgw)) * max(Tb - twp, 0.0)
                    )
                    self._D[cell] = ((pc * ZR * 1000.0) * (sini - fc)) - (
                        tfc * Ep / 24.0
                    )
                    self._ETA[cell] = (1000.0 * ZR * pc * (sini - s)) - self._D[cell]

            elif fc > sini >= sc:
                tfc = 0.0
                tsc = (sini - sc) / nu
                twp = ((sc - wp) / (nu - nuw)) * np.log(nu / nuw) + tsc

                if Tb < tsc:
                    s = sini - nu * Tb
                    self._D[cell] = 0.0
                    self._ETA[cell] = 1000.0 * ZR * pc * (sini - s)

                elif tsc <= Tb < twp:
                    s = wp + (sc - wp) * (
                        (nu / (nu - nuw))
                        * np.exp((-1) * ((nu - nuw) / (sc - wp)) * (Tb - tsc))
                        - (nuw / (nu - nuw))
                    )
                    self._D[cell] = 0
                    self._ETA[cell] = 1000.0 * ZR * pc * (sini - s)

                else:
                    s = hgw + (wp - hgw) * np.exp(
                        (-1) * (nuw / (wp - hgw)) * (Tb - twp)
                    )
                    self._D[cell] = 0.0
                    self._ETA[cell] = 1000.0 * ZR * pc * (sini - s)

            elif sc > sini >= wp:
                tfc = 0
                tsc = 0
                twp = ((sc - wp) / (nu - nuw)) * np.log(
                    1 + (nu - nuw) * (sini - wp) / (nuw * (sc - wp))
                )

                if Tb < twp:
                    s = wp + ((sc - wp) / (nu - nuw)) * (
                        (np.exp((-1) * ((nu - nuw) / (sc - wp)) * Tb))
                        * (nuw + ((nu - nuw) / (sc - wp)) * (sini - wp))
                        - nuw
                    )
                    self._D[cell] = 0.0
                    self._ETA[cell] = 1000.0 * ZR * pc * (sini - s)

                else:
                    s = hgw + (wp - hgw) * np.exp(
                        (-1) * (nuw / (wp - hgw)) * (Tb - twp)
                    )
                    self._D[cell] = 0.0
                    self._ETA[cell] = 1000.0 * ZR * pc * (sini - s)

            else:
                tfc = 0.0
                tsc = 0.0
                twp = 0.0

                s = hgw + (sini - hgw) * np.exp((-1) * (nuw / (wp - hgw)) * Tb)
                self._D[cell] = 0.0
                self._ETA[cell] = 1000.0 * ZR * pc * (sini - s)

            self._water_stress[cell] = min(
                ((max(((sc - (s + sini) / 2.0) / (sc - wp)), 0.0)) ** 4.0), 1.0
            )
            self._S[cell] = s
            self._SO[cell] = s
            self._Sini[cell] = sini

        current_time += (Tb + Tr) / (24.0 * 365.25)
        return current_time
