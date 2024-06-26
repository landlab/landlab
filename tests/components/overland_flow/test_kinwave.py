"""
Unit tests for landlab.components.overland_flow.KinwaveOverlandFlowModel

last updated: 3/14/16
"""

(_SHAPE, _SPACING, _ORIGIN) = ((10, 10), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_KinWaveOF_name(kin_wave_of):
    assert kin_wave_of.name == "KinwaveOverlandFlowModel"


def test_KinWaveOF_input_var_names(kin_wave_of):
    assert kin_wave_of.input_var_names == (
        "topographic__elevation",
        "topographic__gradient",
    )


def test_KinWaveOF_output_var_names(kin_wave_of):
    assert kin_wave_of.output_var_names == (
        "surface_water__depth",
        "water__specific_discharge",
        "water__velocity",
    )


def test_KinWaveOF_var_units(kin_wave_of):
    assert set(kin_wave_of.input_var_names) | set(kin_wave_of.output_var_names) == set(
        dict(kin_wave_of.units).keys()
    )

    assert kin_wave_of.var_units("topographic__elevation") == "m"
    assert kin_wave_of.var_units("topographic__gradient") == "m/m"
    assert kin_wave_of.var_units("surface_water__depth") == "m"
    assert kin_wave_of.var_units("water__velocity") == "m/s"
    assert kin_wave_of.var_units("water__specific_discharge") == "m2/s"


def test_grid_shape(kin_wave_of):
    assert kin_wave_of.grid.number_of_node_rows == _SHAPE[0]
    assert kin_wave_of.grid.number_of_node_columns == _SHAPE[1]


def test_run_one_step():
    import numpy as np

    from landlab import RasterModelGrid
    from landlab.components.overland_flow import KinwaveOverlandFlowModel

    grid = RasterModelGrid((10, 10), xy_spacing=0.5)
    grid.add_zeros("topographic__elevation", at="node", dtype=float)
    grid.add_zeros("topographic__gradient", at="link")

    topo_arr = np.ones(100).reshape(10, 10)
    i = 0
    while i <= 9:
        topo_arr[:, i] = 5 + (0.002 * i)
        i += 1
    topo_arr = topo_arr.flatten()
    grid["node"]["topographic__elevation"] = topo_arr
    KinWaveOF = KinwaveOverlandFlowModel(
        grid, precip_rate=100.0, precip_duration=1.0, roughness=0.02
    )

    KinWaveOF.run_one_step(60)

    # I'll admit this is very non-robust. Solution roughly based on plot #9
    # from Heng et. al, (2009): "Modeling overland flow and soil eroion on
    # non uniform hillslopes: A finite volume scheme." They do not provide the
    # numerical solution but the plots match...
    max_h_mm = max(grid["node"]["surface_water__depth"]) * 1000.0
    np.testing.assert_almost_equal(max_h_mm, 1.66666666667)
