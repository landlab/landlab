import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.backends.backend_pdf import PdfPages

import landlab
from landlab.plot.imshow import _guess_location_from_name
from landlab.plot.imshow import _guess_location_from_size


@pytest.mark.slow
def test_imshow_grid():
    rmg = landlab.RasterModelGrid((4, 5))

    pp = PdfPages("test.pdf")

    values = np.arange(rmg.number_of_nodes)
    landlab.plot.imshow_grid(rmg, values, at="node", limits=(0, 20))
    pp.savefig()

    plt.clf()
    rmg.status_at_node[7] = rmg.BC_NODE_IS_CLOSED
    values = np.arange(rmg.number_of_cells)
    landlab.plot.imshow_grid(rmg, values, at="cell", symmetric_cbar=True)
    pp.savefig()
    pp.close()


def test_imshow_grid_input():
    rmg = landlab.RasterModelGrid((4, 5))
    values = np.arange(rmg.number_of_nodes - 1)
    with pytest.raises(ValueError):
        landlab.plot.imshow_grid(rmg, values, at="node", limits=(0, 20))


def test_imshowhs_grid_array_wrong_size():
    rmg = landlab.RasterModelGrid((4, 5))
    values = np.arange(rmg.number_of_nodes - 1)
    with pytest.raises(ValueError):
        landlab.plot.imshowhs_grid(rmg, values, at="node", limits=(0, 20))


def test_imshowhs_grid_array_wrong_size_with_field_name():
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    values1 = np.arange(mg.number_of_nodes - 1)

    with pytest.raises(ValueError):
        landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=values1,
            plot_type="Drape1",
            var_name="Soil",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            ticks_km=False,
            limits=(0, 2),
        )


def test_imshowhs_grid_array_size_mismatch():
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    values1 = np.arange(mg.number_of_nodes)
    values2 = np.arange(mg.number_of_nodes - 1)

    with pytest.raises(ValueError):
        landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=values1,
            drape2=values2,
            plot_type="Drape2",
            var_name="Soil",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            ticks_km=False,
            limits=(0, 2),
        )


def test_imshowhs_grid_1():
    """Show DEM draped over the shaded topographic relief"""
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        symmetric_cbar=True,
        limits=(0, 10),
    )


def test_imshowhs_grid_2():
    """Show DEM draped over the shaded topographic relief with exaggeration"""
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        vertical_exa=2,
        ticks_km=True,
        symmetric_cbar=True,
        vmin=0,
        vmax=10,
    )


def test_imshowhs_grid_3():
    """Show Hillshade"""
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        plot_type="Hillshade",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        plt_contour=True,
        vmax=10,
        vmin=0,
    )


def test_imshowhs_grid_4a():
    """
    Show Drape1 draped over the shaded topographic relief with exaggeration
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    # Show Soil thickness draped over the shaded topographic relief
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["soil__depth"],
        plot_type="Drape1",
        var_name="Soil",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        limits=(0, 2),
    )


def test_imshowhs_grid_4b():
    """
    Show Drape1 draped over the shaded topographic relief with exaggeration
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    # Show Soil thickness draped over the shaded topographic relief
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["soil__depth"],
        plot_type="Drape1",
        var_name="Soil",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        vmin=0,
        vmax=2,
        plt_contour=True,
    )


def test_imshowhs_grid_4c():
    """
    Show Drape1 draped over the shaded topographic relief with exaggeration
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    # Show Soil thickness draped over the shaded topographic relief
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["soil__depth"],
        plot_type="Drape1",
        var_name="Soil",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        symmetric_cbar=True,
    )


def test_imshowhs_grid_5():
    """
    Show Drape1 draped over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("Layer_1", at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        plot_type="Drape1",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
    )


def test_imshowhs_grid_6a():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_full("Layer_1", 10.0, at="node")
    mg.add_full("Layer_2", 100.0, at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        limits=(0, 200),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
    )


def test_imshowhs_grid_6b():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <10
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_full("Layer_1", 10.0, at="node")
    mg.add_full("Layer_2", 100.0, at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
        vmin=0,
        vmax=9,
    )


def test_imshowhs_grid_6c():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_full("Layer_1", 10.0, at="node")
    mg.add_full("Layer_2", 100.0, at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
        vmin=0,
        vmax=99,
    )


def test_imshowhs_grid_6d():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <100
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_full("Layer_1", 10.0, at="node")
    mg.add_full("Layer_2", 100.0, at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
        vmin=0,
        vmax=999,
    )


def test_imshowhs_grid_6e():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <100
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_full("Layer_1", 10.0, at="node")
    mg.add_full("Layer_2", 100.0, at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
        add_double_colorbar=True,
        vmin=0,
        vmax=99999,
    )


def test_imshowhs_grid_7():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("Layer_1", at="node")
    mg.add_zeros("Layer_2", at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1="topographic__elevation",
        drape2="soil__depth",
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
        thres_drape2=1,
        cmap2=None,
        add_double_colorbar=True,
    )


def test_imshowhs_grid_8():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax >10<100
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_full("Layer_1", 10.0, at="node")
    mg.add_full("Layer_2", 100.0, at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1="topographic__elevation",
        drape2="soil__depth",
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="white",
        thres_drape2=1,
        cmap2=None,
        add_double_colorbar=True,
        vmin=0,
        vmax=99,
    )


def test_imshowhs_grid_9():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax>100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("Layer_1", at="node")
    mg.add_zeros("Layer_2", at="node")
    landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1="topographic__elevation",
        drape2="soil__depth",
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        ticks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed=None,
        thres_drape2=1,
        cmap2=None,
        add_double_colorbar=True,
        vmin=0,
        vmax=99999,
    )
    with pytest.raises(ValueError):
        landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            plot_type="Oops",
        )


def test_imshowhs_grid_10():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("Layer_1", at="node")
    mg.add_zeros("Layer_2", at="node")
    with pytest.raises(ValueError):
        landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=mg.at_node["Layer_1"],
            plot_type="Drape2",
            var_name="Layer 1",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            ticks_km=False,
            limits=(0, 2),
            colorbar_label_y=-55,
            add_label_bbox=True,
            thres_drape1=0.001,
        )


def test_imshowhs_grid_11():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("Layer_1", at="node")
    mg.add_zeros("Layer_2", at="node")
    with pytest.raises(ValueError):
        landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            plot_type="Drape1",
            var_name="Layer 1",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            ticks_km=False,
            limits=(0, 2),
            colorbar_label_y=-55,
            add_label_bbox=True,
            thres_drape1=0.001,
        )


@pytest.mark.parametrize("valid_units", [None, "foo"])
def test_imshowhs_var_units(valid_units):
    """units should be a string or None"""
    mg = landlab.RasterModelGrid((4, 5))
    mg.add_zeros("topographic__elevation", at="node")
    landlab.plot.imshowhs_grid(mg, "topographic__elevation", var_units=valid_units)


def test_hex_grid_not_allowed():
    """Currently no support for hex"""
    mg = landlab.HexModelGrid((5, 3))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(NotImplementedError):
        landlab.plot.imshowhs_grid(mg, "topographic__elevation")


@pytest.mark.parametrize("at", ["link", "patch", "corner", "face", "cell"])
def test_at_anything_but_node(at):
    """Currently no support for anything but node."""
    mg = landlab.RasterModelGrid((5, 3))
    mg.add_empty("topographic__elevation", at=at)
    with pytest.raises(NotImplementedError):
        landlab.plot.imshowhs_grid(mg, "topographic__elevation", at=at)


@pytest.mark.parametrize("at", ["node", "cell"])
def test_imshow_grid_guess_from_name(at):
    grid = landlab.RasterModelGrid((3, 4))
    grid.add_zeros("z", at=at)
    landlab.plot.imshow_grid(grid, "z")


@pytest.mark.parametrize("at", ["node", "cell"])
def test_imshow_grid_guess_from_size(at):
    grid = landlab.RasterModelGrid((3, 4))
    values = grid.zeros(at=at)
    landlab.plot.imshow_grid(grid, values)


def test_imshow_grid_unknown_location():
    expected = "unable to determine location of values, use 'at' keyword"
    grid = landlab.RasterModelGrid((5, 3))
    values = np.empty(grid.number_of_links + 1)
    with pytest.raises(TypeError, match=expected):
        landlab.plot.imshow_grid(grid, values)
    with pytest.raises(TypeError, match=expected):
        landlab.plot.imshow_grid(grid, "foo")


@pytest.mark.parametrize("at", ["link", "patch", "corner", "face"])
def test_imshow_grid_unsupported_location(at):
    expected = (
        "value location, \\'[a-z]+\\', is not supported \\(must be "
        "one of 'node', 'cell'\\)"
    )
    grid = landlab.RasterModelGrid((5, 3))
    grid.add_zeros("z", at=at)
    with pytest.raises(TypeError, match=expected):
        landlab.plot.imshow_grid(grid, "z")
    with pytest.raises(TypeError, match=expected):
        landlab.plot.imshow_grid(grid, grid[at]["z"])


def test_values_at_is_deprecated():
    grid = landlab.RasterModelGrid((5, 3))
    grid.add_zeros("topographic__elevation", at="node")
    with pytest.deprecated_call(
        match="the 'values_at' keyword is deprecated, use the 'at' keyword instead"
    ):
        landlab.plot.imshow_grid(grid, "topographic__elevation", values_at="node")


@pytest.mark.parametrize("at", ["node", "link", "patch", "corner", "face", "cell"])
def test_guess_location(at):
    grid = landlab.RasterModelGrid((3, 4))
    values = grid.add_zeros("z", at=at)

    assert _guess_location_from_name(grid, "z") == at
    assert _guess_location_from_name(grid, "foo") is None

    guess = _guess_location_from_size(grid, np.empty_like(values))
    assert grid.number_of_elements(guess) == values.size
    assert _guess_location_from_size(grid, np.empty(grid.number_of_links + 1)) is None


@pytest.mark.parametrize("actual", ["node", "cell"])
@pytest.mark.parametrize("at", ["link", "patch", "corner", "face"])
def test_location_node_cell_first(actual, at):
    grid = landlab.RasterModelGrid((3, 4))
    values = grid.add_zeros("z", at=actual)
    grid.add_zeros("z", at=at)

    assert _guess_location_from_name(grid, "z") == actual
    assert _guess_location_from_size(grid, np.empty_like(values)) == actual
