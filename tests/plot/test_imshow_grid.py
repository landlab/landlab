import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.backends.backend_pdf import PdfPages

import landlab


@pytest.mark.slow
def test_imshow_grid():
    rmg = landlab.RasterModelGrid((4, 5))

    pp = PdfPages("test.pdf")

    values = np.arange(rmg.number_of_nodes)
    landlab.plot.imshow_grid(rmg, values, values_at="node", limits=(0, 20))
    pp.savefig()

    plt.clf()
    rmg.status_at_node[7] = rmg.BC_NODE_IS_CLOSED
    values = np.arange(rmg.number_of_cells)
    landlab.plot.imshow_grid(rmg, values, values_at="cell", symmetric_cbar=True)
    pp.savefig()
    pp.close()


def test_imshow_grid_input():
    rmg = landlab.RasterModelGrid((4, 5))
    values = np.arange(rmg.number_of_nodes - 1)
    with pytest.raises(ValueError):
        _ = landlab.plot.imshow_grid(rmg, values, values_at="node", limits=(0, 20))


def test_imshowhs_grid_input():
    rmg = landlab.RasterModelGrid((4, 5))
    values = np.arange(rmg.number_of_nodes - 1)
    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(rmg, values, values_at="node", limits=(0, 20))


def test_imshowhs_grid_input_Layer1():
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    values1 = np.arange(mg.number_of_nodes - 1)

    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=values1,
            plot_type="Drape1",
            var_name="Soil",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            thicks_km=False,
            limits=(0, 2),
        )


def test_imshowhs_grid_input_Layer2():
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    values1 = np.arange(mg.number_of_nodes)
    values2 = np.arange(mg.number_of_nodes - 1)

    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=values1,
            drape2=values2,
            plot_type="Drape2",
            var_name="Soil",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            thicks_km=False,
            limits=(0, 2),
        )


def test_imshowhs_grid_1():
    """
    Show DEM draped over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        symmetric_cbar=True,
        limits=(0, 10),
    )


# %%
def test_imshowhs_grid_2():
    """
    Show DEM draped over the shaded topographic relief with exaggeration
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        vertical_exa=2,
        thicks_km=True,
        symmetric_cbar=True,
        vmin=0,
        vmax=10,
    )


def test_imshowhs_grid_3():
    """
    Show Hillshade
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        plot_type="Hillshade",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        plt_contour=True,
        vmax=10,
        vmin=0,
    )


def test_imshowhs_grid_4a():
    """
    Show Drape1 draped over the shaded topographic relief with exaggeration
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    # Show Soil thickness draped over the shaded topographic relief
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["soil__depth"],
        plot_type="Drape1",
        var_name="Soil",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
    )


def test_imshowhs_grid_4b():
    """
    Show Drape1 draped over the shaded topographic relief with exaggeration
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    # Show Soil thickness draped over the shaded topographic relief
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["soil__depth"],
        plot_type="Drape1",
        var_name="Soil",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        vmin=0,
        vmax=2,
        plt_contour=True,
    )


def test_imshowhs_grid_4c():
    """
    Show Drape1 draped over the shaded topographic relief with exaggeration
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    # Show Soil thickness draped over the shaded topographic relief
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["soil__depth"],
        plot_type="Drape1",
        var_name="Soil",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        symmetric_cbar=True,
    )


# %%


def test_imshowhs_grid_5():
    """
    Show Drape1 draped over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    _ = mg.add_zeros("Layer_1", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        plot_type="Drape1",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
    )


def test_imshowhs_grid_6a():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    L1 = mg.add_zeros("Layer_1", at="node")
    L2 = mg.add_zeros("Layer_2", at="node")
    L1[:] += 10
    L2[:] += 100
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 200),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
    )


def test_imshowhs_grid_6b():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <10
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    L1 = mg.add_zeros("Layer_1", at="node")
    L2 = mg.add_zeros("Layer_2", at="node")
    L1[:] += 10
    L2[:] += 100
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        vmin=0,
        vmax=9,
    )


def test_imshowhs_grid_6c():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    L1 = mg.add_zeros("Layer_1", at="node")
    L2 = mg.add_zeros("Layer_2", at="node")
    L1[:] += 10
    L2[:] += 100
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        vmin=0,
        vmax=99,
    )


def test_imshowhs_grid_6d():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    L1 = mg.add_zeros("Layer_1", at="node")
    L2 = mg.add_zeros("Layer_2", at="node")
    L1[:] += 10
    L2[:] += 100
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        vmin=0,
        vmax=999,
    )
    # %%


def test_imshowhs_grid_6e():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax <100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    L1 = mg.add_zeros("Layer_1", at="node")
    L2 = mg.add_zeros("Layer_2", at="node")
    L1[:] += 10
    L2[:] += 100
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        add_double_colorbar=True,
        vmin=0,
        vmax=99999,
    )


# %%


def test_imshowhs_grid_7():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    _ = mg.add_zeros("Layer_1", at="node")
    _ = mg.add_zeros("Layer_2", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1="topographic__elevation",
        drape2="soil__depth",
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        thres_drape2=1,
        cmap2=None,
        add_double_colorbar=True,
    )


# %%
def test_imshowhs_grid_8():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax >10<100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    L1 = mg.add_zeros("Layer_1", at="node")
    L2 = mg.add_zeros("Layer_2", at="node")
    L1[:] += 10
    L2[:] += 100
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1="topographic__elevation",
        drape2="soil__depth",
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        thres_drape2=1,
        cmap2=None,
        add_double_colorbar=True,
        vmin=0,
        vmax=99,
    )


# %%
def test_imshowhs_grid_9():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief, vmax>100
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    _ = mg.add_zeros("Layer_1", at="node")
    _ = mg.add_zeros("Layer_2", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1="topographic__elevation",
        drape2="soil__depth",
        plot_type="Drape2",
        var_name="Layer 1",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
        color_for_closed="red",
        thres_drape2=1,
        cmap2=None,
        add_double_colorbar=True,
        vmin=0,
        vmax=99999,
    )
    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
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
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    _ = mg.add_zeros("Layer_1", at="node")
    _ = mg.add_zeros("Layer_2", at="node")
    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=mg.at_node["Layer_1"],
            plot_type="Drape2",
            var_name="Layer 1",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            thicks_km=False,
            limits=(0, 2),
            colorbar_label_y=-55,
            add_label_bbox=True,
            thres_drape1=0.001,
        )


def test_imshowhs_grid_11():
    """
    Show Layer 1 and Layer 2 over the shaded topographic relief
    """
    # %%
    mg = landlab.RasterModelGrid((4, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    _ = mg.add_zeros("Layer_1", at="node")
    _ = mg.add_zeros("Layer_2", at="node")
    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            plot_type="Drape1",
            var_name="Layer 1",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            thicks_km=False,
            limits=(0, 2),
            colorbar_label_y=-55,
            add_label_bbox=True,
            thres_drape1=0.001,
        )


def test_hex_mfd():
    """
    Currently no support for hex
    """
    # %%
    mg = landlab.HexModelGrid((5, 3))
    _ = mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(NotImplementedError):
        _ = landlab.plot.imshowhs_grid(mg, "topographic__elevation")


# %%
def test_at_cell():
    """
    Currently no support for at cell
    """
    # %%
    mg = landlab.HexModelGrid((5, 3))
    _ = mg.add_field("topographic__elevation", np.zeros((7,)), at="cell")
    with pytest.raises(NotImplementedError):
        _ = landlab.plot.imshowhs_grid(mg, "topographic__elevation", at="cell")


# %%
def test_at_other():
    """
    Currently no support for non at node valley locations
    """
    # %%
    mg = landlab.HexModelGrid((5, 3))
    _ = mg.add_field("topographic__elevation", np.zeros((24,)), at="corner")
    with pytest.raises(TypeError):
        _ = landlab.plot.imshowhs_grid(mg, "topographic__elevation", at="corner")
