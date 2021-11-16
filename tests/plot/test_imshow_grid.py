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


def test_imshowhs_grid():  # Show DEM draped over the shaded topographic relief
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
    )

    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        plot_type="Hillshade",
        var_name="Topo",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
    )

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
    # Show LAyer 1 over the shaded topographic relief
    _ = mg.add_zeros("Layer_1", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        plot_type="Drape1",
        var_name="LS \n erosion",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
    )

    # Show Layer 1 and Layer 2 over the shaded topographic relief
    _ = mg.add_zeros("Layer_2", at="node")
    _ = landlab.plot.imshowhs_grid(
        mg,
        "topographic__elevation",
        drape1=mg.at_node["Layer_1"],
        drape2=mg.at_node["Layer_2"],
        plot_type="Drape2",
        var_name="LS \n erosion",
        var_units=r"m",
        grid_units=("m", "m"),
        cmap="terrain",
        thicks_km=False,
        limits=(0, 2),
        colorbar_label_y=-55,
        add_label_bbox=True,
        thres_drape1=0.001,
    )

    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            plot_type="Oops",
        )

    with pytest.raises(ValueError):
        _ = landlab.plot.imshowhs_grid(
            mg,
            "topographic__elevation",
            drape1=mg.at_node["Layer_1"],
            plot_type="Drape2",
            var_name="LS \n erosion",
            var_units=r"m",
            grid_units=("m", "m"),
            cmap="terrain",
            thicks_km=False,
            limits=(0, 2),
            colorbar_label_y=-55,
            add_label_bbox=True,
            thres_drape1=0.001,
        )
