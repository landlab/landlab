import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import FlowDirectorMFD
from landlab.components.mass_wasting_runout import MassWastingRunout


@pytest.fixture
def example_square_mg():
    """sloped, convergent, irregular surface"""
    dem = np.array(
        [
            [10, 8, 4, 3, 4, 7.5, 10],
            [10, 9, 3.5, 4, 5, 8, 10],
            [10, 9, 6.5, 5, 6, 8, 10],
            [10, 9.5, 7, 6, 7, 9, 10],
            [10, 10, 9.5, 8, 9, 9.5, 10],
            [10, 10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10, 10],
        ]
    )
    dem = np.hstack(dem).astype(float)
    mg = RasterModelGrid((7, 7), 10)
    _ = mg.add_field("topographic__elevation", dem, at="node")
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    mg.at_node["node_id"] = np.hstack(mg.nodes)
    fd = FlowDirectorMFD(mg, diagonals=True, partition_method="slope")
    fd.run_one_step()
    nn = mg.number_of_nodes
    mg.at_node["mass__wasting_id"] = np.zeros(nn).astype(int)
    mg.at_node["mass__wasting_id"][np.array([38])] = 1
    depth = np.ones(nn) * 1
    mg.add_field("node", "soil__thickness", depth)
    np.random.seed(seed=7)
    mg.at_node["particle__diameter"] = np.random.uniform(0.05, 0.25, nn)
    mg.at_node["organic__content"] = np.random.uniform(0.01, 0.10, nn)
    return mg


@pytest.fixture
def example_square_MWRu(example_square_mg):
    slpc = [0.03]
    qsc = 0.01
    k = 0.02
    mofd = 1
    tracked_attributes = ["particle__diameter", "organic__content"]
    example_square_MWRu = MassWastingRunout(
        example_square_mg,
        critical_slope=slpc,
        threshold_flux=qsc,
        erosion_coefficient=k,
        max_flow_depth_observed_in_field=mofd,
        save=True,
        tracked_attributes=tracked_attributes,
    )
    return example_square_MWRu


@pytest.fixture
def example_flat_mg():
    "small, flat surface"
    dem = np.ones(25)
    mg = RasterModelGrid((5, 5), 10)
    _ = mg.add_field("topographic__elevation", dem, at="node")
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    mg.set_watershed_boundary_condition_outlet_id(3, dem)
    mg.at_node["node_id"] = np.hstack(mg.nodes)
    nn = mg.number_of_nodes
    mg.at_node["mass__wasting_id"] = np.zeros(mg.number_of_nodes).astype(int)
    depth = np.ones(nn) * 1
    mg.add_field("node", "soil__thickness", depth)
    return mg


@pytest.fixture
def example_bumpy_mg():
    """sloped, irregular surface"""
    dem = np.ones(25)
    mg = RasterModelGrid((5, 5), 10)
    _ = mg.add_field("topographic__elevation", dem, at="node")
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    mg.set_watershed_boundary_condition_outlet_id(3, dem)
    mg.at_node["node_id"] = np.hstack(mg.nodes)
    mg.at_node["topographic__elevation"][np.array([6, 7, 8, 11, 13, 16, 17, 18])] = (
        np.array([3, 2, 5, 5, 7, 9, 8, 11])
    )
    nn = mg.number_of_nodes
    mg.at_node["mass__wasting_id"] = np.zeros(mg.number_of_nodes).astype(int)
    depth = np.ones(nn) * 1
    mg.add_field("node", "soil__thickness", depth)
    return mg


@pytest.fixture
def example_pile_MWRu():
    """sets up a MWR for modeling runout of a pile of debris"""
    # odd number
    r = 51
    c = r
    dxdy = 5
    ls_h = 5
    w = 5
    slpc = [0.03]
    qsc = 0.01
    k = 0.02
    eta = 0.2  # exponent
    Dp = 0.2  # particle diameter
    ls_h = 5  # landslide (pile) thickness
    hs = 1  # soil thickness
    # create model grid
    mg = RasterModelGrid((r, c), dxdy)
    mg.at_node["node_id"] = np.hstack(mg.nodes)
    # set boundary conditions
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    # soil thickness
    thickness = np.ones(mg.number_of_nodes) * hs
    mg.add_field("node", "soil__thickness", thickness)
    mg.add_field("topographic__elevation", np.ones(r * c) * 1, at="node")
    # set particle diameter
    mg.at_node["particle__diameter"] = np.ones(len(mg.node_x)) * Dp
    # run flow director, add slope and receiving node fields
    fd = FlowDirectorMFD(mg, diagonals=True, partition_method="square_root_of_slope")
    fd.run_one_step()

    # create pile
    # find central point in domain
    x = mg.node_x.max() / 2
    y = mg.node_y.max() / 2
    # find all nodes with radius of central point
    dn = ((mg.node_x - x) ** 2 + (mg.node_y - y) ** 2) ** 0.5
    pile_nodes = np.hstack(mg.nodes)[dn < w * mg.dx]
    mg.at_node["mass__wasting_id"] = np.zeros(mg.number_of_nodes).astype(int)
    mg.at_node["mass__wasting_id"][pile_nodes] = 1
    # set thickness of landslide
    mg.at_node["soil__thickness"][pile_nodes] = ls_h
    mg.at_node["topographic__elevation"][pile_nodes] = mg.at_node[
        "topographic__elevation"
    ][pile_nodes] + (ls_h - hs)

    # profile nodes and initial topography
    pf = mg.nodes[int((r - 1) / 2), :]
    # set up MWR
    MWRu = MassWastingRunout(
        mg,
        critical_slope=slpc,
        threshold_flux=qsc,
        erosion_coefficient=k,
        erosion_exponent=eta,
        effective_qsi=False,
        save=True,
        grain_shear=False,
        settle_deposit=True,
    )
    MWRu.r = r
    MWRu.pf = pf
    return MWRu


@pytest.fixture
def example_flume_MWRu():
    dxdy = 10
    rows = 15  #
    columns = 15  # must be odd number
    ls_width = 3  # number of cells wide, must be odd number
    ls_length = 5  # number of cells long
    slope_above_break = 0.6
    slope_below_break = 0.001
    slope_break = 0.15  # ratio of flume 2 to place the slope break
    soil_thickness = 2
    Dp = 0.2  # particle diameter
    qsi_max = 5
    ls_h = 5

    # run parameters
    qsc = 0.01  # pick qsc
    Sc = 0.05  # critical slope
    k = 0.00328  # erosion coefficient

    mg, lsn, pf, cc = flume_maker(
        rows=rows,
        columns=columns,
        slope_above_break=slope_above_break,
        slope_below_break=slope_below_break,
        slope_break=slope_break,
        ls_width=ls_width,
        ls_length=ls_length,
        dxdy=dxdy,
        double_flume=True,
    )
    # set boundary conditions
    mg.set_closed_boundaries_at_grid_edges(
        True, True, True, True
    )  # close all boundaries
    dem = mg.at_node["topographic__elevation"]
    mg.set_watershed_boundary_condition_outlet_id(cc, dem)

    # add flow direction field
    fd = FlowDirectorMFD(mg, diagonals=True, partition_method="square_root_of_slope")
    fd.run_one_step()

    # soil thickness
    thickness = np.ones(mg.number_of_nodes) * soil_thickness
    mg.add_field("node", "soil__thickness", thickness)

    # define landslide location
    mg.at_node["mass__wasting_id"] = np.zeros(mg.number_of_nodes).astype(int)
    mg.at_node["mass__wasting_id"][lsn] = 1

    # define landslide thickness
    mg.at_node["soil__thickness"][lsn] = ls_h
    # set particle diameter
    mg.at_node["particle__diameter"] = np.ones(len(mg.node_x)) * Dp

    example_flume_MWRu = MassWastingRunout(
        mg,
        critical_slope=[Sc],
        threshold_flux=qsc,
        erosion_coefficient=k,
        max_flow_depth_observed_in_field=qsi_max,
        tracked_attributes=["particle__diameter"],
        grain_shear=True,
        save=True,
    )
    example_flume_MWRu.pf = pf
    return example_flume_MWRu


def flume_maker(
    rows=5,
    columns=3,
    slope_above_break=0.5,
    slope_below_break=0.05,
    slope_break=0.7,
    ls_width=1,
    ls_length=1,
    dxdy=10,
    double_flume=False,
):
    """
    Parameters
    ----------
    rows : integer
        number of rows in model domain. First and last row are used as boundaries
    columns : integer
        number of columns in domain. First and last column are used as boundaries
    slope_above_break : float
        slope of flume above break [m/m]
    slope_below_break : float
        slope of flume below break [m/m]
    slope_break : float
        ratio of length the slope break is placed, measured from the outlet
        slope break will be placed at closest node to slope_break value
    ls_width : odd value, integer
        width of landslide in number of cells, , must be <= than rows-2
    ls_length : integer
        length of landslide in number of cells, must be <= than rows-2
    dxdy : float
        side length of sqaure cell, [m]. The default is 10.
    double_flume : bool
        False: makes just one flume; True: makes two flumes and puts attaches the
        lower end of one to the upper end of the other to make one long flume with
        a slope break in the middle of the flume

    Returns
    -------
    mg : raster model grid
        includes the field topographic__elevation
    lsn : np array
        0-d array of node id's that are the landslide
    pf : np array
        0-d array of node id's along the center of the flume. Used for profile
        plots.
    cc : np array
        0-d array of landslide node column id's'
    """

    def single_flume(slope_break):
        r = rows
        c = columns
        sbr_r = slope_break
        if ls_width == 1:
            cc = int(c / 2) + 1
        elif (ls_width <= c) and (ls_width % 2 == 1):
            cc = []
            for i in range(ls_width):
                dif = -((ls_width) % 2) + i
                cc.append(int(c / 2) + dif + 1)
            cc = np.array(cc)
        mg = RasterModelGrid((r, c + 2), dxdy)
        ycol = np.reshape(mg.node_y, (r, c + 2))[:, 0]
        yn = np.arange(r)
        yeL = []
        sb = mg.node_y.max() * sbr_r
        sbr = yn[ycol >= sb].min()
        sb = ycol[sbr]
        for y in ycol:
            if y < sb:
                ye = y * slope_below_break
            else:
                ye = (y - sb) * slope_above_break + sb * slope_below_break
            yeL.append(ye)
        dem = np.array(
            [
                yeL,
            ]
            * c
        ).transpose()
        wall = np.reshape(dem[:, 0], (len(dem[:, 0]), 1)) + 2 * dxdy
        dem = np.concatenate((wall, dem, wall), axis=1)
        dem = np.hstack(dem).astype(float)
        _ = mg.add_field("topographic__elevation", dem, at="node")
        # profile nodes
        if ls_width == 1:
            pf = mg.nodes[:, cc]
        elif ls_width > 1:
            pf = mg.nodes[:, int(c / 2) + 1]
        # landslide nodes
        lsn = mg.nodes[-(ls_length + 1) : -1, cc]
        return mg, lsn, pf, cc

    if double_flume:
        mg1, lsn1, pf1, cc1 = single_flume(slope_break)
        mg2, lsn2, pf2, cc2 = single_flume(min(slope_break * 2, 0.8))
        mg = RasterModelGrid((rows * 2, columns + 2), dxdy)
        # translate mg1 vertically to the max height of mg2 minus the wall
        # height plus the flume slope * dxdy
        t1 = (
            mg1.at_node["topographic__elevation"]
            + mg2.at_node["topographic__elevation"].max()
            - dxdy * 2
            + dxdy * slope_above_break
        )
        t2 = mg2.at_node["topographic__elevation"]
        topo = np.concatenate((t2, t1))
        _ = mg.add_field("topographic__elevation", topo, at="node")
        nn = len(mg2.node_x)
        pf = np.concatenate((pf2, np.array(pf1 + nn).astype(int)))
        cc = cc1
        lsn = lsn1 + nn
    else:
        mg, lsn, pf, cc = single_flume()
    return mg, lsn, pf, cc
