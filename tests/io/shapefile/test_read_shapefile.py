from io import BytesIO

import numpy as np
import shapefile
from numpy.testing import assert_array_equal
from pytest import approx, raises
from shapefile import ShapefileException

from landlab import ExampleData
from landlab.io.shapefile import read_shapefile


def test_read_methow(tmpdir):
    # test of the big methow network.
    shp_file = "Methow_Network.shp"
    with tmpdir.as_cwd():
        ExampleData("io/shapefile", case="methow").fetch()

        grid = read_shapefile(shp_file)
        assert grid.number_of_nodes == 721
        assert grid.number_of_links == 720

        assert "x_of_polyline" in grid.at_link
        assert "y_of_polyline" in grid.at_link

        fields = [
            "GridID",
            "Length_m",
            "ToLink",
            "usarea_km2",
            "uselev_m",
            "dselev_m",
            "Slope",
        ]
        for field in fields:
            assert field in grid.at_link

        # test for node location.
        assert grid.x_of_node[0] == approx(-1672349.0889982011)
        assert grid.y_of_node[0] == approx(1160800.240247)

        # test for node location.
        assert grid.x_of_node[25] == approx(-1677950.8375069483)
        assert grid.y_of_node[25] == approx(1173162.0319129999)

        # test for field order.
        # link 15
        assert grid.at_link["GridID"][15] == 16
        assert grid.at_link["Length_m"][15] == approx(1593.206627)
        assert grid.at_link["ToLink"][15] == 9
        assert grid.at_link["usarea_km2"][15] == 4573.3896
        assert grid.at_link["uselev_m"][15] == 305.88
        assert grid.at_link["dselev_m"][15] == 295.37
        assert grid.at_link["Slope"][15] == 0.006597

        # link 27
        assert grid.at_link["GridID"][27] == 41
        assert grid.at_link["Length_m"][27] == approx(1576.988871)
        assert grid.at_link["ToLink"][27] == 36
        assert grid.at_link["usarea_km2"][27] == 78.0597
        assert grid.at_link["uselev_m"][27] == 524.54
        assert grid.at_link["dselev_m"][27] == 399.81
        assert grid.at_link["Slope"][27] == 0.079094


def test_read_methow_subbasin(tmpdir):
    # test of the small methow network with
    shp_file = "MethowSubBasin.shp"
    points_shapefile = "MethowSubBasin_Nodes_4.shp"
    with tmpdir.as_cwd():
        ExampleData("io/shapefile", case="methow").fetch()

        grid = read_shapefile(
            shp_file, points_shapefile=points_shapefile, threshold=0.01
        )
        assert grid.number_of_nodes == 30
        assert grid.number_of_links == 29

        assert "x_of_polyline" in grid.at_link
        assert "y_of_polyline" in grid.at_link

        # verify that fields are present.
        node_fields = [
            "GridID",
            "ToLink",
            "usarea_km2",
            "uselev_m",
            "dselev_m",
            "Elev_m",
        ]
        link_fields = [
            "GridID",
            "Length_m",
            "ToLink",
            "usarea_km2",
            "uselev_m",
            "dselev_m",
            "Slope",
        ]
        for field in link_fields:
            assert field in grid.at_link
        for field in node_fields:
            assert field in grid.at_node

        # test for node location.
        assert grid.x_of_node[0] == approx(728390.38284378243)
        assert grid.y_of_node[0] == approx(5368319.8330760002)

        assert grid.x_of_node[22] == approx(725029.95616002998)
        assert grid.y_of_node[22] == approx(5374213.8330760002)

        # verify that fields are mapped correctly. choose two links and two nodes
        # to test.

        # link 1
        assert grid.at_link["GridID"][0] == 267
        assert grid.at_link["Length_m"][0] == approx(1799.167)
        assert grid.at_link["ToLink"][0] == 270
        assert grid.at_link["usarea_km2"][0] == 22.3047
        assert grid.at_link["uselev_m"][0] == 1305.7
        assert grid.at_link["dselev_m"][0] == 1232.77
        assert grid.at_link["Slope"][0] == 0.040535

        # link 16
        assert grid.at_link["GridID"][16] == 300
        assert grid.at_link["Length_m"][16] == approx(1487.342911)
        assert grid.at_link["ToLink"][16] == 263
        assert grid.at_link["usarea_km2"][16] == 75.1707
        assert grid.at_link["uselev_m"][16] == 1070.98
        assert grid.at_link["dselev_m"][16] == 986.74
        assert grid.at_link["Slope"][16] == 0.056638

        # node 1
        assert grid.at_node["GridID"][0] == 267
        assert grid.at_node["ToLink"][0] == 270.0
        assert grid.at_node["usarea_km2"][0] == 22.3047
        assert grid.at_node["uselev_m"][0] == 1305.7
        assert grid.at_node["dselev_m"][0] == 1232.77
        assert grid.at_node["Elev_m"][0] == 1304.24

        # node 29
        assert grid.at_node["GridID"][29] == 339
        assert grid.at_node["ToLink"][29] == 341.0
        assert grid.at_node["usarea_km2"][29] == 15.4314
        assert grid.at_node["uselev_m"][29] == 1534.42
        assert grid.at_node["dselev_m"][29] == 1438.68
        assert grid.at_node["Elev_m"][29] == 1535.58


def test_read_methow_subbasin_with_name_mapping_and_field_subsetting(tmpdir):
    # test of the small methow network with
    shp_file = "MethowSubBasin.shp"
    points_shapefile = "MethowSubBasin_Nodes_4.shp"

    with tmpdir.as_cwd():
        ExampleData("io/shapefile", case="methow").fetch()

        grid = read_shapefile(
            shp_file,
            points_shapefile=points_shapefile,
            node_fields=["usarea_km2", "ToLink", "Elev_m"],
            link_fields=["usarea_km2", "ToLink"],
            link_field_conversion={
                "usarea_km2": "drainage_area",
                "ToLink": "shapefile_to",
            },
            node_field_conversion={
                "usarea_km2": "drainage_area",
                "Elev_m": "topographic__elevation",
                "ToLink": "shapefile_to",
            },
            link_field_dtype={"ToLink": np.int},
            node_field_dtype={"ToLink": np.int},
            threshold=0.01,
        )

        assert grid.number_of_nodes == 30
        assert grid.number_of_links == 29

        assert "x_of_polyline" in grid.at_link
        assert "y_of_polyline" in grid.at_link

        # verify that fields are present and that some fields are not present
        node_fields = ["shapefile_to", "drainage_area", "topographic__elevation"]
        link_fields = [
            "shapefile_to",
            "drainage_area",
        ]

        not_node_fields = [
            "GridID",
            "ToLink",
            "usarea_km2",
            "uselev_m",
            "dselev_m",
            "Elev_m",
        ]
        not_link_fields = [
            "GridID",
            "Length_m",
            "ToLink",
            "usarea_km2",
            "uselev_m",
            "dselev_m",
            "Slope",
        ]

        for field in link_fields:
            assert field in grid.at_link
        for field in node_fields:
            assert field in grid.at_node
        for field in not_link_fields:
            assert field not in grid.at_link
        for field in not_node_fields:
            assert field not in grid.at_node

        # test for node location.
        assert grid.x_of_node[0] == approx(728390.38284378243)
        assert grid.y_of_node[0] == approx(5368319.8330760002)

        assert grid.x_of_node[22] == approx(725029.95616002998)
        assert grid.y_of_node[22] == approx(5374213.8330760002)

        # verify dtype changes.
        assert grid.at_link["shapefile_to"].dtype == np.int
        assert grid.at_node["shapefile_to"].dtype == np.int

        # verify that fields are mapped correctly. choose two links and two nodes
        # to test.

        # link 1
        assert grid.at_link["shapefile_to"][0] == 270
        assert grid.at_link["drainage_area"][0] == 22.3047

        # link 16
        assert grid.at_link["shapefile_to"][16] == 263
        assert grid.at_link["drainage_area"][16] == 75.1707

        # node 1
        assert grid.at_node["shapefile_to"][0] == 270
        assert grid.at_node["drainage_area"][0] == 22.3047
        assert grid.at_node["topographic__elevation"][0] == 1304.24

        # node 29
        assert grid.at_node["shapefile_to"][29] == 341
        assert grid.at_node["drainage_area"][29] == 15.4314
        assert grid.at_node["topographic__elevation"][29] == 1535.58


def test_bad_file(tmpdir):
    shp_file = "bad_file.shp"
    with raises(ShapefileException):
        read_shapefile(shp_file)


def test_points(datadir):
    shp_file = datadir / "points.shp"
    with raises(ValueError):
        read_shapefile(shp_file)


def test_multipart(datadir):
    shp_file = datadir / "multipartpolyline.shp"
    with raises(ValueError):
        read_shapefile(shp_file)


def test_bad_points():
    shp = BytesIO()
    shx = BytesIO()
    dbf = BytesIO()
    w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)
    w.shapeType = 3
    w.field("spam", "N")
    w.line([[[5, 5], [10, 10]]])
    w.record(37)
    w.line([[[5, 0], [5, 5]]])
    w.record(100)
    w.line([[[5, 5], [0, 10]]])
    w.record(239)
    w.close()

    # pass a line shapefile here insted.
    p_shp = BytesIO()
    p_shx = BytesIO()
    p_dbf = BytesIO()
    p_w = shapefile.Writer(shp=p_shp, shx=p_shx, dbf=p_dbf)
    w.shapeType = 3
    p_w.field("spam", "N")
    p_w.line([[[5, 5], [10, 10]]])
    p_w.record(37)
    p_w.line([[[5, 0], [5, 5]]])
    p_w.record(100)
    p_w.line([[[5, 5], [0, 10]]])
    p_w.record(239)
    p_w.close()

    with raises(ValueError):
        read_shapefile(shp, dbf=dbf, points_shapefile=p_shp, points_dbf=p_dbf)


def test_points_but_too_far():
    shp = BytesIO()
    shx = BytesIO()
    dbf = BytesIO()
    w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)
    w.shapeType = 3
    w.field("spam", "N")
    w.line([[[5, 5], [10, 10]]])
    w.record(37)
    w.line([[[5, 0], [5, 5]]])
    w.record(100)
    w.line([[[5, 5], [0, 10]]])
    w.record(239)
    w.close()

    # make a
    p_shp = BytesIO()
    p_shx = BytesIO()
    p_dbf = BytesIO()
    p_w = shapefile.Writer(shp=p_shp, shx=p_shx, dbf=p_dbf)
    p_w.shapeType = 1
    p_w.field("eggs", "N")
    p_w.point(5, 0)
    p_w.record(2)
    p_w.point(5, 5)
    p_w.record(4)
    p_w.point(0, 10)
    p_w.record(8)
    p_w.point(12, 10)
    p_w.record(6)
    p_w.close()

    with raises(ValueError):
        read_shapefile(shp, dbf=dbf, points_shapefile=p_shp, points_dbf=p_dbf)


def test_points_but_not_one_one():
    shp = BytesIO()
    shx = BytesIO()
    dbf = BytesIO()
    w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)
    w.shapeType = 3
    w.field("spam", "N")
    w.line([[[5, 5], [10, 10]]])
    w.record(37)
    w.line([[[5, 0], [5, 5]]])
    w.record(100)
    w.line([[[5, 5], [0, 10]]])
    w.record(239)
    w.close()

    # make a
    p_shp = BytesIO()
    p_shx = BytesIO()
    p_dbf = BytesIO()
    p_w = shapefile.Writer(shp=p_shp, shx=p_shx, dbf=p_dbf)
    p_w.shapeType = 1
    p_w.field("eggs", "N")
    p_w.point(5, 0)
    p_w.record(2)
    p_w.point(5, 5)
    p_w.record(4)
    p_w.point(0, 10)
    p_w.record(8)
    p_w.point(10, 10)
    p_w.record(6)
    p_w.point(10, 10)
    p_w.record(7)
    p_w.close()

    with raises(ValueError):
        read_shapefile(shp, dbf=dbf, points_shapefile=p_shp, points_dbf=p_dbf)


def test_points_but_one_missing():
    shp = BytesIO()
    shx = BytesIO()
    dbf = BytesIO()
    w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)
    w.shapeType = 3
    w.field("spam", "N")
    w.line([[[5, 5], [10, 10]]])
    w.record(37)
    w.line([[[5, 0], [5, 5]]])
    w.record(100)
    w.line([[[5, 5], [0, 10]]])
    w.record(239)
    w.close()

    # make a
    p_shp = BytesIO()
    p_shx = BytesIO()
    p_dbf = BytesIO()
    p_w = shapefile.Writer(shp=p_shp, shx=p_shx, dbf=p_dbf)
    p_w.shapeType = 1
    p_w.field("eggs", "N")
    p_w.point(5, 0)
    p_w.record(2)
    p_w.point(5, 5)
    p_w.record(4)
    p_w.point(0, 10)
    p_w.record(8)
    p_w.close()

    with raises(ValueError):
        read_shapefile(shp, dbf=dbf, points_shapefile=p_shp, points_dbf=p_dbf)


def test_simple_reorder():

    orders = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]

    lines = [
        [[[5, 5], [10, 10]]],
        [[[5, 0], [5, 5]]],
        [[[5, 5], [0, 10]]],
    ]

    records = [37, 100, 239]

    point_orders = [
        (0, 1, 2, 3),
        (0, 2, 3, 1),
        (0, 3, 2, 1),
        (0, 2, 3, 1),
        (1, 0, 2, 3),
        (1, 0, 3, 2),
        (1, 2, 3, 0),
        (1, 3, 2, 0),
    ]

    points = [(5, 0), (5, 5), (0, 10), (10, 10)]
    point_records = [2, 4, 8, 6]

    for order in orders:

        for p_order in point_orders:

            shp = BytesIO()
            shx = BytesIO()
            dbf = BytesIO()

            w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)

            w.shapeType = 3
            w.field("spam", "N")

            for o in order:
                w.line(lines[o])
                w.record(records[o])
            w.close()

            p_shp = BytesIO()
            p_shx = BytesIO()
            p_dbf = BytesIO()
            p_w = shapefile.Writer(shp=p_shp, shx=p_shx, dbf=p_dbf)
            p_w.shapeType = 1
            p_w.field("eggs", "N")
            for po in p_order:
                p_w.point(*points[po])
                p_w.record(point_records[po])
            p_w.close()

            grid = read_shapefile(
                shp, dbf=dbf, points_shapefile=p_shp, points_dbf=p_dbf
            )

            assert_array_equal(grid.nodes, np.array([0, 1, 2, 3]))
            assert_array_equal(grid.x_of_node, np.array([5.0, 5.0, 0.0, 10.0]))
            assert_array_equal(grid.y_of_node, np.array([0.0, 5.0, 10.0, 10.0]))
            assert_array_equal(grid.nodes_at_link, np.array([[0, 1], [2, 1], [1, 3]]))
            assert "spam" in grid.at_link
            assert_array_equal(grid.at_link["spam"], np.array([100, 239, 37]))

            del grid, w, shp, shx, dbf
