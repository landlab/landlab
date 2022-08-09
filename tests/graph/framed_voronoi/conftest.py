from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph


def pytest_generate_tests(metafunc):
    if "layout_graph" in metafunc.fixturenames:
        metafunc.parametrize(
            "layout_graph",
            [HorizontalRectVoronoiGraph],
        )
