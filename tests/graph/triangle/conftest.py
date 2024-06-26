import json

import numpy as np
import pytest


@pytest.fixture(scope="session")
def geojson_concave_polygon(tmp_path_factory):
    path_to_file = tmp_path_factory.mktemp("geojson") / "polygon_concave.geojson"
    with open(path_to_file, "w") as fp:
        json.dump(
            {
                "type": "Feature",
                "properties": {},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [[0, 0], [10, 0], [10, 10], [5, 15], [0, 10], [0, 0]],
                        [[2, 2], [8, 2], [8, 8], [2, 8], [2, 2]],
                    ],
                },
            },
            fp,
        )

    return path_to_file


@pytest.fixture(scope="session")
def geojson_circular_polygon(tmp_path_factory):
    path_to_file = tmp_path_factory.mktemp("geojson") / "polygon_circular.geojson"

    x = np.cos(np.linspace(0.0, 2.0 * np.pi, num=64, endpoint=False))
    y = np.sin(np.linspace(0.0, 2.0 * np.pi, num=64, endpoint=False))

    with open(path_to_file, "w") as fp:
        json.dump(
            {
                "type": "Feature",
                "properties": {},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [np.c_[x, y].tolist() + [[1.0, 0.0]]],
                },
            },
            fp,
        )

    return path_to_file


@pytest.fixture(scope="session")
def geojson_interior_rings(tmp_path_factory):
    path_to_file = (
        tmp_path_factory.mktemp("geojson") / "polygon_two_interior_rings.geojson"
    )
    with open(path_to_file, "w") as fp:
        json.dump(
            {
                "type": "Feature",
                "properties": {},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
                        [[3, 3], [4, 3], [4, 4], [3, 4], [3, 3]],
                        [[6, 6], [7, 6], [7, 7], [6, 7], [6, 6]],
                    ],
                },
            },
            fp,
        )

    return path_to_file
