import json
import pytest


@pytest.fixture(scope="session")
def concave_polygon(tmp_path_factory):
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
