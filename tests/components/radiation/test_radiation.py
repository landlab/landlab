"""
Unit tests for landlab.components.radiation.radiation
"""

import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import Radiation


@pytest.mark.parametrize(
    "name",
    [
        "radiation__extraterrestrial_flux",
        "radiation__clearsky_flux",
        "radiation__incoming_shortwave_flux",
        "radiation__net_flux",
        "radiation__net_longwave_flux",
        "radiation__net_shortwave_flux",
        "radiation__ratio_to_flat_surface",
    ],
)
def test_fields_created_and_initialized(name):
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")
    Radiation(grid)

    assert np.allclose(grid.at_cell[name], 0.0)


@pytest.mark.parametrize(
    "name",
    [
        "radiation__extraterrestrial_flux",
        "radiation__clearsky_flux",
        "radiation__incoming_shortwave_flux",
        "radiation__net_flux",
        "radiation__net_longwave_flux",
        "radiation__net_shortwave_flux",
        "radiation__ratio_to_flat_surface",
    ],
)
def test_fields_updated(name):
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")
    radiation = Radiation(grid)
    radiation.update()

    assert np.all(grid.at_cell[name] > 0.0)


@pytest.mark.parametrize(
    "kwds",
    [
        {"latitude": -100},
        {"latitude": 90.1},
        {"min_daily_temp": 20.0, "max_daily_temp": 15.0},
        {"albedo": -1.0},
        {"albedo": 10.0},
    ],
)
def test_validators(kwds):
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")
    with pytest.raises(ValueError):
        Radiation(grid, **kwds)


def test_slope():
    grid = RasterModelGrid((5, 4), xy_spacing=10.0)
    grid.at_node["topographic__elevation"] = [
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0],
        [2.0, 2.0, 2.0, 2.0],
        [3.0, 4.0, 4.0, 3.0],
        [4.0, 4.0, 4.0, 4.0],
    ]
    radiation = Radiation(grid)
    radiation.update()

    assert np.all(grid.at_cell["radiation__ratio_to_flat_surface"] > 1.0)

    grid = RasterModelGrid((5, 4), xy_spacing=10.0)
    grid.at_node["topographic__elevation"] = [
        [4.0, 4.0, 4.0, 4.0],
        [3.0, 4.0, 4.0, 3.0],
        [2.0, 2.0, 2.0, 2.0],
        [1.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    radiation = Radiation(grid)
    radiation.update()

    assert np.all(grid.at_cell["radiation__ratio_to_flat_surface"] < 1.0)


def test_turbidity():
    grid = RasterModelGrid((5, 5), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")

    radiations = []
    for value in np.linspace(0.0, 1.0):
        radiation = Radiation(
            grid, latitude=0.0, clearsky_turbidity=value, opt_airmass=1.0
        )
        radiation.update()
        radiations.append(grid.at_cell["radiation__net_flux"].max())

    assert np.all(np.diff(radiations) < 0.0)


@pytest.mark.parametrize("time_0,time_1", [(0.0, 1.0), (2000.5, 2010.5), (-0.5, 0.5)])
def test_time(time_0, time_1):
    grid = RasterModelGrid((5, 5), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")
    radiation = Radiation(grid, latitude=0.0, current_time=time_0)
    radiation.update()

    radiation_at_time_0 = grid.at_cell["radiation__net_shortwave_flux"].copy()

    radiation = Radiation(grid, latitude=0.0, current_time=time_1)
    radiation.update()

    radiation_at_time_1 = grid.at_cell["radiation__net_shortwave_flux"]

    assert np.allclose(radiation_at_time_0, radiation_at_time_1)


TIMES = {
    "march_equinox": 1.39 * 365.0 / (2.0 * np.pi),
    "june_solstice": (np.pi / 2.0 + 1.39) * 365 / (2.0 * np.pi),
    "september_equinox": (np.pi + 1.39) * 365.0 / (2.0 * np.pi),
    "december_solstice": (3.0 * np.pi / 2.0 + 1.39) * 365.0 / (2.0 * np.pi),
}


@pytest.mark.parametrize(
    "latitude,solstice",
    [
        (70.0, "december"),
        (80.0, "december"),
        (-70.0, "june"),
        (-80.0, "june"),
    ],
)
def test_solstice(solstice, latitude):
    time = TIMES[f"{solstice}_solstice"] / 365.0

    grid = RasterModelGrid((5, 5), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")
    radiation = Radiation(grid, latitude=latitude, current_time=time)
    radiation.update()

    assert np.allclose(grid.at_cell["radiation__extraterrestrial_flux"], 0.0)


@pytest.mark.parametrize("time", ["march_equinox", "september_equinox"])
@pytest.mark.parametrize("latitude", [10.0, 30.0, 45.0, 60.0])
def test_equinox(time, latitude):
    """Check that radiation at the equinox is symmetric."""
    equinox = TIMES[time] / 365.0

    grid = RasterModelGrid((5, 5), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")

    radiation = Radiation(grid, latitude=latitude, current_time=equinox)
    radiation.update()

    radiation_at_north = grid.at_cell["radiation__ratio_to_flat_surface"].copy()

    radiation = Radiation(grid, latitude=-latitude, current_time=equinox)
    radiation.update()

    radiation_at_south = grid.at_cell["radiation__ratio_to_flat_surface"]

    assert np.allclose(radiation_at_north, radiation_at_south)


@pytest.mark.parametrize("latitude", [10.0, 30.0, 45.0, 60.0])
@pytest.mark.parametrize("hemisphere", ["north", "south"])
def test_latitude(latitude, hemisphere):
    """Check that radiation decreases away from the equator."""
    equinox = 81.0 / 365.0

    latitude *= -1.0 if hemisphere == "south" else 1.0

    grid = RasterModelGrid((5, 5), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")

    radiation = Radiation(grid, latitude=0.0, current_time=equinox)
    radiation.update()
    radiation_at_equator = grid.at_cell["radiation__net_shortwave_flux"].copy()

    radiation = Radiation(grid, latitude=latitude / 2, current_time=equinox)
    radiation.update()
    radiation_at_half_latitude = grid.at_cell["radiation__net_shortwave_flux"].copy()

    radiation = Radiation(grid, latitude=latitude, current_time=equinox)
    radiation.update()
    radiation_at_latitude = grid.at_cell["radiation__net_shortwave_flux"]

    assert np.all(radiation_at_equator > radiation_at_half_latitude)
    assert np.all(radiation_at_half_latitude > radiation_at_latitude)
