import matplotlib.colors as colors
import numpy as np
import pytest
from matplotlib.colors import Normalize

from landlab.plot import plot_network_and_parcels

network_norm = Normalize(-1, 6)
parcel_color_norm = Normalize(0, 1)
parcel_color_norm2 = colors.LogNorm(vmin=0.01, vmax=1)
parcel_size_norm = Normalize(-1, 1)
parcel_size_norm2 = colors.LogNorm(vmin=0.01, vmax=1)


link_color_options = [
    {},  # empty dictionary = defaults
    {
        "network_color": "r",  # specify some simple modifications.
        "network_linewidth": 0.7,
        "parcel_alpha": 0,
    },
    {
        "link_attribute": "sediment_total_volume",  # use a link attribute
        "parcel_alpha": 0,
    },
    {
        "link_attribute": "sediment_total_volume",
        "network_cmap": "jet",  # change colormap
        "network_norm": network_norm,  # and normalize
        "link_attribute_title": "Total Sediment Volume",
        "parcel_alpha": 0,
        "network_linewidth": 3,
    },
]

parcel_color_options = [
    {},  # empty dictionary = defaults
    {"parcel_color": "r", "parcel_size": 10},  # specify some simple modifications.
    {
        "parcel_color_attribute": "D",  # use a parcel attribute.
        "parcel_color_norm": parcel_color_norm,
        "parcel_color_attribute_title": "Diameter [m]",
        "parcel_alpha": 1.0,
    },
    {"parcel_color_attribute": "abrasion_rate", "parcel_color_cmap": "bone"},
]

parcel_size_options = [
    {},  # empty dictionary = defaults
    {"parcel_color": "b", "parcel_size": 5},  # specify some simple modifications.
    {
        "parcel_size_attribute": "D",  # use a parcel attribute.
        "parcel_size_norm": parcel_size_norm2,
        "parcel_size_attribute_title": "Diameter [m]",
        "parcel_alpha": 1.0,
    },
    {
        "parcel_size_attribute": "abrasion_rate",  # an
        "parcel_size_min": 10,
        "parcel_size_max": 100,
    },
]


@pytest.mark.parametrize("arg", ["synthetic", "methow"])
@pytest.mark.parametrize(
    ("l_opts", "pc_opts", "ps_opts"),
    zip(link_color_options, parcel_color_options, parcel_size_options),
)
def test_link_options(arg, l_opts, pc_opts, ps_opts, request):
    nst = request.getfixturevalue(arg)
    grid = nst.grid
    parcels = nst._parcels
    opts = {**l_opts, **pc_opts, **ps_opts}
    plot_network_and_parcels(grid, parcels, parcel_time_index=0, **opts)


@pytest.mark.parametrize("title", ["A random number", None])
@pytest.mark.parametrize("arg", ["synthetic", "methow"])
def test_link_array(arg, title, request):
    nst = request.getfixturevalue(arg)
    grid = nst.grid
    parcels = nst._parcels

    random_link = np.random.randn(grid.size("link"))

    opts = {
        "link_attribute": random_link,  # use an array of size link.
        "network_cmap": "jet",  # change colormap
        "network_norm": network_norm,  # and normalize
        "link_attribute_title": title,
        "parcel_alpha": 0,
        "network_linewidth": 3,
    }
    plot_network_and_parcels(grid, parcels, parcel_time_index=0, **opts)


@pytest.mark.parametrize("arg", ["synthetic", "methow"])
def test_with_filter(arg, request):
    nst = request.getfixturevalue(arg)
    grid = nst.grid
    parcels = nst._parcels
    parcel_filter = np.zeros((parcels.dataset.dims["item_id"]), dtype=bool)
    parcel_filter[::10] = True
    plot_network_and_parcels(
        grid,
        parcels,
        parcel_time_index=0,
        parcel_filter=parcel_filter,
        link_attribute="sediment_total_volume",
        network_norm=network_norm,
        parcel_alpha=1.0,
        parcel_size_attribute="D",
        parcel_color_attribute="D",
        parcel_color_norm=parcel_color_norm2,
        parcel_size_norm=parcel_size_norm,
        parcel_size_attribute_title="D",
    )


def test_double_network_color(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    with pytest.raises(ValueError):
        plot_network_and_parcels(
            grid, parcels, link_attribute="sediment_total_volume", network_color="r"
        )


def test_double_parcel_color(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    with pytest.raises(ValueError):
        plot_network_and_parcels(
            grid, parcels, parcel_color_attribute="D", parcel_color="r"
        )


def test_double_parcel_size(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    with pytest.raises(ValueError):
        plot_network_and_parcels(
            grid, parcels, parcel_size_attribute="D", parcel_size=3
        )


def test_categorical_parcel_color(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    "quartzite"
    with pytest.raises(ValueError):
        plot_network_and_parcels(grid, parcels, parcel_color_attribute="quartzite")


def test_categorical_parcel_size(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    "quartzite"
    with pytest.raises(ValueError):
        plot_network_and_parcels(grid, parcels, parcel_size_attribute="quartzite")


def test_missing_parcel_color(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    "quartzite"
    with pytest.raises(ValueError):
        plot_network_and_parcels(grid, parcels, parcel_color_attribute="not_here")


def test_missing_parcel_size(synthetic):
    grid = synthetic.grid
    parcels = synthetic._parcels
    "quartzite"
    with pytest.raises(ValueError):
        plot_network_and_parcels(grid, parcels, parcel_size_attribute="not_here")
