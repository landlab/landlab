from functools import partial
from itertools import tee

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.patches import Patch
from scipy.interpolate import interp1d


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def plot_layers(
    elevation_at_layer,
    x=None,
    sea_level=0.0,
    color_water=(0.8, 1.0, 1.0),
    color_bedrock=(0.8, 0.8, 0.8),
    color_layer=None,
    layer_start=0,
    layer_stop=-1,
    n_layers=5,
    layer_line_width=0.5,
    layer_line_color="k",
    title=None,
    x_label="Distance",
    y_label="Elevation",
    legend_location="lower left",
):
    elevation_at_layer = np.asarray(elevation_at_layer)
    elevation_at_layer = np.expand_dims(
        elevation_at_layer, axis=tuple(np.arange(2 - elevation_at_layer.ndim))
    )

    n_layers = np.minimum(n_layers, len(elevation_at_layer))
    legend_item = partial(Patch, edgecolor="k", linewidth=0.5)

    if x is None:
        x = np.arange(elevation_at_layer.shape[1])

    top_surface = elevation_at_layer[-1]
    bottom_surface = elevation_at_layer[0]

    x_water, y_water = _insert_shorelines(x, top_surface, sea_level=sea_level)
    if color_water:
        plt.fill_between(
            x_water, y_water, sea_level, where=y_water <= sea_level, fc=color_water
        )
    water_surface = np.full_like(x_water, sea_level, dtype=float)
    water_surface[y_water > sea_level] = np.nan
    plt.plot(x_water, water_surface, color="b")

    plt.fill_between(
        x,
        bottom_surface,
        np.full_like(bottom_surface, bottom_surface.min()),
        color=color_bedrock,
    )
    plt.plot(x, bottom_surface, color="k")
    plt.plot(x, top_surface, color="g")

    if layer_stop < 0:
        layer_stop = len(elevation_at_layer) + layer_stop + 1

    layers_to_plot = np.linspace(layer_start, layer_stop - 1, n_layers, dtype=int)
    if len(layers_to_plot):
        plt.plot(
            x,
            elevation_at_layer[layers_to_plot].T,
            color=layer_line_color,
            linewidth=layer_line_width,
        )

    if color_layer is not None and len(layers_to_plot):
        cmap = cm.get_cmap(color_layer)
        for layer, (lower, upper) in enumerate(
            pairwise(elevation_at_layer[layers_to_plot])
        ):
            plt.fill_between(
                x,
                lower,
                upper,
                fc=cmap(layers_to_plot[layer] * 256 // len(elevation_at_layer)),
            )

    if legend_location:
        items = [
            ("Ocean", color_water),
            ("Bedrock", color_bedrock),
        ]
        legend = [legend_item(label=label, fc=color) for label, color in items if color]
        legend and plt.legend(handles=legend, loc=legend_location)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    title and plt.title(title)

    plt.show()


def _insert_shorelines(x, y, sea_level=0.0):
    """Insert shorelines into x-y arrays.

    Examples
    --------
    >>> from landlab.plot.layers import _insert_shorelines
    >>> _insert_shorelines([0, 1, 2],[2, 1, -1])
    (array([ 0. ,  1. ,  1.5,  2. ]), array([ 2.,  1.,  0., -1.]))
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)

    y_relative_to_sea_level = y - sea_level
    shorelines = _search_shorelines(y_relative_to_sea_level)
    x_of_shoreline = _interp_shorelines(x, y_relative_to_sea_level, shorelines)

    return (
        np.insert(x, shorelines + 1, x_of_shoreline),
        np.insert(y, shorelines + 1, sea_level),
    )


def _search_shorelines(y):
    """Search an array for changes in sign.

    Parameters
    ----------
    y : array-like
        Input array to check for sign changes.

    Returns
    -------
    int
        Indices into *y* where a sign has changed.

    Examples
    --------
    >>> from landlab.plot.layers import _search_shorelines
    >>> _search_shorelines([2, 1, -1])
    array([1])
    >>> _search_shorelines([2, 0, 0, -2])
    array([2])
    >>> _search_shorelines([-2, -2, 1, 2])
    array([1])
    >>> len(_search_shorelines([2, 0, 1])) == 0
    True
    >>> len(_search_shorelines([2, 3, 4])) == 0
    True
    >>> len(_search_shorelines([0, 0, 0])) == 0
    True
    """
    sign = np.sign(y)

    zeros = sign == 0
    if not np.all(zeros):
        while np.any(zeros):
            sign[zeros] = np.roll(sign, 1)[zeros]
            zeros = sign == 0

    return np.where(sign[1:] != sign[:-1])[0]


def _interp_shorelines(x, y, shorelines):
    """Interpolate between adjacent elements to find shoreline positions.

    Parameters
    ----------
    x : array-like
        Distances.
    y : array-like
        Elevations.
    shorelines : array-like of int
        Indices to shoreline elements.

    Returns
    -------
    array of float
        Distances to interpolated shorelines.

    Examples
    --------
    >>> from landlab.plot.layers import _interp_shorelines
    >>> _interp_shorelines([0, 1, 2], [1, -1, -1], [0])
    array([ 0.5])
    >>> _interp_shorelines([0, 1, 2, 3], [1, -1, -1, 4], [0, 2])
    array([ 0.5, 2.2])
    """
    x_of_shoreline = []
    for shoreline in shorelines:
        coast = slice(shoreline, shoreline + 2)
        x_of_shoreline.append(interp1d(y[coast], x[coast])(0.0))

    return np.asarray(x_of_shoreline)
