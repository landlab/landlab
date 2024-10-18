from functools import partial
from itertools import tee

import matplotlib.pyplot as plt
import numpy as np
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
    layer_line_width=0.5,
    layer_line_color="k",
    title=None,
    x_label="Distance",
    y_label="Elevation",
    legend_location="lower left",
):
    """Plot a stack of sediment layers as a cross section.

    Create a plot of the elevation sediment layers, including surfaces for
    sea level and bedrock.

    Parameters
    ----------
    elevation_at_layer : array-like of shape *(n_layers, n_stacks)*
        Elevation to each layer along the profile. Layers are provided
        row-by-row, with the bottom-most layer being the first row.
    x : array-like, optional
        Distance to each stack along the cross-section. If not provided,
        stack number will be used.
    sea_level : float, optional
        Elevation of sea level.
    color_water : tuple of float, optional
        Tuple of *(red, green, blue)* values for water.
    color_bedrock : tuple of float, optional
        Tuple of *(red, green, blue)* values for bedrock.
    color_layer : string, optional
        Colormap to use to color in the layers.
    layer_line_width : float, optional
        Width of line used to plot layer surfaces.
    layer_line_color : string, optional
        Color of the line used to plot layer surfaces.
    title : string, optional
        Text to be used for the graph's title. The default is to not
        include a title.
    x_label : string, optional
        Text to be used for the x (horizontal) axis label.
    y_label : string, optional
        Text to be used for the y (vertical) axis label.
    legend_location : string, optional
        Where to put the legend.
    """
    elevation_at_layer = np.asarray(elevation_at_layer)
    elevation_at_layer = np.expand_dims(
        elevation_at_layer,
        axis=tuple(np.arange(2 - elevation_at_layer.ndim)),
    )

    if len(elevation_at_layer) == 0:
        raise ValueError(
            f"no layers to plot (elevation_at_layer.shape is {np.shape(elevation_at_layer)}"
        )

    if x is None:
        x = np.arange(elevation_at_layer.shape[1])

    top_surface = elevation_at_layer[-1]
    bottom_surface = elevation_at_layer[0]

    if len(elevation_at_layer) > 0:
        _plot_layers(
            x,
            elevation_at_layer,  # [layers_to_plot],
            color=color_layer,
            lc=layer_line_color,
            lw=layer_line_width,
        )
    _plot_water(x, top_surface, sea_level=sea_level, fc=color_water)
    _plot_bedrock(x, bottom_surface, fc=color_bedrock)
    _plot_surface(x, top_surface, sea_level=sea_level)

    legend_location and _plot_legend(
        legend_location, color_water=color_water, color_bedrock=color_bedrock
    )

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    title and plt.title(title)

    plt.show()


def _plot_water(x, y, sea_level=0.0, fc=(0.8, 1.0, 1.0)):
    x_water, y_water = _insert_shorelines(x, y, sea_level=sea_level)
    if fc is not None:
        plt.fill_between(x_water, y_water, sea_level, where=y_water <= sea_level, fc=fc)

    water_surface = np.full_like(x_water, sea_level, dtype=float)
    water_surface[y_water > sea_level] = np.nan
    plt.plot(x_water, water_surface, color="b")


def _plot_bedrock(x, y, fc=(0.8, 0.8, 0.8)):
    if fc is not None:
        plt.fill_between(
            x,
            y,
            np.full_like(y, y.min()),
            color=fc,
        )
    plt.plot(x, y, color="k")


def _plot_surface(x, y, sea_level=0.0):
    under_water = y <= sea_level
    plt.plot(x[~under_water], y[~under_water], color="g")
    plt.plot(x[under_water], y[under_water], color="b")


def _plot_layers(x, layers, color=None, lc="k", lw=0.5):
    if color is not None:
        cmap = plt.colormaps[color] if isinstance(color, str) else color

        for layer, (lower, upper) in enumerate(pairwise(layers)):
            plt.fill_between(
                x,
                lower,
                upper,
                fc=cmap(layer * 256 // len(layers)),
            )
    plt.plot(
        x,
        layers.T,
        color=lc,
        linewidth=lw,
    )


def _plot_legend(legend_location, color_water=None, color_bedrock=None):
    legend_item = partial(Patch, edgecolor="k", linewidth=0.5)
    items = [
        ("Ocean", color_water),
        ("Bedrock", color_bedrock),
    ]
    legend = [legend_item(label=label, fc=color) for label, color in items if color]
    legend and plt.legend(handles=legend, loc=legend_location)


def _insert_shorelines(x, y, sea_level=0.0):
    """Insert shorelines into x-y arrays.

    Examples
    --------
    >>> from landlab.plot.layers import _insert_shorelines
    >>> _insert_shorelines([0, 1, 2], [2, 1, -1])
    (array([0. ,  1. ,  1.5,  2. ]), array([ 2.,  1.,  0., -1.]))
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)

    y_relative_to_sea_level = y - sea_level
    shorelines = _search_zero_crossings(y_relative_to_sea_level)
    x_of_shoreline = _interp_zero_crossings(x, y_relative_to_sea_level, shorelines)

    return (
        np.insert(x, shorelines + 1, x_of_shoreline),
        np.insert(y, shorelines + 1, sea_level),
    )


def _search_zero_crossings(y):
    """Search an array for changes in sign between elements.

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
    >>> from landlab.plot.layers import _search_zero_crossings

    The returned index is to the element before the zero-crossing.

    >>> list(_search_zero_crossings([2, 1, -1]))
    [1]
    >>> list(_search_zero_crossings([-2, -2, 1, 2]))
    [1]
    >>> list(_search_zero_crossings([-2, -2, 1, 2, -1]))
    [1, 3]

    These are not zero-crossings.

    >>> len(_search_zero_crossings([2, 0, 0, -2])) == 0
    True
    >>> len(_search_zero_crossings([2, 0, 1])) == 0
    True
    >>> len(_search_zero_crossings([2, 3, 4])) == 0
    True
    >>> len(_search_zero_crossings([0, 0, 0])) == 0
    True
    """
    sign = np.sign(y)

    # zeros = sign == 0
    # if not np.all(zeros):
    #     while np.any(zeros):
    #         sign[zeros] = np.roll(sign, 1)[zeros]
    #         zeros = sign == 0

    # return np.where(sign[1:] != sign[:-1])[0]
    return np.where(sign[1:] * sign[:-1] < 0)[0]


def _interp_zero_crossings(x, y, shorelines):
    """Interpolate between adjacent elements to find x-locations of zero-crossings.

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
    >>> from landlab.plot.layers import _interp_zero_crossings
    >>> _interp_zero_crossings([0, 1, 2], [1, -1, -1], [0])
    array([0.5])
    >>> _interp_zero_crossings([0, 1, 2, 3], [1, -1, -1, 4], [0, 2])
    array([0.5, 2.2])
    """
    x_of_shoreline = []
    for shoreline in shorelines:
        coast = slice(shoreline, shoreline + 2)

        # for scipy<1.10 interp1d requires x and y to have at least two elements,
        # which is not the case if theshoreline is the last element.
        x_of_shoreline.append(
            interp1d(np.broadcast_to(y[coast], 2), np.broadcast_to(x[coast], 2))(0.0)
        )

    return np.asarray(x_of_shoreline)
