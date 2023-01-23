#! /usr/bin/env python
"""The Landlab.

:Package name: TheLandlab
:Release date: 2018-09-18
:Authors: Greg Tucker, Nicole Gasparini, Erkan Istanbulluoglu, Daniel Hobley,
    Sai Nudurupati, Jordan Adams, Eric Hutton, Katherine Barnhart, Margaux
    Mouchene, Nathon Lyons
:URL: https://landlab.readthedocs.io/en/release/
:License: MIT
"""
import contextlib

from numpy import set_printoptions

from ._registry import registry
from ._version import __version__
from .core.errors import MissingKeyError, ParameterValueError
from .core.model_component import Component
from .core.model_parameter_loader import load_params
from .core.utils import ExampleData
from .field import FieldError
from .grid import (
    FramedVoronoiGrid,
    HexModelGrid,
    ModelGrid,
    NetworkModelGrid,
    RadialModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
    create_grid,
)
from .grid.linkstatus import LinkStatus
from .grid.nodestatus import NodeStatus
from .plot import imshow_grid, imshow_grid_at_node, imshowhs_grid, imshowhs_grid_at_node

with contextlib.suppress(TypeError):
    set_printoptions(legacy="1.13")
del set_printoptions

cite_as = registry.format_citations

__all__ = [
    "__version__",
    "registry",
    "MissingKeyError",
    "ParameterValueError",
    "Component",
    "FieldError",
    "load_params",
    "ExampleData",
    "ModelGrid",
    "HexModelGrid",
    "RadialModelGrid",
    "RasterModelGrid",
    "FramedVoronoiGrid",
    "VoronoiDelaunayGrid",
    "NetworkModelGrid",
    "LinkStatus",
    "NodeStatus",
    "create_grid",
    "imshow_grid",
    "imshow_grid_at_node",
    "imshowhs_grid",
    "imshowhs_grid_at_node",
]
