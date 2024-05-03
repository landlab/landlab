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
from landlab._registry import registry
from landlab._version import __version__
from landlab.core.errors import MissingKeyError
from landlab.core.errors import ParameterValueError
from landlab.core.model_component import Component
from landlab.core.model_parameter_loader import load_params
from landlab.core.utils import ExampleData
from landlab.field.errors import FieldError
from landlab.grid.base import ModelGrid
from landlab.grid.create import create_grid
from landlab.grid.framed_voronoi import FramedVoronoiGrid
from landlab.grid.hex import HexModelGrid
from landlab.grid.icosphere import IcosphereGlobalGrid
from landlab.grid.linkstatus import LinkStatus
from landlab.grid.network import NetworkModelGrid
from landlab.grid.nodestatus import NodeStatus
from landlab.grid.radial import RadialModelGrid
from landlab.grid.raster import RasterModelGrid
from landlab.grid.voronoi import VoronoiDelaunayGrid
from landlab.plot.imshow import imshow_grid
from landlab.plot.imshow import imshow_grid_at_node
from landlab.plot.imshowhs import imshowhs_grid
from landlab.plot.imshowhs import imshowhs_grid_at_node

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
    "IcosphereGlobalGrid",
    "LinkStatus",
    "NodeStatus",
    "create_grid",
    "imshow_grid",
    "imshow_grid_at_node",
    "imshowhs_grid",
    "imshowhs_grid_at_node",
]
