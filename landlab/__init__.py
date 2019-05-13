#! /usr/bin/env python
"""The Landlab

:Package name: TheLandlab
:Release date: 2018-09-18
:Authors: Greg Tucker, Nicole Gasparini, Erkan Istanbulluoglu, Daniel Hobley,
    Sai Nudurupati, Jordan Adams, Eric Hutton, Katherine Barnhart, Margaux
    Mouchene, Nathon Lyons
:URL: http://csdms.colorado.edu/trac/landlab
:License: MIT
"""
from __future__ import absolute_import

from numpy import set_printoptions

from ._registry import registry
from ._version import get_versions
from .core.model_component import Component
from .core.model_parameter_dictionary import (
    MissingKeyError,
    ModelParameterDictionary,
    ParameterValueError,
)
from .core.model_parameter_loader import load_params
from .field.scalar_data_fields import FieldError
from .framework.collections import Arena, NoProvidersError, Palette
from .framework.decorators import Implements, ImplementsOrRaise
from .framework.framework import Framework
from .grid import (
    ACTIVE_LINK,
    BAD_INDEX_VALUE,
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_GRADIENT_BOUNDARY,
    FIXED_LINK,
    FIXED_VALUE_BOUNDARY,
    INACTIVE_LINK,
    LOOPED_BOUNDARY,
    HexModelGrid,
    ModelGrid,
    NetworkModelGrid,
    RadialModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
    create_and_initialize_grid,
    create_grid,
)
from .plot import analyze_channel_network_and_plot, imshow_grid, imshow_grid_at_node

try:
    set_printoptions(legacy="1.13")
except TypeError:
    pass
finally:
    del set_printoptions

cite_as = registry.format_citations

__all__ = [
    "registry",
    "ModelParameterDictionary",
    "MissingKeyError",
    "ParameterValueError",
    "Component",
    "Palette",
    "Arena",
    "NoProvidersError",
    "Implements",
    "ImplementsOrRaise",
    "Framework",
    "FieldError",
    "LandlabTester",
    "load_params",
    "ModelGrid",
    "HexModelGrid",
    "RadialModelGrid",
    "RasterModelGrid",
    "VoronoiDelaunayGrid",
    "NetworkModelGrid",
    "BAD_INDEX_VALUE",
    "CORE_NODE",
    "FIXED_VALUE_BOUNDARY",
    "FIXED_GRADIENT_BOUNDARY",
    "LOOPED_BOUNDARY",
    "CLOSED_BOUNDARY",
    "ACTIVE_LINK",
    "FIXED_LINK",
    "INACTIVE_LINK",
    "create_and_initialize_grid",
    "create_grid",
    "imshow_grid",
    "imshow_grid_at_node",
    "analyze_channel_network_and_plot",
]

__version__ = get_versions()["version"]
del get_versions
