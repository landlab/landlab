#! /usr/bin/env python
"""
The Landlab

:package name: TheLandlab
:Version: 0.1.0
:release date: 2013-03-24
:authors:
  Greg Tucker,
  Nicole Gasparini,
  Erkan Istanbulluoglu,
  Daniel Hobley,
  Sai Nudurupati,
  Jordan Adams,
  Eric Hutton

:url: https://csdms.colorado.edu/trac/landlab

:license: MIT
"""
from landlab.framework.collections import Palette, Arena, NoProvidersError
from landlab.framework.decorators import Implements, ImplementsOrRaise
from landlab.framework.framework import Framework
from landlab.model_grid import (BoundaryCondition, ModelGrid,
                                RasterModelGrid)
