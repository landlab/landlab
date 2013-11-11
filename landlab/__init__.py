#! /usr/bin/env python
"""
The Landlab

:Package name: TheLandlab
:Version: 0.1.0
:Release date: 2013-03-24
:Authors:
  Greg Tucker,
  Nicole Gasparini,
  Erkan Istanbulluoglu,
  Daniel Hobley,
  Sai Nudurupati,
  Jordan Adams,
  Eric Hutton

:URL: http://csdms.colorado.edu/trac/landlab

:License: MIT
"""

import os
if 'DISPLAY' not in os.environ:
    import matplotlib
    matplotlib.use('Agg')

from landlab.model_parameter_dictionary import ModelParameterDictionary
from landlab.framework.collections import Palette, Arena, NoProvidersError
from landlab.framework.decorators import Implements, ImplementsOrRaise
from landlab.framework.framework import Framework
from landlab.grid import *
from landlab.plot import *
from landlab.model_component import Component
