#! /usr/bin/env python
"""The Landlab

:Package name: TheLandlab
:Release date: 2013-03-24
:Authors: Greg Tucker, Nicole Gasparini, Erkan Istanbulluoglu, Daniel Hobley,
    Sai Nudurupati, Jordan Adams, Eric Hutton
:URL: http://csdms.colorado.edu/trac/landlab
:License: MIT
"""

from __future__ import absolute_import

__version__ = '1.0.0-beta.1'


import os
if 'DISPLAY' not in os.environ:
    try:
        import matplotlib
    except ImportError:
        import warnings
        warnings.warn('matplotlib not found', ImportWarning)
    else:
        matplotlib.use('Agg')

from .core.model_parameter_dictionary import ModelParameterDictionary
from .core.model_parameter_dictionary import (MissingKeyError,
                                              ParameterValueError)
from .core.model_component import Component
from .framework.collections import Palette, Arena, NoProvidersError
from .framework.decorators import Implements, ImplementsOrRaise
from .framework.framework import Framework
from .field.scalar_data_fields import FieldError
from .grid import *
from .plot import *

from .testing.nosetester import LandlabTester
test = LandlabTester().test
bench = LandlabTester().bench

__all__ = ['ModelParameterDictionary', 'MissingKeyError',
           'ParameterValueError', 'Component', 'Palette', 'Arena',
           'NoProvidersError', 'Implements', 'ImplementsOrRaise',
           'Framework', 'FieldError', 'LandlabTester']


