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
import os

# from ._info import version as  __version__
from ._registry import registry

__all__ = ['registry']

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
from .core.model_parameter_loader import load_params
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

__all__.extend(['ModelParameterDictionary', 'MissingKeyError',
                'ParameterValueError', 'Component', 'Palette', 'Arena',
                'NoProvidersError', 'Implements', 'ImplementsOrRaise',
                'Framework', 'FieldError', 'LandlabTester', 'load_params'])

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
