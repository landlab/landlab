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

from numpy import set_printoptions
try:
    set_printoptions(legacy='1.13')
except TypeError:
    pass
finally:
    del set_printoptions

from ._registry import registry

cite_as = registry.format_citations

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

__all__.extend(['ModelParameterDictionary', 'MissingKeyError',
                'ParameterValueError', 'Component', 'Palette', 'Arena',
                'NoProvidersError', 'Implements', 'ImplementsOrRaise',
                'Framework', 'FieldError', 'LandlabTester', 'load_params'])

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
