import sys
import importlib
import warnings

from .bmi_bridge import wrap_as_bmi
from ..components import COMPONENTS


def import_landlab_component(mod_name, cls_name):
    try:
        landlab_module = importlib.import_module(mod_name)
    except ImportError:
        warnings.warn('unable to import {mod}'.format(mod=mod_name))
        raise
    else:
        obj = landlab_module.__dict__[cls_name]
        as_bmi = wrap_as_bmi(obj)
        return as_bmi.__name__, as_bmi
        # setattr(sys.modules[__name__], as_bmi.__name__, as_bmi)
        # return as_bmi.__name__


__all__ = []
for cls in COMPONENTS:
    try:
        as_bmi = wrap_as_bmi(cls)
    except TypeError:
        warnings.warn('unable to wrap class {name}'.format(name=cls.__name__))
    else:
        setattr(sys.modules[__name__], as_bmi.__name__, as_bmi)
        __all__.append(as_bmi.__name__)
