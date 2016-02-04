import sys
import importlib
import warnings

from .bmi_bridge import wrap_as_bmi


LANDLAB_COMPONENTS = [
    ('landlab.components.flexure', 'Flexure'),
    ('landlab.components.overland_flow', 'OverlandFlow'),
]


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
for mod, cls in LANDLAB_COMPONENTS:
    try:
        name, obj = import_landlab_component(mod, cls)
    except ImportError:
        pass
    else:
        setattr(sys.modules[__name__], name, obj)
        __all__.append(name)
