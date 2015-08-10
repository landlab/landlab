#!/usr/bin/env python
"""
.. codeauthor:: Eric Hutton <huttone@colorado.edu>

.. sectionauthor:: Eric Hutton <huttone@colorado.edu>
"""

from landlab.components.flexure.flexure import FlexureComponent
from landlab.components.flexure.funcs import (get_flexure_parameter,
                                              subside_point_load,
                                              subside_point_loads)

__all__ = ['FlexureComponent', 'get_flexure_parameter',
           'subside_point_load', 'subside_point_loads']
