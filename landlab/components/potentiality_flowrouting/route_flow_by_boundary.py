# -*- coding: utf-8 -*-
"""
This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router and diffuser. It needs to be broken into routing
and sed mobility methods post hoc.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

import numpy as np
from landlab import Component

class PotentialityFlowRouter(Component):
    """
    """
    _name = 'PotentialityFlowRouter'
    
    def __init__(self, grid, params):
        self.initialize(grid, params)
    
    def initialize(self, grid, params):
        pass
    
    def route_flow(self, dt):
        """
        """
        