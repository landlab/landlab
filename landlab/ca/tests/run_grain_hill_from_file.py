# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 09:13:46 2016

@author: gtucker
"""

import grain_hill_as_class
import sys

try:
    infile = sys.argv[1]
except:
    print('You need to specify the name of an input file')
    raise

params = grain_hill_as_class.get_params_from_input_file(infile)
grain_hill_as_class.main(params)
