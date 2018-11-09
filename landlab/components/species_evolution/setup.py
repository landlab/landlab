#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  5 12:42:59 2018

@author: njlyons
"""

from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(
    ext_modules=cythonize('cfuncs.pyx'),
    include_dirs=[np.get_include()]
)
