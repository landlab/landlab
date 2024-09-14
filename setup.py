#! /usr/bin/env python

import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

compile_args = [] if "LANDLAB_WITHOUT_OPENMP" in os.environ else ["-fopenmp"]

with open("cython-files.txt") as fp:
    cython_files = {fname.strip() for fname in fp.readlines()}

ext_modules = cythonize(
    [
        Extension(
            path[4:-4].replace("/", "."),
            [path],
            extra_compile_args=compile_args,
            extra_link_args=compile_args,
            define_macros=[("NPY_NO_DEPRECATED_API", "1")],
        )
        for path in cython_files
    ],
    compiler_directives={"embedsignature": True, "language_level": 3},
)

setup(
    include_dirs=[np.get_include()],
    ext_modules=ext_modules,
)
