#! /usr/bin/env python
import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

with open("cython-files.txt") as fp:
    cython_files = {fname.strip() for fname in fp.readlines()}

ext_modules = cythonize(
    [
        Extension(
            path[4:-4].replace("/", "."),
            [path],
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
