#! /usr/bin/env python

import os
import pathlib

import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup


def find_extensions(path="."):
    extensions = pathlib.Path(path).rglob("*.pyx")
    return [
        Extension(str(ext.with_suffix("")).replace(os.path.sep, "."), [str(ext)])
        for ext in extensions
    ]


setup(
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize(
        find_extensions("landlab") + find_extensions("tests"),
        compiler_directives={"embedsignature": True, "language_level": 3},
    ),
)
