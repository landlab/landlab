#! /usr/bin/env python

import os
import re

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy


def find_extensions(paths=["."]):
    extensions = []
    for path in paths:
        for root, _dirs, files in os.walk(os.path.normpath(path)):
            _dirs[:] = [d for d in _dirs if not d[0] == "."]
            extensions += [
                os.path.join(root, fname) for fname in files if fname.endswith(".pyx")
            ]
    return [
        Extension(re.sub(re.escape(os.path.sep), ".", ext[: -len(".pyx")]), [ext])
        for ext in extensions
    ]


setup(
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize(
        module_list=find_extensions(["landlab", "tests"]),
        compiler_directives={"embedsignature": True, "language_level": 3},
    ),
)
