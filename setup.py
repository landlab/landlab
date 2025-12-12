#! /usr/bin/env python
import os
import sys

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup
from setuptools.command.build_ext import build_ext

if os.environ.get("LANDLAB_BUILD_RELEASE", "0") == "1":
    extra_compile_args = ["-O3"] if sys.platform != "win32" else ["/O2"]
else:
    extra_compile_args = ["-O0"] if sys.platform != "win32" else ["/Od"]


def _build_cpu_count():
    try:
        return int(os.environ["LANDLAB_BUILD_CPU_COUNT"])
    except KeyError:
        return os.cpu_count() or 1


with open("cython-files.txt") as fp:
    cython_files = {fname.strip() for fname in fp.readlines()}

ext_modules = cythonize(
    [
        Extension(
            path[4:-4].replace("/", "."),
            [path],
            define_macros=[("NPY_NO_DEPRECATED_API", "1")],
            extra_compile_args=extra_compile_args,
        )
        for path in cython_files
    ],
    compiler_directives={"embedsignature": True, "language_level": 3},
    nthreads=_build_cpu_count(),
)


class build_ext_parallel(build_ext):
    def finalize_options(self):
        super().finalize_options()
        if self.parallel is None:
            self.parallel = _build_cpu_count()


setup(
    include_dirs=[np.get_include()],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext_parallel},
)
