#! /usr/bin/env python
import os
import sys

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup
from setuptools.command.build_ext import build_ext

is_windows = sys.platform.startswith("win")
is_macos = sys.platform == "darwin"
is_linux = sys.platform.startswith("linux")

openmp_prefix = os.environ.get("OPENMP_PREFIX")
if openmp_prefix:
    with_openmp = True
else:
    if is_macos:
        with_openmp = os.environ.get("WITH_OPENMP", "0") == "1"
    else:
        with_openmp = os.environ.get("WITH_OPENMP", "1") == "1"

include_dirs: list[str] = []
library_dirs: list[str] = []
libraries: list[str] = []
extra_compile_args: list[str] = []
extra_link_args: list[str] = []

if os.environ.get("LANDLAB_BUILD_RELEASE", "0") == "1":
    extra_compile_args += ["-O3"] if not is_windows else ["/O2"]
else:
    extra_compile_args += ["-O0"] if not is_windows else ["/Od"]

if openmp_prefix:
    include_dirs += [os.path.join(openmp_prefix, "include")]
    library_dirs += [os.path.join(openmp_prefix, "lib")]

if with_openmp and is_linux:
    extra_compile_args += ["-fopenmp"]
    extra_link_args += ["-fopenmp"]

if with_openmp and is_macos:
    extra_compile_args += ["-Xpreprocessor", "-fopenmp"]
    libraries += ["omp"]

if with_openmp and is_windows:
    extra_compile_args += ["/openmp"]


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
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            extra_link_args=extra_link_args,
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
