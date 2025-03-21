from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np


## Previous version
setup(
    name='cfuncs',
    ext_modules=cythonize("cfuncs.pyx"),
    include_dirs=[np.get_include()],
    zip_safe=False,
)



# ##
# ext_modules = [
#     Extension(
#         "raster_divergence",
#         ["raster_divergence.pyx"],
#         extra_compile_args=['-fopenmp'],
#         extra_link_args=['-fopenmp'],
#     )
# ]
#
#
# setup(
#     ext_modules=cythonize(ext_modules),
#     include_dirs=[np.get_include()]
# )
