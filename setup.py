#! /usr/bin/env python

import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

compile_args = [] if "LANDLAB_WITHOUT_OPENMP" in os.environ else ["-fopenmp"]

cython_files = (
    "src/landlab/ca/cfuncs.pyx",
    "src/landlab/components/bedrock_landslider/cfuncs.pyx",
    "src/landlab/components/depression_finder/cfuncs.pyx",
    "src/landlab/components/drainage_density/cfuncs.pyx",
    "src/landlab/components/erosion_deposition/cfuncs.pyx",
    "src/landlab/components/flexure/_ext/flexure1d.pyx",
    "src/landlab/components/flexure/_ext/flexure2d.pyx",
    "src/landlab/components/flexure/_ext/flexure2d_slow.pyx",
    "src/landlab/components/flow_accum/cfuncs.pyx",
    "src/landlab/components/flow_director/cfuncs.pyx",
    "src/landlab/components/flow_router/ext/single_flow/priority_routing/breach.pyx",
    "src/landlab/components/overland_flow/_neighbors_at_link.pyx",
    "src/landlab/components/priority_flood_flow_router/cfuncs.pyx",
    "src/landlab/components/space/ext/calc_qs.pyx",
    "src/landlab/components/space/ext/calc_sequential_ero_depo.pyx",
    "src/landlab/components/stream_power/cfuncs.pyx",
    "src/landlab/components/threshold_eroder/cfuncs.pyx",
    "src/landlab/graph/hex/ext/hex.pyx",
    "src/landlab/graph/hex/ext/perimeternodes.pyx",
    "src/landlab/graph/matrix/ext/at_patch.pyx",
    "src/landlab/graph/matrix/ext/matrix.pyx",
    "src/landlab/graph/object/ext/at_node.pyx",
    "src/landlab/graph/object/ext/at_patch.pyx",
    "src/landlab/graph/quantity/ext/of_element.pyx",
    "src/landlab/graph/quantity/ext/of_link.pyx",
    "src/landlab/graph/quantity/ext/of_patch.pyx",
    "src/landlab/graph/sort/ext/_deprecated_sparse.pyx",
    "src/landlab/graph/sort/ext/argsort.pyx",
    "src/landlab/graph/sort/ext/intpair.pyx",
    "src/landlab/graph/sort/ext/remap_element.pyx",
    "src/landlab/graph/sort/ext/spoke_sort.pyx",
    "src/landlab/graph/structured_quad/ext/at_cell.pyx",
    "src/landlab/graph/structured_quad/ext/at_face.pyx",
    "src/landlab/graph/structured_quad/ext/at_link.pyx",
    "src/landlab/graph/structured_quad/ext/at_node.pyx",
    "src/landlab/graph/structured_quad/ext/at_patch.pyx",
    "src/landlab/graph/voronoi/ext/delaunay.pyx",
    "src/landlab/graph/voronoi/ext/voronoi.pyx",
    "src/landlab/grid/ext/raster_divergence.pyx",
    "src/landlab/grid/ext/raster_gradient.pyx",
    "src/landlab/layers/ext/eventlayers.pyx",
    "src/landlab/utils/_matrix.pyx",
    "src/landlab/utils/ext/jaggedarray.pyx",
    # "tests/components/flow_router/ext/single_flow/priority_routing/test_breach_c.pyx",
)

ext_modules = cythonize(
    [
        Extension(
            path[:-4].replace("/", "."),
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
