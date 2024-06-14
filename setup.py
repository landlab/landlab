#! /usr/bin/env python

import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

compile_args = ["-fopenmp"] if "WITH_OPENMP" in os.environ else []

cython_files = (
    "landlab/ca/cfuncs.pyx",
    "landlab/components/bedrock_landslider/cfuncs.pyx",
    "landlab/components/depression_finder/cfuncs.pyx",
    "landlab/components/drainage_density/cfuncs.pyx",
    "landlab/components/erosion_deposition/cfuncs.pyx",
    "landlab/components/flexure/cfuncs.pyx",
    "landlab/components/flexure/ext/flexure1d.pyx",
    "landlab/components/flow_accum/cfuncs.pyx",
    "landlab/components/flow_director/cfuncs.pyx",
    "landlab/components/flow_router/ext/single_flow/priority_routing/breach.pyx",
    "landlab/components/overland_flow/_neighbors_at_link.pyx",
    "landlab/components/priority_flood_flow_router/cfuncs.pyx",
    "landlab/components/space/ext/calc_qs.pyx",
    "landlab/components/space/ext/calc_sequential_ero_depo.pyx",
    "landlab/components/stream_power/cfuncs.pyx",
    "landlab/components/threshold_eroder/cfuncs.pyx",
    "landlab/data_record/_aggregators.pyx",
    "landlab/graph/hex/ext/hex.pyx",
    "landlab/graph/hex/ext/perimeternodes.pyx",
    "landlab/graph/matrix/ext/at_patch.pyx",
    "landlab/graph/matrix/ext/matrix.pyx",
    "landlab/graph/object/ext/at_node.pyx",
    "landlab/graph/object/ext/at_patch.pyx",
    "landlab/graph/quantity/ext/of_element.pyx",
    "landlab/graph/quantity/ext/of_link.pyx",
    "landlab/graph/quantity/ext/of_patch.pyx",
    "landlab/graph/sort/ext/argsort.pyx",
    "landlab/graph/sort/ext/remap_element.pyx",
    "landlab/graph/sort/ext/spoke_sort.pyx",
    "landlab/graph/structured_quad/ext/at_cell.pyx",
    "landlab/graph/structured_quad/ext/at_face.pyx",
    "landlab/graph/structured_quad/ext/at_link.pyx",
    "landlab/graph/structured_quad/ext/at_node.pyx",
    "landlab/graph/structured_quad/ext/at_patch.pyx",
    "landlab/graph/voronoi/ext/delaunay.pyx",
    "landlab/graph/voronoi/ext/voronoi.pyx",
    "landlab/grid/cfuncs.pyx",
    "landlab/grid/ext/raster_divergence.pyx",
    "landlab/grid/ext/raster_gradient.pyx",
    "landlab/layers/ext/eventlayers.pyx",
    "landlab/utils/_matrix.pyx",
    "landlab/utils/ext/jaggedarray.pyx",
    "tests/components/flow_router/ext/single_flow/priority_routing/test_breach_c.pyx",
)

ext_modules = cythonize(
    [
        Extension(
            path[:-4].replace("/", "."),
            [path],
            extra_compile_args=compile_args,
            extra_link_args=compile_args,
        )
        for path in cython_files
    ],
    compiler_directives={"embedsignature": True, "language_level": 3},
)

setup(
    include_dirs=[np.get_include()],
    ext_modules=ext_modules,
)
