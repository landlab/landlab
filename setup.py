#! /usr/bin/env python

#from ez_setup import use_setuptools
#use_setuptools()

from setuptools import setup, find_packages, Extension
from setuptools.command.install import install
from setuptools.command.develop import develop
from distutils.extension import Extension
import pkg_resources

import versioneer


numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')


ext_modules = [
    Extension('landlab.ca.cfuncs',
              ['landlab/ca/cfuncs.pyx']),
    Extension('landlab.grid.cfuncs',
              ['landlab/grid/cfuncs.pyx']),
    Extension('landlab.components.flexure.cfuncs',
              ['landlab/components/flexure/cfuncs.pyx']),
    Extension('landlab.components.flow_accum.cfuncs',
              ['landlab/components/flow_accum/cfuncs.pyx']),
    Extension('landlab.components.flow_director.cfuncs',
              ['landlab/components/flow_director/cfuncs.pyx']),
    Extension('landlab.components.stream_power.cfuncs',
              ['landlab/components/stream_power/cfuncs.pyx']),
    Extension('landlab.components.space.cfuncs',
              ['landlab/components/space/cfuncs.pyx']),
    Extension('landlab.components.drainage_density.cfuncs',
              ['landlab/components/drainage_density/cfuncs.pyx']),
    Extension('landlab.components.erosion_deposition.cfuncs',
              ['landlab/components/erosion_deposition/cfuncs.pyx']),
    Extension('landlab.utils.ext.jaggedarray',
              ['landlab/utils/ext/jaggedarray.pyx']),
    Extension('landlab.graph.structured_quad.ext.at_node',
              ['landlab/graph/structured_quad/ext/at_node.pyx']),
    Extension('landlab.graph.structured_quad.ext.at_link',
              ['landlab/graph/structured_quad/ext/at_link.pyx']),
    Extension('landlab.graph.structured_quad.ext.at_patch',
              ['landlab/graph/structured_quad/ext/at_patch.pyx']),
    Extension('landlab.graph.structured_quad.ext.at_cell',
              ['landlab/graph/structured_quad/ext/at_cell.pyx']),
    Extension('landlab.graph.structured_quad.ext.at_face',
              ['landlab/graph/structured_quad/ext/at_face.pyx']),
    Extension('landlab.graph.hex.ext.hex',
              ['landlab/graph/hex/ext/hex.pyx']),
    Extension('landlab.graph.sort.ext.remap_element',
              ['landlab/graph/sort/ext/remap_element.pyx']),
    Extension('landlab.graph.sort.ext.argsort',
              ['landlab/graph/sort/ext/argsort.pyx']),
    Extension('landlab.graph.sort.ext.spoke_sort',
              ['landlab/graph/sort/ext/spoke_sort.pyx']),
    Extension('landlab.graph.voronoi.ext.voronoi',
              ['landlab/graph/voronoi/ext/voronoi.pyx']),
    Extension('landlab.graph.voronoi.ext.delaunay',
              ['landlab/graph/voronoi/ext/delaunay.pyx']),
    Extension('landlab.graph.object.ext.at_node',
              ['landlab/graph/object/ext/at_node.pyx']),
    Extension('landlab.graph.object.ext.at_patch',
              ['landlab/graph/object/ext/at_patch.pyx']),
    Extension('landlab.graph.quantity.ext.of_link',
              ['landlab/graph/quantity/ext/of_link.pyx']),
    Extension('landlab.graph.quantity.ext.of_patch',
              ['landlab/graph/quantity/ext/of_patch.pyx']),
    Extension('landlab.graph.matrix.ext.matrix',
              ['landlab/graph/matrix/ext/matrix.pyx']),
    Extension('landlab.grid.structured_quad.cfuncs',
              ['landlab/grid/structured_quad/cfuncs.pyx']),
    Extension('landlab.grid.structured_quad.c_faces',
              ['landlab/grid/structured_quad/c_faces.pyx']),
    Extension('landlab.layers.ext.eventlayers',
              ['landlab/layers/ext/eventlayers.pyx']),
]


def register(**kwds):
    import httplib, urllib

    data = urllib.urlencode(kwds)
    header = {"Content-type": "application/x-www-form-urlencoded",
              "Accept": "text/plain"}
    conn = httplib.HTTPConnection('csdms.colorado.edu')
    conn.request('POST', '/register/', data, header)


def register_landlab():
    try:
        from sys import argv
        import platform
        data = {
            'name': 'landlab',
            'version': __version__,
            'platform': platform.platform(),
            'desc': ';'.join(argv),
        }
        register(**data)
    except Exception:
        pass


class install_and_register(install):
    def run(self):
        install.run(self)
        register_landlab()


class develop_and_register(develop):
    def run(self):
        develop.run(self)
        register_landlab()


import os

#cython_pathspec = os.path.join('landlab', 'components','**','*.pyx')
#ext_modules = cythonize(cython_pathspec)


setup(name='landlab',
      version=versioneer.get_version(),
      author='Eric Hutton',
      author_email='eric.hutton@colorado.edu',
      url='https://github.com/landlab',
      description='Plugin-based component modeling tool.',
      long_description=open('README.rst').read(),
      install_requires=['scipy>=0.12',
                        'nose>=1.3',
                        'matplotlib',
                        'sympy',
                        'pandas',
                        'six',
                        'pyyaml',
                        'netCDF4',
                        'xarray',
                       ],
      #                  'Cython>=0.22'],
      setup_requires=['cython'],
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Programming Language :: Cython',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: Implementation :: CPython',
          'Topic :: Scientific/Engineering :: Physics'
      ],
      packages=find_packages(),
      package_data={'': ['tests/*txt', 'data/*asc', 'data/*nc',
                         'preciptest.in']},
      test_suite='nose.collector',
      cmdclass=versioneer.get_cmdclass({
          'install': install_and_register,
          'develop': develop_and_register,
      }),
      entry_points={
          'console_scripts': [
              'landlab=landlab.cmd.landlab:main',
          ]
      },
      include_dirs = [numpy_incl, ],
      ext_modules = ext_modules,
     )
