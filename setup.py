#! /usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

from distutils.extension import Extension

import sys
if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

from Cython.Build import cythonize
import numpy as np

from landlab import __version__


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

cython_pathspec = os.path.join('landlab', 'components','**','*.pyx')
ext_modules = cythonize(cython_pathspec)


setup(name='landlab',
      version=__version__,
      author='Eric Hutton',
      author_email='eric.hutton@colorado.edu',
      url='https://csdms.colorado.edu/trac/landlab',
      description='Plugin-based component modeling tool.',
      long_description=open('README.md').read(),
      #setup_requires=['numpy>=1.7',
      #                'scipy>=0.12',
      #                'nose>=1.0',
      #                'numpydoc' ],
      #install_requires=['numpy>=1.7',
      #                  'scipy>=0.12',
      #                  'nose>=1.0',
      #                  'numpydoc'],
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'craters = landlab.components.craters:main',
              'landlab = landlab.cmd:main',
              'landlab_ex_1 = landlab.examples.diffusion2Dtest:main',
              'landlab_ex_2 = landlab.examples.diffusion2Dtest2:main',
              'landlab_ex_3 = landlab.examples.diffusion2Dtest3:main',
          ]
      },
      package_data={'': ['data/*asc', 'preciptest.in']},
      test_suite='nose.collector',
      cmdclass={
          'install': install_and_register,
          'develop': develop_and_register,
      },

      include_dirs = [np.get_include()],
      ext_modules = ext_modules,
     )
