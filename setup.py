#! /usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup
import multiprocessing

from setuptools.command.install import install
from setuptools.command.develop import develop


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


setup(name='TheLandlab',
      version='0.1.0',
      author='Eric Hutton',
      author_email='eric.hutton@colorado.edu',
      url='https://csdms.colorado.edu/trac/landlab',
      description='Plugin-based component modeling tool.',
      long_description=open('README.rst').read(),
      setup_requires=['numpy>=1.7',
                      'scipy>=0.12',
                      'nose>=1.0',
                      'memory_profiler'],
      install_requires=['numpy>=1.7',
                        'scipy>=0.12',
                        'nose>=1.0',
                        'memory_profiler'],
      packages=['landlab',
                'landlab.components',
                'landlab.grid',
                'landlab.io',
                'landlab.tests',
                'landlab.utils',
                'landlab.examples',
                'landlab.framework',
               ],
      entry_points={
          'console_scripts': [
              'craters = landlab.components.craters:main',
              'landlab = landlab.cmd:main',
              'landlab_ex_1 = landlab.examples.diffusion2Dtest:main',
              'landlab_ex_2 = landlab.examples.diffusion2Dtest2:main',
              'landlab_ex_3 = landlab.examples.diffusion2Dtest3:main',
          ]
      },
      package_data={'': ['data/*asc']},
      test_suite='nose.collector',
      cmdclass={
          'install': install_and_register,
          'develop': develop_and_register,
      },
     )
