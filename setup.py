#! /usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup


setup(name='TheLandlab',
      version='0.1.0',
      author='Eric Hutton',
      author_email='eric.hutton@colorado.edu',
      url='https://csdms.colorado.edu/trac/landlab',
      description='Plugin-based component modeling tool.',
      long_description=open ('README.rst').read (),
      install_requires=['CmtBasicModelingInterface'],
      packages=['landlab',
                'landlab.plugins',
                'landlab.tests',
                'landlab.examples',
               ],
      entry_points={
          'console_scripts': [
              'landlab = landlab.cmd:main',
              'landlab_ex_1 = landlab.examples.diffusion2Dtest:main',
              'landlab_ex_2 = landlab.examples.diffusion2Dtest2:main',
              'landlab_ex_3 = landlab.examples.diffusion2Dtest3:main',
          ]
      },
      test_suite='landlab.tests',
     )
