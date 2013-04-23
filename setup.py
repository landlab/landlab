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
      packages=['landlab', 'landlab.plugins', 'landlab.tests'],
      entry_points={
          'console_scripts': ['landlab = landlab.cmd:main', ]
      },
      test_suite='landlab.tests',
     )
