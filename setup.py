# -*- coding: utf-8 -*-
# @Author: Zebedee Nicholls
# @Date:   2017-04-10 13:42:11
# @Last Modified by:   Chris Smith
# @Last Modified time: 2018-01-11 19:17:00

from setuptools import setup
from setuptools import find_packages
import versioneer

# README #
def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='fair',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Python package to perform calculations with the FAIR simple climate model',
      long_description=readme(),
      keywords='simple climate model temperature response carbon cycle emissions forcing',
      url='https://github.com/OMS-NetZero/FAIR',
      author='OMS-NetZero, Chris Smith, Richard Millar, Zebedee Nicholls, Myles Allen',
      author_email='c.j.smith1@leeds.ac.uk, richard.millar@physics.ox.ac.uk',
      license='Apache 2.0',
      packages=find_packages(exclude=['tests*']),
      package_data={'': ['*.csv']},
      install_requires=[
          'numpy>=1.11.3',
          'scipy>=0.19.0',
      ],
      zip_safe=False,
)
