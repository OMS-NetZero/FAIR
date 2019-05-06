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
      description='Python package to perform calculations with the FaIR simple climate model',
      long_description=readme(),
      keywords='simple climate model temperature response carbon cycle emissions forcing',
      url='https://github.com/OMS-NetZero/FAIR',
      author='OMS-NetZero, Chris Smith, Richard Millar, Zebedee Nicholls, Myles Allen',
      author_email='c.j.smith1@leeds.ac.uk, richard.millar@physics.ox.ac.uk',
      license='Apache 2.0',
      packages=find_packages(exclude=['tests*','docs*']),
      package_data={'': ['*.csv']},
      include_package_data=True,
      install_requires=[
          'matplotlib',
          'numpy>=1.11.3',
          'scipy>=0.19.0',
      ],
      zip_safe=False,
      extras_require={'docs': ['sphinx>=1.4', 'nbsphinx'],
                      'dev' : ['notebook', 'wheel', 'twine'],
                      'test': ['pytest>=4.0', 'nbval', 'pytest-cov', 'codecov']}
)
