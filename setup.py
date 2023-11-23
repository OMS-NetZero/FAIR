import pathlib

from setuptools import find_packages, setup

import versioneer

AUTHORS = [
    ("Chris Smith", "c.j.smith1@leeds.ac.uk"),
    ("Nicholas Leach", "nicholas.leach@stx.ox.ac.uk"),
    ("Stuart Jenkins", "stuart.jenkins@wadham.ox.ac.uk"),
    ("Richard Millar", "richard.millar@ouce.ox.ac.uk"),
    ("Zeb Nicholls", "zebedee.nicholls@climate-energy-college.org"),
    ("Myles Allen", "myles.allen@ouce.ox.ac.uk"),
]

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# using climate-assessment as a template here
REQUIREMENTS_INSTALL = [
    "matplotlib",  # not a requirement of fair, but needed for binder examples
    "numpy",
    "pandas",
    "pooch",
    "scipy",
    "tqdm",
    "xarray",
]
REQUIREMENTS_NOTEBOOKS = ["nbstripout", "jupyter", "ipywidgets", "ipython"]
REQUIREMENTS_TESTS = [
    "codecov",
    "nbmake",
    "netCDF4",
    "pytest-cov",
    "pytest-console-scripts",
    "pytest",
]
REQUIREMENTS_DOCS = ["ipython", "pandoc", "sphinx==6.2.1", "sphinx_rtd_theme==1.2.0"]
REQUIREMENTS_DEPLOY = [
    "build",
    "twine",
    "setuptools",
    "wheel",
]  # plus conda
REQUIREMENTS_STYLE = [
    "bandit",
    "black",
    "flake8",
    "isort",
    "pydocstyle",
]

requirements_dev = [
    *REQUIREMENTS_DOCS,
    *REQUIREMENTS_NOTEBOOKS,
    *REQUIREMENTS_TESTS,
    *REQUIREMENTS_DEPLOY,
    *REQUIREMENTS_STYLE,
]

requirements_dev_nodocs = [
    *REQUIREMENTS_NOTEBOOKS,
    *REQUIREMENTS_TESTS,
    *REQUIREMENTS_DEPLOY,
    *REQUIREMENTS_STYLE,
]

requirements_extras = {
    "docs": REQUIREMENTS_DOCS,
    "tests": REQUIREMENTS_TESTS,
    "deploy": REQUIREMENTS_DEPLOY,
    "dev": requirements_dev,
    "dev-nodocs": requirements_dev_nodocs,
}

setup(
    name="fair",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Finite-amplitude Impulse Response (FaIR) simple climate model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/OMS-NetZero/FAIR",
    author=", ".join([author[0] for author in AUTHORS]),
    author_email=", ".join([author[1] for author in AUTHORS]),
    license="Apache 2.0",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords="simple, climate, model, temperature, forcing, emissions, emulator",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={"": ["*.csv"]},
    python_requires=">=3.8, <4",
    install_requires=REQUIREMENTS_INSTALL,
    extras_require=requirements_extras,
)
