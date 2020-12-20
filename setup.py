from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand

import versioneer

PACKAGE_NAME = "fair"
DESCRIPTION = (
    "Python package to perform calculations with the FaIR simple climate model"
)
KEYWORDS = [
    "simple climate model",
    "temperature response",
    "carbon cycle",
    "emissions",
    "forcing",
]

AUTHORS = [
    ("John Broadbent", "johngeoffreybroadbent@gmail.com"),
    ("Nicholas Leach", "nicholas.leach@stx.ox.ac.uk"),
    ("Chris Smith", "c.j.smith1@leeds.ac.uk"),
    ("Zeb Nicholls", "zebedee.nicholls@climate-energy-college.org"),
]

URL = "https://github.com/OMS-NetZero/FAIR"
PROJECT_URLS = {
    "Bug Reports": "https://github.com/OMS-NetZero/FAIR/issues",
    "Source": "https://github.com/OMS-NetZero/FAIR",
}
LICENSE = "Apache 2.0"
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: Apache License",
    "Intended Audience :: Developers",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
]

PYTHON_REQUIREMENTS = ">=3.6, <4"
REQUIREMENTS = [
    "matplotlib",
    "numexpr",
    "numpy",
    "pandas",
    "scipy",
    "pyam-iamc",
]
REQUIREMENTS_TESTS = ["pytest>=4.0", "nbval", "pytest-cov", "codecov"]
REQUIREMENTS_DOCS = ["sphinx>2.1", "sphinx_rtd_theme"]
REQUIREMENTS_DEPLOY = ["twine>=1.11.0", "setuptools>=41.2", "wheel>=0.31.0"]
REQUIREMENTS_NOTEBOOKS = ["notebook"]

REQUIREMENTS_DEV = [
    *["black==19.10b0", "flake8", "isort>=5",],
    *REQUIREMENTS_TESTS,
    *REQUIREMENTS_DEPLOY,
    *REQUIREMENTS_DOCS,
    *REQUIREMENTS_NOTEBOOKS,
]

REQUIREMENTS_EXTRAS = {
    "tests": REQUIREMENTS_TESTS,
    "deploy": REQUIREMENTS_DEPLOY,
    "dev": REQUIREMENTS_DEV,
    "docs": REQUIREMENTS_DOCS,
    "notebooks": REQUIREMENTS_NOTEBOOKS,
}

SOURCE_DIR = "fair"

PACKAGES = find_packages(exclude=["tests*", "docs*"])
PACKAGE_DATA = {"": ["*.csv"]}

README = "README.rst"

# Get the long description from the README file
with open(README, "r") as f:
    README_LINES = ["FaIR", "====", ""]
    add_line = False
    for line in f:
        if line.strip() == ".. sec-begin-long-description":
            add_line = True
        elif line.strip() == ".. sec-end-long-description":
            break
        elif add_line:
            README_LINES.append(line.strip())

if len(README_LINES) < 3:
    raise RuntimeError("Insufficient description given")


class Fair(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest

        pytest.main(self.test_args)


cmdclass = versioneer.get_cmdclass()
cmdclass.update({"test": Fair})

setup(
    name=PACKAGE_NAME,
    version=versioneer.get_version(),
    cmdclass=cmdclass,
    description=DESCRIPTION,
    long_description="\n".join(README_LINES),
    long_description_content_type="text/x-rst",
    author=", ".join([author[0] for author in AUTHORS]),
    author_email=", ".join([author[1] for author in AUTHORS]),
    url=URL,
    project_urls=PROJECT_URLS,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords=KEYWORDS,
    packages=PACKAGES,
    package_data=PACKAGE_DATA,
    include_package_data=True,
    install_requires=REQUIREMENTS,
    extras_require=REQUIREMENTS_EXTRAS,
    python_requires=PYTHON_REQUIREMENTS,
    zip_safe=False,
)
