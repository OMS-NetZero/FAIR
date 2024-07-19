import fair
from fair._version import get_versions
import sphinx_rtd_theme

import sys
import os

html_context = {}  # correct?

# Set canonical URL from the Read the Docs Domain
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    html_context["READTHEDOCS"] = True

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'IPython.sphinxext.ipython_console_highlighting',
]

project = 'fair'
copyright = '2024, fair development team'
author = 'fair development team'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
version = get_versions()["version"]
release = version

# Be strict about any broken references
nitpicky = True
nitpick_ignore_regex = [
    (r'py:.*', 'np.ndarray'),
    (r'py:class', 'attr'),
]

# Number figures
numfig = True

# readthedocs being a pain
master_doc = 'index'
