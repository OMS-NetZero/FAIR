import fair
from fair._version import get_versions
import sphinx_rtd_theme

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'IPython.sphinxext.ipython_console_highlighting',
]

project = 'FaIR'
copyright = '2023, FaIR Development Team'
author = 'FaIR Development Team'

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
